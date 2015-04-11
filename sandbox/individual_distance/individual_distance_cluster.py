#!/usr/bin/python

import os, numpy as np, scipy as sp, nibabel.freesurfer as fs
from sklearn.cluster import KMeans
from sklearn.utils.arpack import eigsh  

# Set defaults:
n_components_embedding = 25
comp_min = 2
comp_max = 20
subjects = []
base_dir = '/scr/liberia1/LEMON_LSD/LSD_rest_surf'
output_base_dir = '/scr/liberia1'

subjects = [26410]
    
for subject in subjects:
    for hemi in ['lh', 'rh']:
        cort = np.sort(fs.io.read_label(('/scr/ilz2/LEMON_LSD/freesurfer/fsaverage5/label/%s.cortex.label' % hemi))) 
        dataNorm = list([])
        for session in ['rest1a', 'rest1b', 'rest2a', 'rest2b']:    
            # load data
            dataAll = np.array(fs.mghformat.MGHImage.load(
                            os.path.join(base_dir,'lsd_%s_%s_preprocessed_fsaverage5_%s.mgz' 
                            % (session, str(subject), hemi))
                            ).get_data().squeeze())
            # remove noncortex:
            data = dataAll[cort,:]
            # norm and concatenate
            dataNorm.extend((data.transpose() - data.mean(axis=1)) / data.std(axis=1))
            del data

        # get full ength of surface dataset for padding cluster results at end 
        fullsize = len(dataAll)
        del dataAll
        
        # correlate
        dataCorr = np.corrcoef(np.transpose(np.array(dataNorm)))
        del dataNorm
        dataCorr[np.isnan(dataCorr)] = 0
        # prep for embedding
        K = (dataCorr + 1) / 2.  
        del dataCorr
        v = np.sqrt(np.sum(K, axis=1)) 
        A = K/(v[:, None] * v[None, :])  
        del K
        A = np.squeeze(A * [A > 0])
        # diffusion embedding
        lambdas, vectors = eigsh(A, k=n_components_embedding+1)  
        del A
        lambdas = lambdas[::-1]  
        vectors = vectors[:, ::-1]  
        psi = vectors/vectors[:, 0][:, None]  
        lambdas = lambdas[1:] / (1 - lambdas[1:])  
        embedding = psi[:, 1:(n_components_embedding + 1 + 1)] * lambdas[:n_components_embedding+1][None, :]  
        # kmeans clustering
        results = []
        for n_components in xrange(comp_min, comp_max+1):
            est = KMeans(n_clusters=n_components, n_jobs=-1, init='k-means++', n_init=300)
            est.fit_transform(embedding)
            labels = est.labels_
            clust = labels.astype(np.float)
            # reinsert zeros:
            padded = np.zeros(fullsize)
            padded[cort] = clust + 1                      
            results.append(padded)
            
        # Save to matfile:
        sp.io.savemat((os.path.join(output_base_dir,'clusters_%s_%s_em%s_%s_%s.mat' % 
                (str(subject), hemi, str(n_components_embedding), str(comp_min), str(comp_max)))), 
                {'results':results})
                