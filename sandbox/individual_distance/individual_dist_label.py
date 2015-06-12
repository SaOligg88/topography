#!/usr/bin/python

import os, numpy as np, scipy as sp
import nibabel as nib
from surfer import Brain
import h5py
from sklearn import preprocessing
from sklearn.utils.arpack import eigsh  

# Set defaults:
dataDir = '/afs/cbs.mpg.de/projects/mar005_lsd-lemon-surf/probands'
fsDir = '/afs/cbs.mpg.de/projects/mar004_lsd-lemon-preproc/freesurfer'
output_base_dir = '/scr/liberia1'
subject = ['26410']


def vizBrain(data, subject_id='fsaverage5', hemi='lh', surface='pial', filename='brain.png'):
    brain = Brain(subject_id, hemi, surface)
    dmin = data.min()#+(data.std()/2)
    dmax = data.max()#-(data.std()/2)
    brain.add_data(data, dmin, dmax, colormap="hot", alpha=0.7)
    brain.save_montage(filename, order=['lat', 'med'], orientation='h', border_size=10)

def DoFiedler(conn):
    # prep for embedding
    K = (conn + 1) / 2.  
    v = np.sqrt(np.sum(K, axis=1)) 
    A = K/(v[:, None] * v[None, :])  
    del K
    A = np.squeeze(A * [A > 0])

    # diffusion embedding
    n_components_embedding = 2
    lambdas, vectors = eigsh(A, k=n_components_embedding+1)  
    del A
    lambdas = lambdas[::-1]  
    vectors = vectors[:, ::-1]  
    psi = vectors/vectors[:, 0][:, None]  
    lambdas = lambdas[1:] / (1 - lambdas[1:])  
    embedding = psi[:, 1:(n_components_embedding + 1 + 1)] * lambdas[:n_components_embedding+1][None, :]  

    return embedding



    
for subject in subjects:
    for hemi in ['lh', 'rh']:
        
        # read in data
        cort = np.sort(fs.io.read_label('%s/fsaverage5/label/%s.cortex.label' % (fsDir, hemi)))
        dataCorr = np.load('%s/%s/correlation_maps/%s_lsd_corr_1ab_fsa5_%s.npy' % (dataDir, subject, subject, hemi))
        fullsize = len(dataCorr)
        
        distFile = '%s/%s/distance_maps/%s_%s_geoDist_fsa5.mat' % (dataDir, subject, subject, hemi)
        with h5py.File(distFile, 'r') as f:
                dist = f['dataAll'][()]
        dist = dist[cort, :][:, cort]
        min_max_scaler = preprocessing.MinMaxScaler()
        dist_scaled = min_max_scaler.fit_transform(dist)
        del dist        
        distmat = np.zeros((fullsize, fullsize))
        distmat[np.ix_(cort, cort)] = dist_scaled
        del dist_scaled
        
        
        # embedding
        embedding = DoFiedler(dataCorr[cort, :][:, cort]) # see below for details
        del dataCorr
        # reinsert zeros:
        fiedler = np.zeros(fullsize)
        fiedler[cort] = embedding[:,1] # before this was embedding[:,0], 
        # but changed to 1 since fiedler vector is the vector belonging second smallest eigenvalue
       
         
        # TODO: use individual-specific annot file, need to be transformed from individual space to fsa5 space
        fs_annot = fs.io.read_annot('/afs/cbs.mpg.de/projects/mar004_lsd-lemon-preproc/freesurfer/fsaverage5/label/lh.aparc.a2009s.annot')#'%s/%s/label/%s.aparc.a2009s.annot' % (fsDir, subject, hemi))
        index = [i for i, s in enumerate(list(fs_annot[2])) if 'G_pariet_inf-Angular' in s]
        label_parietal =  fs_annot[0] == index
        masked_fiedler = fiedler * label_parietal
        parietal_index = np.where(masked_fiedler == max(masked_fiedler))
        # changed np.mean(fiedler) to np.mean(np.nonzero(fiedler)), is that correct?
        if np.mean(np.nonzero(masked_fiedler)) > np.mean(np.nonzero(fiedler)):
            parietal_index = np.where(masked_fiedler == max(masked_fiedler))
        else:
            parietal_index = np.where(masked_fiedler == min(masked_fiedler))  
        
        # TODO: add this as second overlay to surface instead of changing value in masked_fiedler
        masked_fiedler[parietal_index] = 1
        vizBrain(masked_fiedler) # TODO: save to disc for qc
        
        
        label_list = ['S_calcarine'] # TODO: add other labels
        g = 0
        indices = [i for i, s in enumerate(list(fs_annot[2])) if label_list[g] in s]
        label_dist = np.min(distmat[np.where(fs_annot[0] == indices),:], axis=1).squeeze()
        # TODO: add this as second overlay to surface instead of changing value in label_dist
        label_dist[parietal_index] = 1
        vizBrain(label_dist) # TODO: save to disc
                
        
        
        
        
        
        
        
        
        # save out anat_dist for subject / hemi / anat label
        # also create images for quality control: fiedler, masked_fiedler



