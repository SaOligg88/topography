# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

#!/usr/bin/python

import sys, getopt, scipy, os, h5py 
import numpy as np, pylab as pl
from sklearn.utils.arpack import eigsh  
from sklearn.cluster import KMeans
# import matplotlib
# from mayavi.mlab import *
# import nibabel.gifti.giftiio as gio
# from IPython.core.display import Image as im

def main(argv):
    
    # Set defaults:
    n_components_embedding = 25
    comp_min = 2
    comp_max = 20
    
    try:
        opts, args = getopt.getopt(argv,"hi:o:e:c",["subject=","filename="])
    except getopt.GetoptError:
        print 'embedding.py -i <input matrix> -o <output filebasename>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'embedding.py -i <input matrix> -o <output filebasename>'
            sys.exit()
        elif opt in ("-o", "--output"):
            filename = arg
        elif opt in ("-i", "--input"):
            sub = arg
        elif opt in ("-e", "--embedding_components"):
            n_components_embedding = arg
        elif opt in ("-c", "--components"):
            comp_min = arg[0]
            comp_max = arg[1] + 1
            
    # Import files
    f = h5py.File(('%s' % sub),'r')
    dataCorr = np.array(f.get('connData'))
    cortex = np.array(f.get('cortex')) - 1

    # fix cortex:
    cort = [0] * len(cortex)
    count = 0
    for i in cortex:
        cort[count] = int(i)
        count = count + 1
    
    # Run embedding and kmeans
    K = (dataCorr + 1) / 2.  
    v = np.sqrt(np.sum(K, axis=1)) 
    A = K/(v[:, None] * v[None, :])  
    del K
    A = np.squeeze(A * [A > 0])

    embedding = runEmbed(A, n_components_embedding)

    for n_components in xrange(comp_min,comp_max):   
        if n_components == 2:
            results = recort(np.squeeze(kmeans(embedding, n_components)), np.squeeze(np.sort(cort)))
        else:
            results = np.vstack((results, recort(np.squeeze(kmeans(embedding, n_components)), np.squeeze(np.sort(cort)))))
    
    scipy.io.savemat(('%s.mat' % filename), {'results':results})

def runEmbed(data, n_components):
    lambdas, vectors = eigsh(data, k=n_components)   
    lambdas = lambdas[::-1]  
    vectors = vectors[:, ::-1]  
    psi = vectors/vectors[:, 0][:, None]  
    lambdas = lambdas[1:] / (1 - lambdas[1:])  
    embedding = psi[:, 1:(n_components + 1)] * lambdas[:n_components][None, :]  
    #embedding_sorted = np.argsort(embedding[:], axis=1)
    return embedding

def kmeans(embedding, n_components):
    est = KMeans(n_clusters=n_components, n_jobs=-1, init='k-means++', n_init=300)
    est.fit_transform(embedding)
    labels = est.labels_
    data = labels.astype(np.float)
    return data
    
def recort(data, cortex):
    d = ([0] * np.shape(data)) 
    count = 0
    for i in cortex:
        d[i] = data[count] + 1
        count = count + 1
    return d

if __name__ == "__main__":
    main(sys.argv[1:])
    

