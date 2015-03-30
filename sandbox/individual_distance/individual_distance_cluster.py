#!/usr/bin/python

import sys, scipy, os, h5py, functions
import numpy as np, pylab as pl, nibabel as nib
from sklearn.cluster import KMeans

# Set defaults:
n_components_embedding = 25
comp_min = 2
comp_max = 20
subjects = []

# load cortex:
cort = ???

for sub in subjects:     
    dataAll = []  
    for run in [1,2,3,4]: # Import four files and norm
        data = nib.nifti.load()
        dataNorm = (data - data.mean(axis=0)) / data.stdev(axis=0)
        dataAll.hstack(dataNorm)
        
    # Correlate 
    dataCorr = np.corr(dataAll)

    # remove NaNs
    dataCorr.isnan() = 0
    
    # Embed and cluster
    embedding = runEmbed(prepMat(dataCorr), n_components_embedding)
    results = []
    for n_components in xrange(comp_min,comp_max):   
        results.append(recort(KMeans(embedding, n_components), cort))

    scipy.io.savemat(('%s.mat' % filename), {'results':results})

