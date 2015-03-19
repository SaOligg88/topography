#!/usr/bin/python

import sys, os, h5py, scipy, numpy as np
from sklearn.utils.arpack import eigsh
from sklearn.cluster import KMeans
from scipy.io.matlab import savemat

def main(argv):
    
    # Set defaults:
    n_components_embedding = 25
    comp_min = 2
    comp_max = 20 + 1
    varname = 'data'
    filename = './test'
    
    # Import files
    f = h5py.File(('%s.mat' % filename),'r')
    dataCorr = np.array(f.get('%s' % varname))

    # Prep matrix
    K = (dataCorr + 1) / 2.  
    v = np.sqrt(np.sum(K, axis=1)) 
    A = K/(v[:, None] * v[None, :])  
    del K
    A = np.squeeze(A * [A > 0])

    # Run embedding
    lambdas, vectors = eigsh(A, k=n_components_embedding)   
    lambdas = lambdas[::-1]  
    vectors = vectors[:, ::-1]  
    psi = vectors/vectors[:, 0][:, None]  
    lambdas = lambdas[1:] / (1 - lambdas[1:])  
    embedding = psi[:, 1:(n_components_embedding + 1)] * lambdas[:n_components_embedding][None, :]

    # Run kmeans clustering

    def kmeans(embedding, n_components):
        est = KMeans(n_clusters=n_components, n_jobs=-1, init='k-means++', n_init=300)
        est.fit_transform(embedding)
        labels = est.labels_
        data = labels.astype(np.float)
        return data

    results = list()
    for n_components in xrange(comp_min,comp_max):   
        results.append(kmeans(embedding, n_components))

    savemat(('%s_results.mat' % filename), {'results':results})

if __name__ == "__main__":
    main(sys.argv[1:])

