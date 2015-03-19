#!/usr/bin/python

import sys, getopt, scipy, os, h5py 
import numpy as np, pylab as pl
from sklearn.utils.arpack import eigsh  
from sklearn.cluster import KMeans

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
            print 'embedding.py -i <input matrix> <variable name> -o <output filebasename> \n 
            Example: \n
                python embedding.py -i ./test.mat corr -o ./test_output -e 25 -c 2 20'
            sys.exit()
        elif opt in ("-o", "--output"):
            filename = arg
        elif opt in ("-i", "--input"):
            sub = arg[0]
            varname = arg[1]
        elif opt in ("-e", "--embedding_components"):
            n_components_embedding = arg
        elif opt in ("-c", "--components"):
            comp_min = arg[0]
            comp_max = arg[1] + 1
            
    # Import files
    f = h5py.File(('%s' % sub),'r')
    dataCorr = np.array(f.get('%s' % varname))

    # Run embedding and kmeans
    K = (dataCorr + 1) / 2.  
    v = np.sqrt(np.sum(K, axis=1)) 
    A = K/(v[:, None] * v[None, :])  
    del K
    A = np.squeeze(A * [A > 0])

    lambdas, vectors = eigsh(A, k=n_components_embedding)   
    lambdas = lambdas[::-1]  
    vectors = vectors[:, ::-1]  
    psi = vectors/vectors[:, 0][:, None]  
    lambdas = lambdas[1:] / (1 - lambdas[1:])  
    embedding = psi[:, 1:(n_components_embedding + 1)] * lambdas[:n_components_embedding][None, :]

    for n_components in xrange(comp_min,comp_max):   
        if n_components == 2:
            results = np.squeeze(kmeans(embedding, n_components))
        else:
            results = np.vstack(results, np.squeeze(kmeans(embedding, n_components)))
    
    scipy.io.savemat(('%s.mat' % filename), {'results':results})

def kmeans(embedding, n_components):
    est = KMeans(n_clusters=n_components, n_jobs=-1, init='k-means++', n_init=300)
    est.fit_transform(embedding)
    labels = est.labels_
    data = labels.astype(np.float)
    return data

if __name__ == "__main__":
    main(sys.argv[1:])
