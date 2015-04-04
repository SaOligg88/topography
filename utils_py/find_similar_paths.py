#get_ipython().magic(u'matplotlib inline')
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
import scipy.io as sio
import h5py
import networkx as nx
import numpy as np
import gdist

def calcPaths(num):
    length = nx.all_pairs_dijkstra_path(G, num)
    length_paths = []
    for node in length:
        for target in length[node]:
            if len(length[node][target]) == num:
                length_paths.append(length[node][target])
    labeled_paths = labels[length_paths]
    same_labels = (squareform(pdist(labeled_paths)) < 1e-10).sum(axis=1)
    return length_paths, labeled_paths, same_labels

def uniqueRows(labeled_paths, same_labels, cutoff):
    a = labeled_paths[same_labels == cutoff]
    uRows = np.unique(a.view(np.dtype((np.void, a.dtype.itemsize*a.shape[1])))).view(a.dtype).reshape(-1, a.shape[1])
    return uRows
    print uRows
    
def removePaths(labeled_paths, same_labels, vals, row):
    ind = np.in1d(labeled_paths[:,row], vals).reshape(labeled_paths[:,row].shape)
    labeled_paths_new = labeled_paths[ind]
    same_labels_new = same_labels[ind]
    return ind, labeled_paths_new, same_labels_new

def removePathsInverse(labeled_paths, same_labels, vals, row):
    ind = np.in1d(labeled_paths[:,row], vals, invert=True).reshape(labeled_paths[:,row].shape)
    labeled_paths_new = labeled_paths[ind]
    same_labels_new = same_labels[ind]
    return ind, labeled_paths_new, same_labels_new

def printAll(labeled_paths,same_labels):
    print (('1 times:\n %s \n' % uniqueRows(labeled_paths,same_labels,1)))
    print (('2 times:\n %s \n' % uniqueRows(labeled_paths,same_labels,2)))
    print (('3 times:\n %s \n' % uniqueRows(labeled_paths,same_labels,3)))
    print (('4 times:\n %s \n' % uniqueRows(labeled_paths,same_labels,4)))
    print (('5 times:\n %s \n' % uniqueRows(labeled_paths,same_labels,5)))
    print (('6 times:\n %s \n' % uniqueRows(labeled_paths,same_labels,6)))

def labelHist(num,labeled_paths):
    a = np.zeros([num,18])
    for i in xrange(0,num):
        a[i] = np.histogram(labeled_paths[:,i],18, range=(1,18))[0]
    return a.transpose()

def calcPlace(b):
    c = np.zeros([len(b),1])
    for i in xrange(0,len(b)):
        if np.sum(b[i]) == 0:
            c[i] = 0
        else:
            c[i] = np.average(np.array(xrange(0,num)) + 1, weights=b[i])
    return c


'''load data:'''
fp = h5py.File('data/clus.mat')
fp.keys()
adj = fp['clus']['edge'][:]
labels = fp['clus']['edgeNet'][:].flatten()
G = nx.from_numpy_matrix(adj)

num=6
length_paths, labeled_paths, same_labels = calcPaths(num)
ind1, labeled_paths1, same_labels1 = removePaths(labeled_paths, same_labels, [7,12,14,15,16], 0)
ind2, labeled_paths2, same_labels2 = removePathsInverse(labeled_paths1, same_labels1, [7,12,14,15,16], 1)
ind3, labeled_paths3, same_labels3 = removePathsInverse(labeled_paths2, same_labels2, [7,12,14,15,16], 2)
ind4, labeled_paths4, same_labels4 = removePaths(labeled_paths3, same_labels3, [17,3,8], num-2)
ind5, labeled_paths5, same_labels5 = removePaths(labeled_paths4, same_labels4, [13], num-1)
print '\nNum = %s' % num
printAll(labeled_paths1, same_labels1)

'''histograms'''
indX, labeled_pathsX, same_labelsX = removePathsInverse(labeled_paths3, same_labels3, [7,12,14,15,16], 3)

num=5
length_paths, labeled_paths, same_labels = calcPaths(num)
ind1, labeled_paths1, same_labels1 =     removePaths(labeled_paths, same_labels, [7,12,14,15,16], 0)
a = labelHist(num,labeled_paths1)
print a
order = np.argsort(calcPlace(a), axis=0) + 1
fl = np.floor(np.sort(calcPlace(a), axis=0))
print '\nNum = %s' % num
print np.hstack((order, fl))

num=6
length_paths, labeled_paths, same_labels = calcPaths(num)
ind1, labeled_paths1, same_labels1 =     removePaths(labeled_paths, same_labels, [7,12,14,15,16], 0)
a = labelHist(num,labeled_paths1)
order = np.argsort(calcPlace(a), axis=0) + 1
fl = np.floor(np.sort(calcPlace(a), axis=0))
print '\nNum = %s' % num
print np.hstack((order, fl))

np.array(fl[np.argsort(order, axis=0)]).transpose()

np.array(length_paths)[same_labels == same_labels.max()]
