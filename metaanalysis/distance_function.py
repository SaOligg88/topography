#!/bin/python

# Initialize dataset
from neurosynth import Dataset
import pandas as pd
import numpy as np
import sys
sys.path.insert(0,'/scr/litauen1/toro/dist_meta/neurovault/neuro2/neurosynth')
import ply
import numpy as np
import neurosynth as ns
from neurosynth.base.dataset import Dataset
from neurosynth.analysis import decode

# Load the Dataset and add topic-based features--point to right location
# in your local environment.

# dataset = Dataset('database.txt')
# dataset.add_features('features.txt')
# dataset.save("dataset.pkl", keep_mappables = True)
pickled_dataset = '/scr/litauen1/toro/dist_meta/neurovault/neurosynth/dataset.pkl'
dataset = Dataset.load(pickled_dataset)

features = pd.read_csv('v3-topics-50.txt', sep='\t', index_col=0)
# Remove junk topics and replace with more sensible names. Check the topic-keys file
# to see what terms each topic loads on. These topics are identical to those at
# http://neurosynth.org/analyses/topics/v3-topics-50/
topics_to_keep = [1, 4, 6, 14, 16, 17, 18, 20, 21, 22, 23, 25, 27, 29, 30, 31,
                  33, 35, 36, 37, 38, 41, 44, 45, 46, 48, 49] # 43
labels = ['faces', 'comprehension', 'cues', 'WM', 'autonomic', 'reading', 'autobiographical',
          'motion', 'numerical', 'olfactory', 'inhibition', 'motor', 'reward', 'attention',
          'polymodal', 'object', 'eye_movement', 'action', 'auditory', 'lexical', 'pain',
          'memory', 'categorical', 'emotion', 'learning', 'control', 'social'] # 'stress'
features = features.iloc[:, topics_to_keep]
features.columns = labels
dataset.add_features(features, append=False)
decoder = decode.Decoder(dataset, method = 'roi')

# Load binned data
data = decoder.decode(['/scr/litauen1/toro/dist_meta/distDMN_MNI_2mm_0_5.nii.gz',
			'/scr/litauen1/toro/dist_meta/distDMN_MNI_2mm_5_10.nii.gz',
			'/scr/litauen1/toro/dist_meta/distDMN_MNI_2mm_10_15.nii.gz', 
                       '/scr/litauen1/toro/dist_meta/distDMN_MNI_2mm_15_20.nii.gz',
			'/scr/litauen1/toro/dist_meta/distDMN_MNI_2mm_20_25.nii.gz',
			'/scr/litauen1/toro/dist_meta/distDMN_MNI_2mm_25_30.nii.gz',
                       '/scr/litauen1/toro/dist_meta/distDMN_MNI_2mm_30_35.nii.gz', 
                       '/scr/litauen1/toro/dist_meta/distDMN_MNI_2mm_35_40.nii.gz',
                       '/scr/litauen1/toro/dist_meta/distDMN_MNI_2mm_40_45.nii.gz', 
                       '/scr/litauen1/toro/dist_meta/distDMN_MNI_2mm_45_50.nii.gz', 
                       '/scr/litauen1/toro/dist_meta/distDMN_MNI_2mm_50_55.nii.gz', 
                       '/scr/litauen1/toro/dist_meta/distDMN_MNI_2mm_55_60.nii.gz', 
                       '/scr/litauen1/toro/dist_meta/distDMN_MNI_2mm_60_65.nii.gz',
                       '/scr/litauen1/toro/dist_meta/distDMN_MNI_2mm_65_70.nii.gz', 
                       '/scr/litauen1/toro/dist_meta/distDMN_MNI_2mm_70_75.nii.gz',
                       '/scr/litauen1/toro/dist_meta/distDMN_MNI_2mm_75_80.nii.gz', 
                       '/scr/litauen1/toro/dist_meta/distDMN_MNI_2mm_80_85.nii.gz',
			'/scr/litauen1/toro/dist_meta/distDMN_MNI_2mm_85_90.nii.gz', 
                       ], save='decoding_results.txt')

# Write out results
import csv
lol = list(csv.reader(open('decoding_results.txt', 'rb'), delimiter='\t'))
lola = np.array(lol)

inputAll = []
inputNorm = []
namesAll = []
names = lola[1:,0]
size = np.shape(lola)
for i in xrange(1,size[1]):
    input = lola[1:,i]
    ord = np.argsort(input)
    inputAll.append(input[ord][:])
    namesAll.append(names[ord][:])
    inputNorm.append(input)
    
print np.array(namesAll)[:,-5:]

print np.array(inputAll)[:,-6:]

# Save files for matlab
import h5py
f = h5py.File('../../neurosynth_meta.mat','w')
f['inputAll'] = np.asarray(inputAll)
f['namesAll'] = namesAll
f['names'] = names
f['inputNorm'] = inputNorm
f.close()


