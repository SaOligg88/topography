#!/bin/python

# Initialize dataset
from neurosynth import Dataset
import pandas as pd
import numpy as np

# Load the Dataset and add topic-based features--point to right location
# in your local environment.
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

# Load binned data


# Run analysis


# Write out results


