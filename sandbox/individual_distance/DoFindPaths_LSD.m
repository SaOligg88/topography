function[dist2mask] = DoFindPaths_LSD(subID)
% subID = subjectID in quotes

% set variables:
addpath(genpath('../../utils'));
thresh = 2;		
hemi   = 'lh'; 	% 'lh' or 'rh'
dir1   = ''; 	% locations of freesurfer directory containing subID

% Load surface:
filename = ['' hemi ''];
surf = SurfStatRead(filename);
% Load cluster results: (output from individual_distance_cluster.py)
filename = [];
label = load();
% Decide on which number cluster solution:
label = label.results(:,XXX);

%%%%%%%%% Begin %%%%%%%%%
clus = pathsFindHCP(label, thresh, hemi, surf);

%% Plot outputs of clus:
[matrix, s, ind, x, y] = threshGraph(clus);

% Transform fsaverage space to individual space
save clus solution
!mri_surf2surf # check ft package for transform...
load clus_individual_space

% Get cluster of interest in parietal from DMN

% Get distances from DMN:
source = find(clus_individual_space == clus_#); 
[dist, ~] = distExactGeodesic(source, 'freesurfer', hemi, 'distance', [dir1 subID]);

% Get min dist value to masks from primary (using aparc)
mask = [] % label values from aparc
for i = 1:length(mask)
	dist2mask(i) = min(dist(find(aparc == mask(i))));
end
