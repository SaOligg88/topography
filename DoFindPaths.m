function[clus] = DoFindPaths()
% n = number of permutations
% set to 0 to not run any purmutations
%

addpath(genpath('./utils'));
%% Load yeo preconfigured:
[thresh, surf, surfi, surfm, label, edg, lab] = config_yeo_L_17rsns();

%%%%%%%%% Begin %%%%%%%%%
load ../yeoTopo/
surf_gii = gifti('data/Q1-Q6_R440.L.midthickness.32k_fs_LR.surf.gii');
surf.coord = surf_gii.vertices'; surf.tri = surf_gii.faces;
thresh = 2;
hemi = 'L';
clus = pathsFindHCP(label, thresh, hemi, surf);

% save('data/clus.mat', '-v7.3', 'clus');

%% Plot outputs of clus:
[matrix, s, ind, x, y] = threshGraph(clus);

%% clustering
k = 8;
clusters = clustering_kmeans(k, x, y, surf);

%% permute: 
steps.s{1} = [find(clusters == 4); find(clusters == 5); find(clusters == 8)];
steps.s{2} = find(clusters == 1);
steps.s{3} = find(clusters == 3);
steps.s{4} = find(clusters == 7);
steps.s{5} = find(clusters == 6);
steps.s{6} = find(clusters == 2);
networks = [1:length(unique(nonzeros(clus.label)))];

[randPath] = pathsPermute(clus, 100, steps, networks)
