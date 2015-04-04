function[clus] = DoFindPaths()
% n = number of permutations
% set to 0 to not run any purmutations

addpath(genpath('./utils'));
addpath(genpath('./config'));

%% Load yeo preconfigured:
[thresh, surf, surfi, surfm, label, edg, lab] = config_yeo_L_17rsns();

%%%%%%%%% Begin %%%%%%%%%
results = load('/scr/litauen1/Dropbox/misc/embedding/embedding_hcp_2_30.mat');
label = results.results(17,:); % for 18 cluster solution
thresh = 2;
hemi = 'L';

clus = pathsFindHCP(label, thresh, hemi, surf);
save('data/clus.mat', '-v7.3', 'clus');

%% Plot outputs of clus:
numClus = length(unique(nonzeros(clus.regions)));
[matrix, s, ind, x, y] = threshGraph(clus, surf, ['results/label' num2str(numClus)]);

%% Save: 
%colors = [0 0 0; cbrewer('qual','Set1',numClus)];
colors = [0 0 0; csvread('results/colormaps/region_colors.csv')];
WriteSurfMap(clus.regions, ['results/label' num2str(numClus)], colors, surfm);

%% clustering
k = 8;
clusters = clustering_kmeans(k, x, y, surf);

%% permute: 
% Todo: set steps by order of ind:
steps.s{1} = [find(clusters == 4); find(clusters == 5); find(clusters == 8)];
steps.s{2} = find(clusters == 1);
steps.s{3} = find(clusters == 3);
steps.s{4} = find(clusters == 7);
steps.s{5} = find(clusters == 6);
steps.s{6} = find(clusters == 2);
networks = [1:length(unique(nonzeros(clus.label)))];

[randPath] = pathsPermute(clus, 100, steps, networks);
