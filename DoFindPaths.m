function[clus, randpath] = DoPaths(n)
% n = number of permutations
% set to 0 to not run any purmutations
%

addpath(genpath('./utils'));
%% Load yeo preconfigured:
[thresh, surf, surfi, surfm, label, edg, lab] = config_yeo_L_17rsns();

%%%%%%%%% Begin %%%%%%%%%

clus = pathsFindHCP(label, thresh, hemi);

% save('data/clus.mat', '-v7.3', 'clus');

if n > 0
	randPath = pathsPermute(clus, n);
else
	randPath = [];
end
