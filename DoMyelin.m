function[] = DoMyelin()

addpath(genpath('./utils'));

%% variables: 
hemi = 'L'; % or: 'R'

%% Load data:
distDMN = loadDistDMN_group(hemi);
myelin = loadMyelin_group(hemi);
[surf, surfi, surfm] = loadHCPsurf_group(hemi);
[cortex, noncortex] = loadCortex(hemi);

h = figure;
plot(myelin(cortex), distDMN(cortex), '.');
