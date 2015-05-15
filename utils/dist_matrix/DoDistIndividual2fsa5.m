function [] = DoDistIndividual2fsa5(sub, hemi)

% set variables and paths
addpath(genpath('/scr/liberia1/workspace/topography/utils'));
addpath(pathdef);
%sub = num2str(sub);

% run for each individual:
% Run prep surface:
prepSurf('freesurfer', hemi, ['/scr/ilz2/LEMON_LSD/freesurfer/' sub], 'remove', ['/scr/liberia1/data/lsd/surface/mesh/' sub '_' hemi]);

surf_ind_s = SurfStatReadSurf(['/scr/ilz2/LEMON_LSD/freesurfer/' sub '/surf/' hemi '.sphere']);
surf_fsa5_s = SurfStatReadSurf(['/scr/ilz2/LEMON_LSD/freesurfer/fsaverage5/surf/' hemi '.sphere']);

incld_ind = load(['/scr/liberia1/data/lsd/surface/mesh/' sub '_' hemi '_incld.mat']); incld_ind = incld_ind.incld;
s_ind_coord = surf_ind_s.coord(:,incld_ind);

% include hemisphere in below:
incld_fsa = load(['/scr/liberia1/data/lsd/surface/mesh/fsaverage5_' hemi '_incld.mat']); incld_fsa = incld_fsa.incld;
s_fsa_coord = surf_fsa5_s.coord(:,incld_fsa);

% find matching vertex in individ from fsa5
clear minVal indVal indVal0
for i = 1:length(s_fsa_coord)
    % figure out way to compute using matrix operations, rather than
    % iteratively
    [minVal(i) indVal(i)] = min(sum(abs(bsxfun(@minus,s_fsa_coord(:,i),s_ind_coord)).^2).^0.5);
    % disp(i);
end
% subtract one from indices for c++
indVal0 = indVal-1;
fid = fopen(['/scr/liberia1/data/lsd/surface/mesh/' hemi '_' sub '_indices.txt'], 'w');
fprintf(fid,'%d\n',indVal0');
fclose(fid);

% For testing whether projection worked:
% surf_ind_i = SurfStatReadSurf(['/scr/ilz2/LEMON_LSD/freesurfer/' sub '/surf/' hemi '.inflated']);
% a = zeros(length(surf_ind_i.coord),1); a(indVal) = indVal;
% figure; SurfStatView(a,surf_ind_i);
% % save fig option...
% 
% surf_fsa5_i = SurfStatReadSurf('/scr/ilz2/LEMON_LSD/freesurfer/fsaverage5/surf/lh.inflated');
% b = zeros(length(surf_fsa5_i.coord),1); b = indVal;
% figure; SurfStatView(b,surf_fsa5_i);
% % save fig option...

%% Actually run distance calculation:
% notes on condor etc...
disp(['./exactGeodesicMatrix /scr/liberia1/data/lsd/surface/mesh/' sub '_' hemi '.asc /scr/liberia1/data/lsd/surface/temp/ ' sub '_' hemi  ' /scr/liberia1/data/lsd/surface/mesh/' hemi '_' sub '_indices.txt 7426']);



