function [surf, data] = loadDistCorrData(sub)

file = ['/scr/kaiser2/corticalDist/' num2str(sub) '.L.dist.32k.mat'];
dist = load(file);
data.dist = dist.dist;
% import surface:
filename = ['/scr/dattel2/' num2str(sub) '/MNINonLinear/fsaverage_LR32k/' num2str(sub) ...
    '.L.midthickness.32k_fs_LR.surf.gii'];
surf_gii = gifti(filename);
surf.coord = surf_gii.vertices';
surf.tri = surf_gii.faces;
clear surf_gii;
filename = ['/scr/dattel2/' num2str(sub) '/MNINonLinear/fsaverage_LR32k/' num2str(sub) ...
    '.L.midthickness.32k_fs_LR.surf.gii'];
surf_gii = gifti(filename);
surf.coord = surf_gii.vertices';
surf.tri = surf_gii.faces;

filename = ['/scr/dattel2/' num2str(sub) '/MNINonLinear/fsaverage_LR32k/' num2str(sub) ...
    '.L.inflated.32k_fs_LR.surf.gii'];
surf_gii = gifti(filename); surfi.coord = surf_gii.vertices'; surfi.tri = surf_gii.faces;
filename = ['/scr/dattel2/' num2str(sub) '/MNINonLinear/fsaverage_LR32k/' num2str(sub) ...
    '.L.white.32k_fs_LR.surf.gii'];
surf_gii = gifti(filename); surfw.coord = surf_gii.vertices'; surfw.tri = surf_gii.faces;
clear surf_gii;

surfm.tri = surf.tri;
surfm.coord = (surfi.coord .* 0.7 + surfw.coord .* 0.3);

aparc = gifti(['/scr/dattel2/' num2str(sub) '/MNINonLinear/fsaverage_LR32k/' num2str(sub) ...
    '.L.aparc.a2009s.32k_fs_LR.label.gii']);
aparc = aparc.cdata;

% Make noncortex mask
noncortex = zeros([1 length(surf.coord)]);
noncortex(find(aparc == 0)) = 1;
cortex = zeros([1 length(surf.coord)]);
cortex(find(noncortex == 0)) = 1;
data.cort = find(cortex);

dim = 32492;
%filename = ['/scr/murg2/HCP_Q3_glyphsets_left-only/' num2str(sub) '/rfMRI_REST_left_corr_avg.gii.data'];
%file = fopen(filename,'r');
%data.conn = fread(file,[dim,dim],'float32');

% Group:
conn = load('/scr/murg2/HCP_new/HCP_Q1-Q6_GroupAvg_Related440_Unrelated100_v1/HCP_Q1-Q6_R468_rfMRI_groupAvg_left_corr_smoothed2_toR_nifti.mat');
data.conn = conn.dataCorr;
data.conn(isnan(data.conn)) = 0;
data.conn(find(noncortex), find(noncortex)) = 0;

[Yh,Ih] = sort(data.conn,1);
[Yh1,Ih1]=sort(Ih,1);
data.corrMath = zeros([32492 32492]);
data.corrMath(find(Ih1 > (32492-(32492*.1)))) = 1;


dataPcorr = zeros([1 32492]);
dataP = gifti('100307.L.geoDistPrimary.32k.shape.gii');
data.P = dataP.cdata;
dataR = gifti('100307.L.geoDist.32k.shape.gii');
data.R = dataR.cdata;





