function[cortex, noncortex] = loadCortex(hemi, surf)

sub = [100307];

aparc = gifti(['/scr/dattel2/' num2str(sub) '/MNINonLinear/fsaverage_LR32k/' num2str(sub) '.L.aparc.a2009s.32k_fs_LR.label.gii']);
aparc = aparc.cdata;

% Make noncortex mask
noncortex = zeros([1 length(surf.coord)]);
noncortex(find(aparc == 0)) = 1;
cortex = zeros([1 length(surf.coord)]);
cortex(find(noncortex == 0)) = 1;