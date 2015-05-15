function[dataAll] = ReassembleDistMat(sub, hemi)

% fix for right and left...
% r: 9354
% l: 9361
if hemi == 'lh'
    meshlen = 9354;
elseif hemi == 'rh'
    meshlen = 9361;
end
len = meshlen; 

data = zeros(len, meshlen);
for i = 1:len
    vec = load(['/scr/liberia1/data/lsd/surface/temp/' sub '_' hemi ...
        num2str(i-1) '.txt']);
    data(i,:) = [zeros(meshlen-length(vec),1); vec];
    % disp(i)
end

data = data + data';

% insert incld back into data
incld_fsa = load(['/scr/liberia1/data/lsd/surface/mesh/fsaverage5_' hemi '_incld.mat']); incld_fsa = incld_fsa.incld;
dataAll = zeros(10242);
dataAll(incld_fsa,incld_fsa) = data;

save(['/afs/cbs.mpg.de/projects/mar005_lsd-lemon-surf/probands/' sub '/distance_maps/' sub '_' hemi '_geoDist_fsa5.mat'],'dataAll','-v7.3')


