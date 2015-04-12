function[dataAll] = ReassembleDistMat(outputPrefix)

meshlen = 29931;
len = meshlen; 

data = zeros(len, meshlen);
for i = 1:len
    vec = load(['/scr/litauen2/projects/distance/condor/' ...
        num2str(i-1) '.txt']);
    data(i,:) = [zeros(meshlen-length(vec),1); vec];
    disp(i);
end

data = data + data';

% insert incld back into data
load incld.mat
dataAll = zeros(32492);
dataAll(incld,incld) = data;
