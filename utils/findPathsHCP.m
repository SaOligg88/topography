function[clus] = findPathsHCP(label, thresh, hemi)
% 'label' is a vector of the length of the surface with integers denoting
% groups. 
% 'thresh' is equal to the minimum number of adjacent nodes that can
% consistiute a cluster. This was put in to remove single-node clusters. 
% I generally set it somewhere between 0 and 2.

addpath(genpath('./utils/surfstat'));

%hemi = 'L';
dir = ['/scr/murg2/HCP_new/HCP_Q1-Q6_GroupAvg_Related440_Unrelated100_v1/'];
surf_gii = gifti([dir 'Q1-Q6_R440.' hemi '.midthickness.32k_fs_LR.surf.gii']);
surfm.coord = surf_gii.vertices'; surfm.tri = surf_gii.faces;
surf_gii = gifti([dir 'Q1-Q6_R440.' hemi '.very_inflated.32k_fs_LR.surf.gii']);
surfi.coord = surf_gii.vertices'; surfi.tri = surf_gii.faces;
surf.tri = surfm.tri; surf.coord = (surfi.coord .* 0.8 + surfm.coord .* 0.2);

edg = SurfStatEdg(surf);
lab = nonzeros(unique(label));
clus.label = zeros([1 length(surf.coord)]);
countC = 0;
clus.label = zeros([1 32492]);
clus.network = zeros([1 32492]);
for i = 1:length(lab)    
    a = zeros([length(surf.coord) 1]);
    a(find(label == lab(i))) = 1; 
    slm = struct();
    slm.tri = surf.tri';
    slm.t = a';
    [cluster,clusid] = SurfStatPeakClus(slm,ones([length(surf.coord) 1]),0.05, ones(1,length(surf.coord)), edg);  
    clear slm
    for j = 1:length(clusid.clusid)
        nodes = cluster.vertid(find(cluster.clusid == j));
        if length(nodes) < thresh
            disp(['cut ' num2str(i) ' ' num2str(j) ' ' num2str(length(nodes))]);
        else        
            countC = countC + 1;
            clus.label(nodes) = countC;   
            clus.network(nodes) = lab(i);       
        end
    end
    clear cluster clusid
end

clus.edge = zeros(countC);
for i = 1:countC
    a = [edg(find(ismember(edg(:,1),find(clus.label == i))),2); ...
        edg(find(ismember(edg(:,2),find(clus.label == i))),1)];
    
    clus.edge(i,...
        nonzeros(unique(clus.label(a))))...
        = 1;
end
for i = 1:length(clus.edge)
    clus.edgeNet(i) = unique(clus.network(clus.label == i));
end

clus.netScore = zeros(length(lab));
for i =  1:length(lab)
    for j =  1:length(lab)
        clus.netScore(i,j) = length(find(sum(clus.edge(unique(clus.label(find(clus.network == j))),...
            unique(clus.label(find(clus.network == i))))))) ...
            / length(unique(clus.label(find(clus.network == i))));
    end
end

% get probability of network adjacency
count = 1;
for i = 1:length(lab)-1
    for j = i+1:length(lab)
        score = ((clus.netScore(lab(i),lab(j)) + clus.netScore(lab(j),lab(i))) / 2);
        if score ~= 0
            clus.score(count,:) = [lab(i) lab(j) score];
            count = count + 1;
        end
    end
end

% Write out csv file for reading into gephi
filename = ['./results/results.' num2str(length(lab)) '.csv'];
csvwrite(filename,[clus.score [1:length(clus.score)]'])

% find peak nodes for each cluster
for i = 1:max(clus.label)   
   c = find(clus.label == i); 
   while ~isempty(c)
       d = unique(edg(find(ismember(edg(:,1),c) + ismember(edg(:,2),c) == 1),:));
       c = setdiff(c,d);
   end
   clus.peaks(i) = d(1);
end

clus.regions = label;
