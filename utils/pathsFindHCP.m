function[clus] = pathsFindHCP(label, thresh, hemi, surf)
% 'label' is a vector of the length of the surface with integers denoting
% groups. 
% 'thresh' is equal to the minimum number of adjacent nodes that can
% consistiute a cluster. This was put in to remove single-node clusters. 
% I generally set it somewhere between 0 and 2.

addpath(genpath('./utils'));

lab = nonzeros(unique(label));
% [surf, surfi, surfm] = loadHCPsurf_group(hemi);

edg = SurfStatEdg(surf);
clus.label = zeros(1,length(surf.coord));
countC = 0;
clus.label = zeros(1,length(surf.coord));
clus.network = zeros(1,length(surf.coord));
for i = 1:length(lab)    
    a = zeros(length(surf.coord),1);
    a(find(label == lab(i))) = 1; 
    slm = struct();
    slm.tri = surf.tri';
    slm.t = a';
    [cluster,clusid] = SurfStatPeakClus(slm,ones(length(surf.coord),1),0.05, ones(1,length(surf.coord)), edg);  
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
	mem1 = find(ismember(edg(:,1),find(clus.label == i)));
	mem2 = find(ismember(edg(:,2),find(clus.label == i)));
    a = [edg(mem1,2); edg(mem2,1)];
    clus.edge(i,nonzeros(unique(clus.label(a)))) = 1;
end
clus.edge(find(eye(size(clus.edge)))) = 0;

for i = 1:length(clus.edge)
    clus.edgeNet(i) = unique(clus.network(clus.label == i));
end

clus.netScore = zeros(length(lab));
for i =  1:length(lab)
    for j =  1:length(lab)
		netlabs1 = unique(clus.label(find(clus.network == i)));
		netlabs2 = unique(clus.label(find(clus.network == j)));
		% should this be sum or prod?:
        clus.netScore(i,j) = length(find(sum(clus.edge(netlabs2,netlabs1)))) / length(netlabs1);
    end
end

% get probability of network adjacency
count = 1;
for i = 1:length(lab)-1
    for j = i+1:length(lab)
        score = ((clus.netScore(i,j) + clus.netScore(j,i)) / 2);
        if score ~= 0
            clus.score(count,:) = [lab(i) lab(j) score];
            count = count + 1;
        end
    end
end

% Write out csv file for reading into gephi
%filename = ['./results/results.' num2str(length(lab)) '.csv'];
%csvwrite(filename,[clus.score [1:length(clus.score)]'])

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
