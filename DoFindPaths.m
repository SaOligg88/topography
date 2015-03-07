function[clus, countC, randpath] = DoPaths(n)
% n = number of permutations
% set to 0 to not run any purmutations
%

%% Load yeo preconfigured:
thresh=2;
[surf, surfi, surfm, label] = DoHCPsurf_group('L', 2);

edg = SurfStatEdg(surf);
lab = [1 2 3 4 6 7 8 9 10 11 13 14 15 16 17]; % skipping 5 and 12

%%%%%%%%% Begin %%%%%%%%%
addpath('./utils/surfstat/*');
slm = struct();
slm.tri = surf.tri';

clus.label = zeros([1 length(surf.coord)]);
countC = 0;
clus.label = zeros([1 32492]);
clus.network = zeros([1 32492]);
for i = 1:length(lab)    
    a = zeros([length(surf.coord) 1]);
    a(find(label == lab(i))) = 1; 
    slm.t = a';
    [cluster,clusid] = SurfStatPeakClus(slm,ones([length(surf.coord) 1]),0.5, ones(1,length(surf.coord)), edg);  
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
end

clus.edge = zeros(countC);
for i = 1:countC
    a = [edg(find(ismember(edg(:,1),find(clus.label == i))),2); ...
        edg(find(ismember(edg(:,2),find(clus.label == i))),1)];
    clus.edge(i,...
        nonzeros(unique(clus.label(a))))...
        = 1;
end

clus.netScore = zeros(17);
for i = lab
    for j = lab
        clus.netScore(i,j) = length(find(sum(clus.edge(unique(clus.label(find(clus.network == j))),...
            unique(clus.label(find(clus.network == i))))))) ...
            / length(unique(clus.label(find(clus.network == i))));
    end
end

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

if n > 0
	addpath('./utils');
	randPath = pathsPermute(clus, n);
else
	randPath = [];
end


        
