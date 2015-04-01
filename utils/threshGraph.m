function[matrix, s, ind, x, y] = threshGraph(clus, surf, filename)

addpath(genpath('./utils'));

% graph = clus.netScore;
graph = edgeL2adj(clus.score) + edgeL2adj(clus.score)';

thresh = [0:0.001:1];
for i = 1:length(thresh)
	graphThreshed = (graph > thresh(i)) .* graph;
	S = isconnected(graphThreshed);
	if S == 0
		matrix = (graph > thresh(i - 1)) .* graph;
		break; 
	end
	% disp(thresh(i));
end

% graphViz4Matlab('-adjMat', tril(matrix), '-undirected', '[true]');
% graphViz4Matlab('-adjMat', matrix);
% drawNetwork(matrix);
 
% The vector corresponding to the second smallest eigenvalue of the Laplacian matrix:
[V,D]=eig(laplacian_matrix(graph));
[ds,Y]=sort(diag(D));
fv=V(:,Y(2));
[s,ind] = sort(fv);
% [s,ind] = sort(fiedler_vector(edgeL2adj(clus.score) + edgeL2adj(clus.score)'));
figure; plot(s, '+');
set(gca,'XTick',[1:length(s)]);
set(gca,'XTickLabel',num2str(ind));
xlim([0 length(s)+1]);

%% plot smallest
label = clus.regions;
for j = 2;
    fv=V(:,Y(j));
    % boxcox transform
    % fv = boxcox(fv + 1);
    [s,ind] = sort(fv);

    a = zeros(length(surf.coord),1);
    for i = 1:length(unique(nonzeros(label)))
        a(find(label == i)) = (find(ind == i));
    end
    figure; SurfStatView(a, surf, [num2str(j)]);
end

%% Plot two smallest eigenvectors in scatterplot
h = figure;
x = boxcox(V(:,Y(2))+1); 
y = abs(boxcox(V(:,Y(3))+1));
%x = abs(boxcox(V(:,Y(2))+1)); y = abs(boxcox(V(:,Y(3))+1));
scatter(x,y,'.')
a = [1:length(unique(nonzeros(label)))]'; b = num2str(a); c = cellstr(b);
dx = 0.005; dy = 0.00; % displacement so the text does not overlay the data points
text(x+dx, y+dy, c);
xlim([min(x)-0.05 max(x)+0.05]);
ylim([min(y)-0.05 max(y)+0.05]);

%% write out order
[pcccoeff, pcvec] = pca([x y]');
[~,ind] = sort(pcvec(:,1));
csvwrite([filename '.labelorder'], ind);
csvwrite([filename '.labels'], label);


