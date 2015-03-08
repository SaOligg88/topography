function[matrix, s, ind] = threshGraph(clus)

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

