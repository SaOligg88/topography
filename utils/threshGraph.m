function[matrix] = threshGraph(graph)

addpath(genpath('./utils'));

thresh = [0:0.05:1];
for i = 1:length(thresh)
	graphThreshed = (graph > thresh(i)) .* graph;
	S = isconnected(graphThreshed);
	if S == 0
		matrix = (graph > thresh(i - 1)) .* graph;
		break; 
	end
	disp(thresh(i));
end
