function[matrix] = threshGraph(graph)

thresh = [0:0.05:1];
matrix = zeros(size(graph));
for i = 1:length(thresh)
	graphThreshed = zeros(size(graph));
	graphThreshed = graph(find(graph > thresh(i)));
	S = isConnected(graphThreshed);
	if S == 0
		matrix = graph(find(graph > thresh(i-1)));
		break; 
	end
	disp(i);
end
