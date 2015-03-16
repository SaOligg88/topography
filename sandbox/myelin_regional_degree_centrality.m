function[s] = myelin_regional_degree_centrality(mask, dt)
% Set mask as nodes to run analysis on...
% dt = 10; % distance threshold for 'regional' aspect in mm.
% outputs 's' 
% example:
% s = myelin_regional_degree_centrality(10115 10116 10117], 10);

% load distance matrix	
dist = load('/scr/litauen1/dist....');
dist = dist.dist;

% set mask of area to extract value from.
%mask = []; % insert values here : [10015 10204 etc...]

% threshold dist
dist_thresh = {};
for i = 1:length(mask)
	dist_thresh{mask(i}) = find(dist(mask(i),:) > dt);
end

% grab subject list
List = dir('/a/documents/connectome/_all/');
for i = 1:length(List)-2
   subList(i,:) = List(i+2).name;
end

% grab resting state data
dim = 32492;
dirName = []'/a/documents/connectome/_all/'];
for i = 1:length(subList)
	count = 1;
	for sess = {'1', '2'};
		for pe = ['LR' 'RL'];
			filename = [dirName num2str(subList(i)) '/rfMRI_REST_left_corr_avg.gii.data']; % Change file name to rs data...
			data = cifti_read(filename);
			data = data.....; % find extension
			
			for m = 1:length(mask) 
				r = corr([data(mask(m),:) data(dist_thresh{mask(m)},:)]; ...
				% r to z transform TODO:
				z(count, m, :) = r(1,:);  
				count = count + 1;
			end
		end
	end 
	% set DIM
	s(i, :) = mean(z,DIM);
end

boxplot(s);

			
	
	

