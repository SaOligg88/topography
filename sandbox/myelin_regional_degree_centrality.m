function[s] = myelin_regional_degree_centrality(mask, dt)
% Set mask as nodes to run analysis on...
% dt = 10; % distance threshold for 'regional' aspect in mm.
% outputs 's' 
% example:
% s = myelin_regional_degree_centrality(10115 10116 10117], 10);

% load distance matrix	
dist = load('/scr/litauen1/dist.hcp.lh.mat');
dist = dist.data;

% set mask of area to extract value from.
%mask = []; % insert values here : mask = [10015 10204];

% threshold dist
dist_thresh = {};
for i = 1:length(mask)
	dist_thresh{mask(i)} = find(dist(mask(i),:) < dt);
end

% grab subject list
List = dir('/a/documents/connectome/_all/');
for i = 1:length(List)-2
   subList(i,:) = List(i+2).name;
end

% grab resting state data
dim = 32492;
dir1 = ['/a/documents/connectome/_all/'];
for i = 1:length(subList)
	count = 1;
	for sess = ['1', '2'];
		for pe = ['LR' 'RL'];
			filename = [dir1 num2str(subList(i,:)) ...
                '/MNINonLinear/Results/rfMRI_REST' sess '_' pe '/rfMRI_REST' sess '_' pe '_Atlas_hp2000_clean.dtseries.nii']; % Change file name to rs data...
			data = ft_read_cifti(filename);
			data = data.dtseries(1:32492,:); % find extension
			
			for m = 1:length(mask) 
                input = [data(mask(m),:); data(dist_thresh{mask(m)},:)]';
				r = corr(input);
				% r to z transform TODO:
				% z(count, m, :) = r(1,:); 
                z(count, m, :) = mean(r(1,isfinite(r(1,:)),:))
				count = count + 1;
			end
		end
	end 
	% set DIM
	s(i, :) = mean(z,DIM);
end

boxplot(s);

			
	
	

