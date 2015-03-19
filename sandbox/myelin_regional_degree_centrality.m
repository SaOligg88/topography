function[S] = myelin_regional_degree_centrality(mask, dt)
% Set mask as nodes to run analysis on...
% dt = 10; % distance threshold for 'regional' aspect in mm.
% outputs 'S' 
% example:
% S = myelin_regional_degree_centrality(10115 10116 10117], 10);

% load distance matrix	
dist = load('/scr/litauen1/dist.hcp.lh.mat');
dist = dist.data;

% set mask of area to extract value from.
%mask = []; % insert values here : mask = [10015 10208];

% threshold dist
dist_thresh = {};
for i = 1:length(mask)
	dist_thresh{mask(i)} = find(dist(mask(i),:) < dt & dist(mask(i),:) ~= 0);
end

% grab subject list
List = dir('/a/documents/connectome/_all/');
for i = 1:length(List)-2
   subList(i,:) = List(i+2).name;
end

% grab resting state data
dim = 32492;
sess = [1 2];
pe = ['LR'; 'RL'];
dir1 = ['/a/documents/connectome/_all/'];
clear S;
for i = 1:5%length(subList)
	count = 1;
	for s = 1:length(sess);
		for p = 1:length(pe);
            % read in data:
			filename = [dir1 num2str(subList(i,:)) ...
                '/MNINonLinear/Results/rfMRI_REST' ...
                num2str(sess(s)) '_' pe(p,:) '/rfMRI_REST' ...
                num2str(sess(s)) '_' pe(p,:) '_Atlas_hp2000_clean.dtseries.nii'];
			disp(filename);
            data = ft_read_cifti(filename);
			data = data.dtseries(1:32492,:); % for Left Hemi
			
			for m = 1:length(mask) 
                input = [data(mask(m),:); data(dist_thresh{mask(m)},:)]';
				r = corr(input);
                % take only r-values from surrounds nodes:
				R=r(1,find(r(1,:) ~= 1));
                % r to z tranform:
                z=.5.*log((1+R)./(1-R));
                % mean across local nodes:
                Z(count, m) = mean(z);
            end
            count = count + 1;
		end
	end 
	% mean across four runs within individual:
	S(i, :) = mean(Z,1);
end

h = figure;
boxplot(S);

			
	
	

