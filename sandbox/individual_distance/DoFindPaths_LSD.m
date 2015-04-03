function[dist2mask] = DoFindPaths_LSD(subjects)

% set variables:
addpath(genpath('../../utils'));
thresh = 2;		
hemi   = ['lh']; 	% ['lh','rh']
% locations of freesurfer directory containing subjects:
dir1   = '/scr/ilz2/LEMON_LSD/freesurfer/';

for s = 1:length(subjects)
    subject = num2str(subjects(s));
    for h = 1:length(hemi)
        % load standard fsaverage5
        surf = SurfStatReadSurf([dir1 'fsaverage5/surf/' hemi(h) '.inflated']); 
        % Load cluster results: (output from individual_distance_cluster.py)
        results = load([dir1 'cluster_' subject '_' hemi(h) '_em25_2_20.mat']);
        [~, labelannot, colortable] = read_annotation([dir1 subject '/label/' hemi '.aparc.a2009s.annot']);

        %%%%%%%%% Begin %%%%%%%%%
        % Decide on which number cluster solution:
        clusFound = 0;
        % order to check clust number:
        clust = [15 16 17 18 19 14 13 12 11 10 9 8 7 6 5 4 3 2 1];
        c = 1;
        while clusFound == 0
            disp(['trying ' num2str(clust(c)+1) ' cluster solution');
            label = results.results(clust(c),:);
            clus = pathsFindHCP(label, thresh, hemi, surf);
            graph = edgeL2adj(clus.score) + edgeL2adj(clus.score)';
            thresh = [0:0.001:1];
            for i = 1:length(thresh)
                graphThreshed = (graph > thresh(i)) .* graph;
                S = isconnected(graphThreshed);
                if S == 0                
                    break; exact
                end           
            end
            % The vector corresponding to the second smallest eigenvalue of the Laplacian matrix:
            [V,D]=eig(laplacian_matrix(graph));
            [~,Y]=sort(diag(D));
            [~,ind] = sort(V(:,Y(2)));

            % Get DMN
            if sum(clus.edgeNet == ind(end)) > sum(clus.edgeNet == ind(1))
                clusDMN = ind(end);
            else
                clusDMN = ind(1);
            end

            % Get cluster of interest in parietal         
            % 26 and 27 denote parietal % CHECK!!!
            clusDMNpar = unique(nonzeros(ismember(labelannot, colortable.table([26 27],end)) .* clus.label .* find(clus.network == clusDMN)));

            if ~isempty(clusDMNpar)
                clusFound = 1;
            else
                c = c+1; 
            end
        end
        % Transform region from fsaverage space to individual space
            % Where freesurfer grabs the data, does it do so from reg space of individual?
        surf_sphere = SurfStatReadSurf([dir1 'fsaverage5/surf/' hemi(h) '.sphere']); 
        surf_sphere_ind = SurfStatReadSurf([dir subject '/surf/' hemi '.sphere.reg']);
        coords = SurfStatInd2Coord(find(clus.label == clusDMNpar), surf_sphere);
        source = unique(SurfStatCoord2Ind(coords, surf_sphere_ind));

        % Get distances from DMN:
        [distanceMap, ~] = distExactGeodesic(source, 'freesurfer', hemi, 'distance', [dir1 subject]);

        % Get min dist value to masks from primary (using aparc)
        mask = []; % label values from aparc, OR from other end of clust path?
        for i = 1:length(mask)
            % min or mean or median?   
            dist2mask(s,h,i) = min(distanceMap(ismember(labelannot, colortable.table(mask(i),end))));
        end
    end
end
