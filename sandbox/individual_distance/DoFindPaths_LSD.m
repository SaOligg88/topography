function[dist2mask] = DoFindPaths_LSD(subjects)

% set variables:
addpath(genpath('../../utils'));
thresh = 2;		
hemi   = {'lh'}; 	% {'lh','rh'}
% locations of freesurfer directory containing subjects:
dir1   = '/scr/ilz2/LEMON_LSD/freesurfer/';
dist2mask = struct();

for s = 1:length(subjects)
    subject = num2str(subjects(s));
    for h = 1:length(hemi)
        % load standard fsaverage5
        surf = SurfStatReadSurf([dir1 'fsaverage5/surf/' hemi{h} '.inflated']); 
        % Load cluster results: (output from individual_distance_cluster.py)
        % Check!!!:
        results = load([dir1 'clusters_' subject '_' hemi{h} '_em25_2_20.mat']);
        [~, labelannot, colortable] = read_annotation([dir1 subject '/label/' hemi{h} '.aparc.a2009s.annot']);

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
            dmnMasks = [27 ... G_pariet_inf-Supramar
                        57 ... S_interm_prim-Jensen
                        75 ... S_temporal_sup
                        ];
            clusDMNpar = unique(nonzeros(ismember(labelannot, colortable.table([dmnMasks],end)) .* clus.label .* find(clus.network == clusDMN)));
            
            if ~isempty(clusDMNpar)
                clusFound = 1;
            elseif length(clusDMNpar) > 1
                % adjudicate between conflicting clusters by taking
                % furthest posterior
                for i = 1:length(clusDMNpar)
                    realDMN(i) = mean(surf.coord(2,clus.label == clusDMNpar));
                end               
                [~,indRealDMN] = min(realDMN);
                clusDMNpar = clusDMNpar(indRealDMN);
                c = c+1; 
            end
        end
        dist2mask.clusDMN(s,h) = clusDMN;
        dist2mask.clusDMNpar(s,h) = clusDMNpar;
        dist2mask.clusterNum(s,h,:) = clust(c);
        % Transform region from fsaverage space to individual space
            % Where freesurfer grabs the data, does it do so from reg space of individual?
        surf_sphere = SurfStatReadSurf([dir1 'fsaverage5/surf/' hemi{h} '.sphere']); 
        surf_sphere_ind = SurfStatReadSurf([dir subject '/surf/' hemi '.sphere.reg']);
        coords = SurfStatInd2Coord(find(clus.label == clusDMNpar), surf_sphere);
        source = unique(SurfStatCoord2Ind(coords, surf_sphere_ind));

        % Get distances from DMN:
        [distanceMap, ~] = distExactGeodesic(source, 'freesurfer', hemi, 'distance', [dir1 subject]);
        dist2mask.distanceMap(s,h,:) = distanceMap;
        % Get min dist value to masks from primary (using aparc)
        % OR from other end of clust path?
        mask = [45 ... % S_calcarine
                76 ... % S_temporal_transverse
                ];             
        for i = 1:length(mask)
            % min or mean or median?   
            dist2mask.dist(s,h,i) = min(distanceMap(ismember(labelannot, colortable.table(mask(i),end))));
        end
        dist2mask.hemi{h} = hemi{h};
    end
    dist2mask.subjects(s) = str2double(subject);
end


       