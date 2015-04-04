function [dist2mask] = DoFindPaths_LSD(subjects, varargin)

% The inputs are:
% subject:      can be single subject 123456, or [123456, 124325, 23423452]
% The output of dist2mask is  
% .clusDMN:     identifier of network identifier used
% .clusDMNpar:  identifier of DMN region used
% .clusterNum:  cluster number used 
% .distanceMap: full distance map
% .dist:        specific min distances to masked regions
% .distlabels:  names of masked regions
% .hemi:        corresponding hemisphere to dist
% .subjects:    subject id

% set variables:
addpath(genpath('../../utils'));

p = inputParser;
defaultThresh = 2; % remove clusters smallest than this number of nodes
defaultHemi   = {'lh','rh'}; 	% {'lh','rh'}
% locations of freesurfer directory containing subjects:
dir1   = '/scr/ilz2/LEMON_LSD/freesurfer/';
dir2   = '/scr/liberia1/';
defaultFigures = false; % can be changed on the commandline with ...,'figures',true)

addRequired(p,'subects',@isnumeric);
addParameter(p,'thresh',defaultThresh,@isnumeric);
addParameter(p,'figures', defaultFigures);
addParameter(p,'hemi', defaultHemi);

parse(p,subjects,varargin{:}); 
thresh = p.Results.thresh;
hemi = p.Results.hemi;

dist2mask = struct();
for s = 1:length(subjects)
    subject = num2str(subjects(s));
    for h = 1:length(hemi)
        % load standard fsaverage5
        surf = SurfStatReadSurf([dir1 'fsaverage5/surf/' hemi{h} '.inflated']); 
        % Load cluster results: (output from individual_distance_cluster.py)
        % Check!!!:
        results = load([dir2 'clusters_' subject '_' hemi{h} '_em25_2_20.mat']);
        [~, labelannot, colortable] = read_annotation([dir1 'fsaverage5/label/' hemi{h} '.aparc.a2009s.annot']);

        %%%%%%%%% Begin %%%%%%%%%
        % Decide on which number cluster solution:
        clusFound = 0;
        % order to check clust number:
        clust = [size(results.results,1):-1:1];%[10 11 12 13 14 15 16 17 18 19  9 8 7 6 5 4 3 2 1];      
        c = 1;
        while clusFound == 0
            disp(['trying ' num2str(clust(c)+1) ' cluster solution']);
            label = results.results(clust(c),:);
            clus = pathsFindHCP(label, thresh, hemi, surf);
            graph = edgeL2adj(clus.score) + edgeL2adj(clus.score)';
            thr = [0:0.001:1];
            for i = 1:length(thr)
                graphThreshed = (graph > thr(i)) .* graph;
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
            dmnMasks = [75 ... % S_temporal_sup 57 ... % S_interm_prim-Jensen                         
                        27 ... % G_pariet_inf-Supramar                       
                        ];
            clusDMNpar = unique(nonzeros(ismember(labelannot, colortable.table(dmnMasks,end))' ...
                .* (clus.label .* (clus.network == clusDMN))));
            
            if length(clusDMNpar) == 1
                clusFound = 1;
            elseif length(clusDMNpar) > 1
                % adjudicate between conflicting clusters by taking
                % furthest dorsal
                realDMN = [];
                for i = 1:length(clusDMNpar)
                    realDMN(i) = mean(surf.coord(3,clus.label == clusDMNpar(i)));
                end               
                [~,indRealDMN] = max(realDMN);
                clusDMNpar = clusDMNpar(indRealDMN);
                if length(clusDMNpar) == 1
                    clusFound = 1;         
                end
            end
            c = c+1; 
        end
        dist2mask.clusDMN(s,h) = clusDMN;
        dist2mask.clusDMNpar(s,h) = clusDMNpar;
        dist2mask.clusterNum(s,h,:) = clust(c);
        dist2mask.labelNetwork(s,h,:) = clus.network;
        dist2mask.labelClus(s,h,:) = clus.label;
        
        % Transform region from fsaverage space to individual space
            % Where freesurfer grabs the data, does it do so from reg space of individual?
        surf_sphere = SurfStatReadSurf([dir1 'fsaverage5/surf/' hemi{h} '.sphere']); 
        surf_sphere_ind = SurfStatReadSurf([dir1 subject '/surf/' hemi{h} '.sphere.reg']);
        coords = SurfStatInd2Coord(find(clus.label == clusDMNpar), surf_sphere);
        source = unique(SurfStatCoord2Ind(coords', surf_sphere_ind));

        % Get distances from DMN:
        [distanceMap, ~] = distExactGeodesic(source, 'freesurfer', hemi{h}, 'distance', [dir1 subject]);
        dist2mask.distanceMap{s,h,:} = distanceMap;
        % Get min dist value to masks from primary (using aparc)
        % OR from other end of clust path?
        [~, labelannot, colortable] = read_annotation([dir1 subject '/label/' hemi{h} '.aparc.a2009s.annot']);
        mask = [46 ... % S_calcarine
                76 ... % S_temporal_transverse
                ];             
        for i = 1:length(mask)
            % min or mean or median?   
            dist2mask.dist{s,h,i} = min(distanceMap(ismember(labelannot, colortable.table(mask(i),end))));
            dist2mask.distlabels{s,h,i} = colortable.struct_names{mask(i)};
        end       
        dist2mask.hemi{h} = hemi{h};        
    
    end    
    dist2mask.subjects(s) = str2double(subject);
    
    if true(p.Results.figures)
        surf_ind = SurfStatReadSurf({[dir1 subject '/surf/lh.inflated'], [dir1 subject '/surf/rh.inflated']});
        % look at results: 
        
        h = figure; SurfStatView([dist2mask.distanceMap{s,1,:}; dist2mask.distanceMap{s,2,:}], surf_ind);
        saveas(h, [subject '.distanceMap.DMN.parietal.' clust(c-1) '.png']);
        clf
        surf = SurfStatReadSurf({[dir1 'fsaverage5/surf/lh.inflated'], [dir1 'fsaverage5/surf/rh.inflated']}); 
        title = ['left: ' num2str(dist2mask.clusDMN(s,1)) ' right: ' num2str(dist2mask.clusDMN(s,2))];
        SurfStatView([squeeze(dist2mask.labelNetwork(s,1,:)); squeeze(dist2mask.labelNetwork(s,2,:))], surf, title);
        saveas(h, [subject '.clusters.' clust(c-1) '.png']); 
        close all
    end
end
save('dist2mask.mat','-v7.3','dist2mask');
