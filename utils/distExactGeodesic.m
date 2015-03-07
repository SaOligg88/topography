function [dist, zone] = exactGeodesic(source, surfType, hemi, analysisType)
% source = source node ids
% surfType = '32' or '164'
% hemi = 'L' or 'R'
% analysisType = 'distance' 'zones' 'all'

addpath('/scr/litauen1/Dropbox/misc/yeoTopo/lme/geodesic/')

if surfType == '32'
	surf = gifti(['/scr/murg2/HCP_new/HCP_Q1-Q6_GroupAvg_Related440_Unrelated100_v1/Q1-Q6_R440.' ...
		hemi '.midthickness.32k_fs_LR.surf.gii']);
	sub = num2str(100307);
	aparc = gifti(['/a/documents/connectome/_all/' sub '/MNINonLinear/fsaverage_LR32k/' sub ...
		'.' hemi '.aparc.a2009s.32k_fs_LR.label.gii']);
	aparc = aparc.cdata;
	noncortex = find(aparc == 0);

	% To remove mesh problem
	if hemi == 'R'
    	noncortex = sort([noncortex; 27806; 27810; 128; 129; 27779; 27792; 27791; 27798; ...
    		27785; 27786; 27793; 27799; 27804; 27809; 27805; 27799; 27800; 27806; 27810; 27810 ]);
	end
	
elif surfType == '164'
	filename = ['/a/documents/connectome/_all/' sub '/MNINonLinear/' sub '.L.midthickness.164k_fs_LR.surf.gii'];
	surf = gifti(filename);
	filename = ['/a/documents/connectome/_all/' sub '/MNINonLinear/' sub '.L.very_inflated.164k_fs_LR.surf.gii'];
	surf_gii = gifti(filename);
	surfvi.coord = surf_gii.vertices'; surfvi.tri = surf_gii.faces;

	aparc = gifti(['/a/documents/connectome/_all/' sub '/MNINonLinear/' sub ...
		'.L.aparc.a2009s.164k_fs_LR.label.gii']);
	aparc = aparc.cdata;
	noncortex = find(aparc == 0);
end

index1 = ismember(surf.faces(:,1),noncortex);
index2 = ismember(surf.faces(:,2),noncortex);
index3 = ismember(surf.faces(:,3),noncortex);
index = index1 + index2 + index3;
surfN.faces = surf.faces(index == 0,:);
[incld r1 r2] = unique(surfN.faces);
surfN.vertices = surf.vertices(incld,:);
new = 1:length(incld);
surfN.faces = reshape(new(r2), size(surfN.faces));
vertices=surfN.vertices;
faces=surfN.faces;
N = length(vertices);

disp('Ready to roll');
%% Run 
global geodesic_library;                
geodesic_library = 'libgeodesic';
mesh = geodesic_new_mesh(vertices,faces);
algorithm = geodesic_new_algorithm(mesh, 'exact');
source_points = {};
s = find(ismember(incld, source))';
for vertex_id = s;
    source_points{end+1} = geodesic_create_surface_point('vertex',vertex_id, vertices(vertex_id,:));
end;
geodesic_propagate(algorithm, source_points);

switch analysisType
	case 'distance'
		[~, distances] = geodesic_distance_and_source(algorithm);
		dist = zeros(length(surf.vertices),1);
		dist(incld) = distances;
		zone = [];
		
	case 'zones'
		distances = zeros(N,1); source_ids = zeros(N,1);
		for i=1:N   
		    q = geodesic_create_surface_point('vertex',i,vertices(i,:));
		    [source_ids(i), distances(i)] = geodesic_distance_and_source(algorithm, q);    
		end;
		dist = zeros(length(surf.vertices),1);
		dist(incld) = distances;
		zone = zeros(length(surf.vertices),1);
		zone(incld) = source_ids;
		
	case 'all'
		distances = zeros([N length(source)]);
		source_ids = zeros([N length(source)]);
		for i = 1:length(s)
		    q = {geodesic_create_surface_point('vertex',s(i), vertices(s(i),:))};
		    geodesic_propagate(algorithm, q);
		    [source_ids(:,i), distances(:,i)] = geodesic_distance_and_source(algorithm);
		end
		dist = zeros(length(surf.vertices),length(source));
		dist(incld,:) = distances;
		zone = zeros(length(surf.vertices),length(source));
		zone(incld,:) = source_ids;
end
geodesic_delete;
