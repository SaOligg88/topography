function[surf, surfi, surfvi] = LoadHCPsurf_individual(sub, hemi, res)
% res = '32' or '164'

% List=dir('/a/documents/connectome/_all/');

% for i = 1:length(List)-2
%    ListSub(i,:) = List(i+2).name;
% end
% sub = ListSub(num);

dir = ['/a/documents/connectome/_all/'];
% import surface:
switch res
	case '32'

		filename = [dir num2str(sub) '/MNINonLinear/fsaverage_LR32k/' num2str(sub) '.' hemi '.midthickness.32k_fs_LR.surf.gii']);
		surf_gii = gifti(filename); 
		surf.coord = surf_gii.vertices';
		surf.tri = surf_gii.faces;

		filename = [dir num2str(sub) '/MNINonLinear/fsaverage_LR32k/' num2str(sub) '.' hemi '.inflated.32k_fs_LR.surf.gii']);
		surf_gii = gifti(filename); 
		surfi.coord = surf_gii.vertices';
		surfi.tri = surf_gii.faces;

		filename = [dir num2str(sub) '/MNINonLinear/fsaverage_LR32k/' num2str(sub) '.' hemi '.very_inflated.32k_fs_LR.surf.gii']);
		surf_gii = gifti(filename); 
		surfvi.coord = surf_gii.vertices';
		surfvi.tri = surf_gii.faces;

	case '164'

		filename = [dir num2str(sub) '/MNINonLinear/' sub '.' hemi '.midthickness.164k_fs_LR.surf.gii']);
		surf_gii = gifti(filename); 
		surf.coord = surf_gii.vertices';
		surf.tri = surf_gii.faces;

		filename = [dir num2str(sub) '/MNINonLinear/' sub '.' hemi '.inflated.164k_fs_LR.surf.gii']);
		surf_gii = gifti(filename); 
		surfi.coord = surf_gii.vertices';
		surfi.tri = surf_gii.faces;

		filename = [dir num2str(sub) '/MNINonLinear/' sub '.' hemi '.very_inflated.164k_fs_LR.surf.gii']);
		surf_gii = gifti(filename); 
		surfvi.coord = surf_gii.vertices';
		surfvi.tri = surf_gii.faces;
end