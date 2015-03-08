surf_gii = gifti('data/Q1-Q6_R440.L.midthickness.32k_fs_LR.surf.gii');
surf.coord = surf_gii.vertices'; surf.tri = surf_gii.faces;
[label] = loadRSNs('L', 2);
thresh = 2;
lab = [1 2 3 4 6 7 8 9 10 11 13 14 15 16 17]; % skipping 5 and 12

lab = nonzeros(unique(label));
lab = lab(2:length(lab));
lab = lab([1 2 3 4 6 7 8 9 10 11 13 14 15 16 17]);
label_new = zeros(length(surf.coord),1);
for i = 1:length(lab)
	label_new(find(label == lab(i))) = i; 
end
label = label_new;