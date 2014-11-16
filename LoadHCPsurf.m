function[surf] = LoadHCPsurf(num)

List=dir('/a/documents/connectome/_all/');
for i = 1:length(List)-2
    ListSub(i,:) = List(i+2).name;
end

sub = ListSub(num);

% import surface:
filename = ['/scr/dattel2/' num2str(sub) '/T1w/Native/' num2str(sub) ...
    '.L.midthickness.native.surf.gii'];
surf_gii = gifti(filename);
surf.coord = surf_gii.vertices';
surf.tri = surf_gii.faces;
clear surf_gii;