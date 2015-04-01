function[vtx,fac] = WriteSurfMap(data, filename, rgb, surf)

vtx = surf.coord';
fac = surf.tri;

% Prepare vertices, colours and faces
rgb = double(uint8(255*rgb));
vtx2ply = [vtx rgb(data+1,1) rgb(data+1,2) rgb(data+1,3)];
fac2ply = [ones(size(fac,1),1)*3 fac-1];

% Save surface (PLY)
plyfile = sprintf('%s.ply',filename);
fid = fopen(plyfile,'w');
fprintf(fid,'ply\n');
fprintf(fid,'format ascii 1.0\n');
fprintf(fid,'element vertex %s\n',num2str(size(vtx,1)));
fprintf(fid,'property float x\n');
fprintf(fid,'property float y\n');
fprintf(fid,'property float z\n');
fprintf(fid,'property uchar red\n');
fprintf(fid,'property uchar green\n');
fprintf(fid,'property uchar blue\n');
fprintf(fid,'element face %s\n',num2str(size(fac,1)));
fprintf(fid,'property list uchar int vertex_index\n');
fprintf(fid,'end_header\n');
fprintf(fid,'%g %g %g %d %d %d\n',vtx2ply');
fprintf(fid,'%d %d %d %d\n',fac2ply');
fclose(fid);

