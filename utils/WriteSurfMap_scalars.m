function[vtx,fac] = WriteSurfMap(data, filename, surf, range)

vtx = surf.coord';
fac = surf.tri;
srfwrite(vtx,fac,'tmp.srf');
dpxwrite('tmp.dpv', data, vtx, 1:length(vtx));
dpx2map('tmp.dpv', 'tmp.srf',   filename, 'jet', range, ['[0.0001 ' num2str(max(data)) ']']);
!rm tmp.dpv tmp.srf

