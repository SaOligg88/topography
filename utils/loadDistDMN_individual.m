function[distances, zones, surfi_164] = loadDistDMN_individual(sub, hemi)

[surf, surfi, surfvi] = loadHCPsurf_individual(sub, hemi, '32')
[surf_164, surfi_164, surfvi_164] = loadHCPsurf_individual(sub, hemi, '164');

dmn = eigenvector centrality map...;

%% threshold:
edg = SurfStatEdg(surf);
slm = struct();
slm.tri = surf.tri';
slm.t = dmn;
[peak,clus,clusid] = SurfStatPeakClus(slm,(vert(:,1) ==3),0.001, ones(1,length(surf.coord)), edg);
% former mask: ones([length(surf.coord) 1])

peaks = [];
for i = 1:length(clus.clusid)
    ind = find(peak.clusid == i);
    [m im] = max(peak.t(ind));
    peaks(i) = peak.vertid(ind(im));
end
% a = zeros(32492,1); a(peaks) = 1; figure; SurfStatView(a, surf);

c = SurfStatInd2Coord(peaks, surfvi);
inds = SurfStatCoord2Ind(c', surfvi_164);
    
[distances,zones] = distExactGeodesic(inds, '164', hemi, 'zones', sub);

figure; SurfStatView(distances, surfvi_164);
figure; SurfStatView(zones, surfvi_164);
