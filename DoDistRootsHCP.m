function[geoDistL, geoDistPrimaryL] = DoDistRootsHCP(sub)

addpath('/scr/litauen1/Dropbox/misc/surfstat/');
addpath('/scr/litauen1/Dropbox/misc/topology');

% import surface:
filename = ['/scr/dattel2/' num2str(sub) '/T1w/Native/' num2str(sub) ...
    '.L.midthickness.native.surf.gii'];
surf_gii = gifti(filename);
surf.coord = surf_gii.vertices';
surf.tri = surf_gii.faces;

filename = ['/scr/dattel2/' num2str(sub) '/T1w/Native/' num2str(sub) ...
    '.L.inflated.native.surf.gii'];
surf_gii = gifti(filename); surfi.coord = surf_gii.vertices'; surfi.tri = surf_gii.faces;
filename = ['/scr/dattel2/' num2str(sub) '/T1w/Native/' num2str(sub) ...
    '.L.white.native.surf.gii'];
surf_gii = gifti(filename); surfw.coord = surf_gii.vertices'; surfw.tri = surf_gii.faces;
clear surf_gii;

surfm.tri = surf.tri;
surfm.coord = (surfi.coord .* 0.7 + surfw.coord .* 0.3);

aparc = gifti(['/scr/dattel2/' num2str(sub) '/MNINonLinear/Native/' num2str(sub) ...
    '.L.aparc.a2009s.native.label.gii']);
aparc = aparc.cdata;

% Make noncortex mask
noncortex = zeros([1 length(surf.coord)]);
noncortex(find(aparc == 0)) = 1;

% % Make mask with anterior insula label
% labels = [18];
% mask = zeros([1 length(surf.coord)]);
% for i = 1:length(labels)
%     mask(find(aparc == labels(i))) = 1;
% end
% % add noncortex border
% edg = SurfStatEdg(surf);
% noncortEdg = noncortex(edg);
% border = unique([edg(find(eq(noncortEdg(:,1), noncortEdg(:,2)) == 0), 1); ...
%     edg(find(eq(noncortEdg(:,1), noncortEdg(:,2)) == 0), 2)]);
% mask(border) = 1;
% 
% % run distance
% geoDist = surfGeoDist(surf, mask, noncortex);
% geoDist(find(noncortex)) = 0;

%% Run for line:
% Make mask with anterior insula label
labels = [18];
maskIns = zeros([1 length(surf.coord)]);
for i = 1:length(labels)
    maskIns(find(aparc == labels(i))) = 1;
end
% add noncortex border
maskCing = zeros([1 length(surf.coord)]);
edg = SurfStatEdg(surf);
noncortEdg = noncortex(edg);
border = unique([edg(find(eq(noncortEdg(:,1), noncortEdg(:,2)) == 0), 1); ...
    edg(find(eq(noncortEdg(:,1), noncortEdg(:,2)) == 0), 2)]);
maskCing(border) = 1;
cortex = zeros([1 length(surf.coord)]);
cortex(find(noncortex == 0)) = 1;
maskCing = maskCing .* cortex;

% run distance
[geoDistIns surf] = surfGeoDist(surf, maskIns, noncortex);
[geoDistCing surf] = surfGeoDist(surf, maskCing, noncortex);
geoDist = ([geoDistIns; geoDistCing]);
[s ind] = sort(geoDist, 1);
lab = ind(1,:)';
%lab(find(noncortex)) = 0;

edg = SurfStatEdg(surf);
labE = lab(edg);
border = unique([edg(find(eq(labE(:,1), labE(:,2)) == 0), 1); edg(find(eq(labE(:,1), labE(:,2)) == 0), 2)]);
mask = zeros([length(surf.coord) 1]);
mask(border) = 1;

geoDist = s(1,:)';
geoDistL = s(1,:)';
%geoDistL(find(noncortex)) = 0;
geoDistL(find(mask)) = 0;

% visulize
% figure; SurfStatView(geoDist, surf);
% 
% %% Run for primary:
% Make mask with anterior insula label
% aparc(find(aparc == 29)) = 46; % include precentral gyrus in motor mask
labels = [46 45 75];
maskPrimary = zeros([1 length(surf.coord)]);
j = 1;
for i = 1:length(labels)
    maskPrimary(find(aparc == labels(i))) = 1;
    if sum(maskPrimary.*mask') > 0; % Removing from roots border
        maskPrimary1 = maskPrimary;
        maskPrimary1(find(lab == 2)) = 0;
        [geoDistP(j,:) surf] = surfGeoDist(surf, maskPrimary1, noncortex);
        maskPrimary2 = maskPrimary;
        maskPrimary2(find(lab == 1)) = 0;
        [geoDistP(j+1,:) surf] = surfGeoDist(surf, maskPrimary2, noncortex);
        j = j + 2;
    else
        [geoDistP(j,:) surf] = surfGeoDist(surf, maskPrimary, noncortex);
        j = j + 1;
    end
end

[s ind] = sort(geoDistP, 1);
labP = ind(1,:)';
edg = SurfStatEdg(surf);
labE = labP(edg);
border = unique([edg(find(eq(labE(:,1), labE(:,2)) == 0), 1); edg(find(eq(labE(:,1), labE(:,2)) == 0), 2)]);
mask = zeros([length(surf.coord) 1]);
mask(border) = 1;

geoDistPrimary = s(1,:)';
geoDistPrimaryL = s(1,:)';
%geoDistPrimaryL(find(noncortex)) = 0;
geoDistPrimaryL(find(mask)) = 0;
% % run distance
% %geoDistPrimary = surfGeoDist(surf, maskPrimary, noncortex)';
% 
% geoDistAll = [geoDistP; geoDistIns; geoDistCing];
% [s ind] = sort(geoDistAll, 1);
% lab = ind(1,:)';
% edg = SurfStatEdg(surf);
% labE = lab(edg);
% border = unique([edg(find(eq(labE(:,1), labE(:,2)) == 0), 1); edg(find(eq(labE(:,1), labE(:,2)) == 0), 2)]);
% mask = zeros([length(surf.coord) 1]);
% mask(border) = 1;
% geoDistA = s(1,:)';
% geoDistAL = s(1,:)';
% geoDistAL(find(mask)) = 0;
% 
% % take pic
% h = figure('visible','off'); SurfStatView(geoDistL, surfm); saveas(h, [num2str(sub) '.geoDist.png']); close all
% h = figure('visible','off'); SurfStatView(geoDistPrimaryL, surfm); saveas(h, [num2str(sub) '.geoDistPrimary.png']); close all
% h = figure('visible','off'); SurfStatView(labP, surfm); saveas(h, [num2str(sub) '.geoDistPrimaryClus.png']); close all
% h = figure('visible','off'); SurfStatView(lab, surfm); saveas(h, [num2str(sub) '.geoDistAll.png']); close all
% 
% % run patches analysis: 
% geoDistRoot = [geoDistIns; geoDistCing];
% [sR indR] = sort(geoDistRoot, 1);
% labR = indR(1,:)';
% 
% labels = [46 45 75];
% maskPrimary = zeros([1 length(surf.coord)]);
% for i = 1:length(labels)
%     maskPrimary(find(aparc == labels(i))) = 1;
%     [geoDistP(:,i) surf] = surfGeoDist(surf, maskPrimary, noncortex);
% end
% [sP indP] = sort(geoDistP', 1);
% labC = zeros([1 length(surf.coord)]);
% k = 1;
% for i = 1:length(labels)
%     for j = 1:length(labels)
%         labC(intersect(find(indP(1,:) == i),find(indP(2,:) == j))) = k;
%         k = k + 1;
%     end
% end
% labCR = labC' + ((labR-1)*k);
% 
% h = figure('visible','off'); SurfStatView(labP, surfm); saveas(h, [num2str(sub) '.geoDistPatches.png']); close all
% 
% 
% [surf max_degree] = surfGetNeighbors(surf);
% noncortexIDs = find(noncortex);
% surf.nbr(:,noncortexIDs) = 0;
% for i = 1:length(noncortexIDs)
%     surf.nbr(find(surf.nbr == noncortexIDs(i))) = 0;
% end
% % geoDistIns - maskCing
% 
% labels = [46 45 75];
% maskPrimary = zeros([length(labels) length(surf.coord)]);
% for i = 1:length(labels)
%     maskPrimary(find(aparc == labels(i)),i) = 1;
% end
% 
% a = maskCing .* geoDistIns;
% a = maskPrimary(:,1) 
% b = find(a == min(nonzeros(a)));
% edgeIC = b;
% while a(b) ~= 0 
%     c = find(surf.nbr(:,b) == min(nonzeros(surf.nbr(:,b))));
%     %b = surf.nbr(min(surf.nbr(b)),b);
%     if intersect(surf.nbr(c,b), edgeIC)
%         a(b) = 0;
%     else
%         edgeIC = [edgeIC; surf.nbr(c,b)];
%         disp(a(b));
%         b = surf.nbr(c,b); 
%     end
% end

