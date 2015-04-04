function [distDMN, zonesDMN] = makeDMNdist_yeo()

addpath(genpath('./utils'));
rsnL = gifti('/scr/litauen1/Dropbox/misc/yeoTopo/RSN_L.gii');
rsnR = gifti('/scr/litauen1/Dropbox/misc/yeoTopo/RSN_R.gii');
rsn = [rsnL.cdata(:,1:2); rsnR.cdata(:,1:2)];
rsn7 = rsn(:,1);
labs = unique(rsn7);
for i = 1:length(labs)
    rsn7(find(rsn7 == labs(i))) = i-1;
end

[surfL, surfiL, surfmL] = loadHCPsurf_group('L');
[surfR, surfiR, surfmR] = loadHCPsurf_group('R');

[clusL7] = pathsFindHCP(rsn7(1:32492),3,'L', surfL);
[clusR7] = pathsFindHCP(rsn7(32493:32492*2),3,'R', surfR);

rsn17 = rsn(:,2);
labs = unique(rsn17);
for i = 1:length(labs)
    rsn17(find(rsn17 == labs(i))) = i-1;
end

[clusL17] = pathsFindHCP(rsn17(1:32492),3,'L', surfmL);
[clusR17] = pathsFindHCP(rsn17(32493:32492*2),3,'R', surfmR);

% hippocampal DMN node:
DMNpeaksL = [clusL17.peaks(find(clusL17.edgeNet == 2)) clusL7.peaks(28)];
DMNpeaksR = [clusR17.peaks(find(clusR17.edgeNet == 2)) clusR7.peaks(19)];

[distDMNL, zonesDMNL] = distExactGeodesic(DMNpeaksL, '32', 'L', 'all', '100307');
[distDMNR, zonesDMNR] = distExactGeodesic(DMNpeaksR, '32', 'R', 'all', '100307');

% make maps
DMNclusL  = zeros(32492,1);
toFind = find(clusL17.edgeNet == 2);
for  i = 1:length(toFind)
    DMNclusL(find(clusL17.label == toFind(i))) = i;
end
DMNclusL(find(clusL7.label == 28)) = 7;
DMNclusL([clusL17.peaks(toFind) clusL7.peaks(28)]) = 8;
change = [1 3 2 6 7 5 4 ];
DMNclusL = DMNclusL .* 10;
for i = 1:7
    DMNclusL(find(DMNclusL == i * 10)) = change(i);
end
DMNclusL(find(DMNclusL == 80)) = 8;
DoFlatten(DMNclusL', 'hcp.zones.clus.DMN.L', 'label_col_left', '32');

DMNclusR = zeros(32492,1);
toFind = find(clusR17.edgeNet == 2);
for  i = 1:length(toFind)
    DMNclusR(find(clusR17.label == toFind(i))) = i;
end
DMNclusR(find(clusR7.label == 19)) = 7;
DMNclusR([clusR17.peaks(toFind) clusR7.peaks(19)]) = 8;
change = [3 1 2 6 7 5 4 ];
DMNclusR = DMNclusR .* 10;
for i = 1:7
    DMNclusR(find(DMNclusR == i * 10)) = change(i);
end
DMNclusR(find(DMNclusR == 80)) = 8;
DoFlatten(DMNclusR', 'hcp.zones.clus.DMN.R', 'label_col_right', '32');

distDMN = min([distDMNL; distDMNR], [], 2);
zonesDMN = min([zonesDMNL; zonesDMNR], [], 2);

