function[thresh, surf, surfi, surfm, label, edg, lab] = config_yeo_L_17rsns();

thresh=2;
[surf, surfi, surfm] = loadHCPsurf_group('L');
label = [];
edg = SurfStatEdg(surf);
lab = [1 2 3 4 6 7 8 9 10 11 13 14 15 16 17]; % skipping 5 and 12
