function[thresh, surf, surfi, surfm, label, edg, lab] = config();

thresh=2;
[surf, surfi, surfm, label] = loadHCPsurf_group('L', 2);

edg = SurfStatEdg(surf);
lab = [1 2 3 4 6 7 8 9 10 11 13 14 15 16 17]; % skipping 5 and 12
