function[surf, surfi, surfm] = DoLSDsurf_individual(hemi)

dir = '/scr/murg2/HCP_new/HCP_Q1-Q6_GroupAvg_Related440_Unrelated100_v1/Q1-Q6_R440.';
surf_gii = gifti([dir hemi '.midthickness.32k_fs_LR.surf.gii']);
surf.coord = surf_gii.vertices'; surf.tri = surf_gii.faces;
surf_gii = gifti([dir hemi '.very_inflated.32k_fs_LR.surf.gii']);
surfi.coord = surf_gii.vertices'; surfi.tri = surf_gii.faces;
surfm.tri = surf.tri; surfm.coord = (surfi.coord .* 0.8 + surf.coord .* 0.2);

