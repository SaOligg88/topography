function[myelin] = loadMyelin_group(hemi)

dir = ['/scr/murg2/HCP_new/HCP_Q1-Q6_GroupAvg_Related440_Unrelated100_v1/'];
myelin = gifti([dir 'Q1-Q6_R440.MyelinMap_BC.32k_fs_LR.' hemi '.gii']);
myelin = myelin.cdata;
