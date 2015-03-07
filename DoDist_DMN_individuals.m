function [peaks clus clusid surf dm1 dm2] = DoIndividDMNdist(sub)

sub = 120111; 
sub = num2str(sub);

surf_gii = gifti(['/a/documents/connectome/_all/' sub '/MNINonLinear/fsaverage_LR32k/' sub '.L.inflated.32k_fs_LR.surf.gii']);
surf.coord = surf_gii.vertices'; surf.tri = surf_gii.faces;
%surf_gii = gifti(['/scr/murg2/HCP_new/HCP_Q1-Q6_GroupAvg_Related440_Unrelated100_v1/Q1-Q6_R440.' hemi '.very_inflated.32k_fs_LR.surf.gii']);
%surfi.coord = surf_gii.vertices'; surfi.tri = surf_gii.faces;
%surfm.tri = surf.tri; surfm.coord = (surfi.coord .* 0.8 + surf.coord .* 0.2);

%% RSN
hemi = 'L';
maskRSN = [1 2];
lab = gifti(['../yeoTopo/RSN_' hemi '.gii']);
lab = lab.cdata(1:32492,[1 2]); % This is where the map gets set (there are 4)
vert = zeros([length(surf.coord) length(maskRSN)]);
for j = 1:length(maskRSN)
    vals = unique(lab(:,j));
    for i = 1:length(vals)   
        vert(find(lab(:,j) == vals(i)), j) = i-1;
    end
end

%% Checking group
tasks = ft_read_cifti('/scr/murg2/HCP_new/HCP_Q1-Q6_GroupAvg_Related440_Unrelated100_v1/HCP_Q1-Q6_R440_tfMRI_ALLTASKS_level3_zstat1_hp200_s2.dscalar.nii');
a = fieldnames(tasks);
b = {};
count = 1;
for i = 1:length(a)
    if findstr(char(a(i)),'tfmri')
        b{count} = a(i);
        count = count + 1;
    end
end
b = b';
for i = 1:length(b)
    data = tasks.(char(b{i}))(1:32492*2);
    dataAll(i,:) = data;
    dataL = data(1:32492);
    dataR = data(32493:32492*2);
    dataL(isinf(dataL)) = 0;
    dataL(isnan(dataL)) = 0;
    vals(i) = mean(dataL(find(vert(:,1) == 3)));
end
for i = 1:7
    dataL = mean(dataAll(find(vals > i),1:32492));
     dataL(isinf(dataL)) = 0;
    dataL(isnan(dataL)) = 0;
    valsG(i) = mean(dataL(find(vert(:,1) == 3)));
end
b{find(vals > 5)}; % leaves 5 tasks
figure; SurfStatView(mean(dataAll(find(vals > 5), 1:32492)),surf);

%% Func data for individual
Tasks={'WM' 'RELATIONAL' 'SOCIAL' 'MOTOR' 'EMOTION' 'LANGUAGE' 'GAMBLING'};
clear data
for i = 1:length(Tasks)
    task = Tasks{i};
    file = ['/a/documents/connectome/_all/' sub '/MNINonLinear/Results/tfMRI_' task '/tfMRI_' task '_hp200_s12_level2.feat/' sub '_tfMRI_' task '_level2_hp200_s12.dscalar.nii'];
    data{i} = ft_read_cifti(file);
end

Cons = {['x' sub '_tfmri_wm_level2_neg_body_hp200_s12'],...
['x' sub '_tfmri_wm_level2_neg_face_hp200_s12'],...
['x' sub '_tfmri_wm_level2_neg_place_hp200_s12'],...
['x' sub '_tfmri_wm_level2_neg_tool_hp200_s12'],...
['x' sub '_tfmri_wm_level2_neg_0bk_hp200_s12'],...
['x' sub '_tfmri_relational_level2_neg_match_hp200_s12'],...
['x' sub '_tfmri_relational_level2_neg_rel_hp200_s12'],...
['x' sub '_tfmri_social_level2_neg_tom_hp200_s12'],...
['x' sub '_tfmri_social_level2_neg_random_hp200_s12']};
match = [1 1 1 1 1 2 2 3 3];

Cons = {['x' sub '_tfmri_motor_level2_neg_cue_hp200_s12'],...
['x' sub '_tfmri_motor_level2_avg_cue_hp200_s12'],...
['x' sub '_tfmri_wm_level2_avg_face_hp200_s12'],...
['x' sub '_tfmri_relational_level2_neg_match_hp200_s12'],...
['x' sub '_tfmri_social_level2_neg_random_hp200_s12']};
match = [4 4 1 2 3 ];

Cons = {['x' sub '_tfmri_motor_level2_neg_cue_hp200_s12'],...
['x' sub '_tfmri_motor_level2_avg_cue_hp200_s12']};
match = [4 4 ];

Cons = {['x' sub '_tfmri_motor_level2_neg_cue_hp200_s12'],...
['x' sub '_tfmri_motor_level2_avg_cue_hp200_s12'],...
['x' sub '_tfmri_wm_level2_avg_face_hp200_s12'],...
...%['x' sub '_tfmri_relational_level2_neg_match_hp200_s12'],...
['x' sub '_tfmri_social_level2_neg_random_hp200_s12'],...
['x' sub '_tfmri_wm_level2_avg_tool_hp200_s12'],...
['x' sub '_tfmri_language_level2_neg_math_hp200_s12'],...
['x' sub '_tfmri_language_level2_neg_story_hp200_s12']};%,...
match = [4 4 1 3 1 6 6];

Cons = {['x' sub '_tfmri_motor_level2_neg_cue_hp200_s12'],...
['x' sub '_tfmri_motor_level2_avg_cue_hp200_s12'],...
['x' sub '_tfmri_wm_level2_neg_2bk_hp200_s12'],...
['x' sub '_tfmri_wm_level2_neg_0bk_hp200_s12'],...
['x' sub '_tfmri_social_level2_neg_random_hp200_s12'],...
['x' sub '_tfmri_wm_level2_avg_tool_hp200_s12'],...
['x' sub '_tfmri_language_level2_neg_math_hp200_s12'],...
['x' sub '_tfmri_language_level2_neg_story_hp200_s12']};%,...
match = [4 4 1 1 3 1 6 6];

% Currently using from mean of group within DMN mask
Cons = {['x' sub '_tfmri_gambling_level2_neg_punish_hp200_s12'],...
['x' sub '_tfmri_gambling_level2_neg_reward_hp200_s12'],...
['x' sub '_tfmri_language_level2_story_math_hp200_s12'],...
['x' sub '_tfmri_language_level2_neg_math_hp200_s12'],...
['x' sub '_tfmri_emotion_level2_neg_shapes_hp200_s12']};%,...
match = [7 7 6 6 5];

clear d;
for i = 1:length(Cons)
    disp(i);
    myexp = ['data{' num2str(match(i)) '}.' Cons{i}];
    d(:,i) = evalin('caller',myexp);
end  
d = d(1:32492,:);
d(find(isnan(d)))= 0;
dm = mean(d,2);
dm1 = dm; dm2 = dm; dm0 = dm;
dm0(find(dm0 < 0)) = 0; % Thresh > 3.1
%dm0(find(vert(:,1) ~= 3)) = 0;
%dm1(find(vert(:,1) ~= 3)) = 0;
%dm1(find(dm1 < 0)) = 0;
%dm2(find(vert(:,2) ~= 2)) = 0;
%dm2(find(dm2 < 0)) = 0;

%ds = SurfStatSmooth(dm', surf, 1);

%% threshold:
edg = SurfStatEdg(surf);
slm = struct();
slm.tri = surf.tri';
slm.t = dm0';
[peak,clus,clusid] = SurfStatPeakClus(slm,(vert(:,1) ==3),0.001, ones(1,length(surf.coord)), edg);
% former mask: ones([length(surf.coord) 1])
cl = clusid;
cl(clusid > 7) = 0;
cl(clusid > 7 & clusid ~= 12 & clusid ~= 8 | clusid == 2 | clusid == 1) = 0;

ids = unique(nonzeros(cl));
clear peaks
for i = 1:length(ids)%1:7%length(clus.clusid)
    ind = find(peak.clusid == ids(i));
    [m im] = max(peak.t(ind));
    peaks(i) = peak.vertid(ind(im));
end
a = zeros(32492,1);
a(peaks) = 1;
figure; SurfStatView(a, surf);


clear peaks
for i = 1:7%length(clus.clusid)
    ind = find(peak.clusid == i);
    [m im] = max(peak.t(ind));
    peaks(i) = peak.vertid(ind(im));
end
a = zeros(32492,1);
a(peaks) = 1;
figure; SurfStatView(a, surf);

surf_gii = gifti(['/a/documents/connectome/_all/' sub '/MNINonLinear/fsaverage_LR32k/' sub '.L.very_inflated.32k_fs_LR.surf.gii']);
surfo.coord = surf_gii.vertices'; surfo.tri = surf_gii.faces;
surf_gii = gifti(['/a/documents/connectome/_all/' sub '/MNINonLinear/' sub '.L.very_inflated.164k_fs_LR.surf.gii']);
surfi.coord = surf_gii.vertices'; surfi.tri = surf_gii.faces;

c = SurfStatInd2Coord(peaks, surfo);
inds = SurfStatCoord2Ind(c', surfi);
    
% %% Conn data:
% dim = 32492;
% filename = ['/scr/murg2/HCP_Q3_glyphsets_left-only/' sub '/rfMRI_REST_left_corr_avg.gii.data'];
% file = fopen(filename,'r');
% conn = fread(file,[dim,dim],'float32');
% conn(isnan(conn)) = 0;
% 
% %% Perhaps take peak meandistConn instead to get DMN?
 
[distances,zones] = DoExactGeodesicClust164(inds, sub);
% 
%surf_gii = gifti(['/a/documents/connectome/_all/' sub '/MNINonLinear/' sub '.L.very_inflated.164k_fs_LR.surf.gii']);
%surfi.coord = surf_gii.vertices'; surfi.tri = surf_gii.faces; 
surf_gii = gifti(['/a/documents/connectome/_all/' sub '/MNINonLinear/' sub '.L.midthickness.164k_fs_LR.surf.gii']);
surfm.coord = surf_gii.vertices'; surfm.tri = surf_gii.faces;
figure; SurfStatView(zones, surfm);

figure; SurfStatView(distances, surfi);
figure; SurfStatView(zones, surfi);



