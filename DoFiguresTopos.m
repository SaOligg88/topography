function[surfm, vert] = DoFiguresTopos()

addpath('./utils');
addpath('./utils/cbrewer');
addpath('./utils/surfstat');

%% variables: 
hemi = 'L';
maskRSN = [1 2 3 4]; 

%% Load data:
surf_gii = gifti(['/scr/murg2/HCP_new/HCP_Q1-Q6_GroupAvg_Related440_Unrelated100_v1/Q1-Q6_R440.' hemi '.midthickness.32k_fs_LR.surf.gii']);
surf.coord = surf_gii.vertices'; surf.tri = surf_gii.faces;
surf_gii = gifti(['/scr/murg2/HCP_new/HCP_Q1-Q6_GroupAvg_Related440_Unrelated100_v1/Q1-Q6_R440.' hemi '.very_inflated.32k_fs_LR.surf.gii']);
surfi.coord = surf_gii.vertices'; surfi.tri = surf_gii.faces;
surfm.tri = surf.tri; surfm.coord = (surfi.coord .* 0.8 + surf.coord .* 0.2);
% Include masks:
	lab = gifti(['RSN_' hemi '.gii']);
lab = lab.cdata(1:32492,maskRSN); % This is where the map gets set (there are 4)

vert = zeros([length(surf.coord) length(maskRSN)]);

for j = 1:length(maskRSN)
    vals = unique(lab(:,j));
    for i = 1:length(vals)   
        vert(find(lab(:,j) == vals(i)), j) = i-1;
    end

    % remove noise areas:
    %for i = 1:length(noiseVals{j})
    %    vert(find(vert(:,j) == noiseVals{j}(i))) = 0;
    %end
%     colors = [0,0,0; cbrewer('qual','Set1',length(vals-1))];
%     h = figure; set(h, 'Position', [10 20 1200 900])
%     SurfStatView(vert(:,j), surfm, []); SurfStatColormap(colors); saveas(h, ['RSN.' num2str(j) '.png']);
end
% 
% %% lines of symmetry:
% cortex = find(data.dist(10015,:));
% distStd = zeros([32492, length(unique(vert(:,2)))-1]);
% for i = 1:length(unique(vert(:,2)))-1
%    distStd(cortex,i) = std(data.dist(cortex,find(vert(:,2) == i)),0,2);
% end
% figure;
% for i = 1:17
%     subplot(4,5,i); hist(distStd(:,i),1000); ylim([0 500]); hold on;
% end
% 
% incld = [2 3 4 6 7 8 9 13 15 16 17];
% 
% 
% noise = [5 12];
% primary = [1 10 11 14];
% vertLim = zeros([1 32492]);
% for i = [incld primary]
%     vertLim(find(vert(:,2) == i)) = i;
% end
% figure; SurfStatView(vertLim, surfm);figure; 
% 
% figure;
% for i = incld2
%     subplot(4,5,i); hist(distStd(:,i),1000); ylim([0 500]); hold on;
% end
% 
% %%
incld2 = [2 3 4 7 8 15 16];
distStdNorm = zeros([32492 1]);
for i = 1:length(incld2)
    distStdNorm(find(distStd(:,i)),i) = (distStd(find(distStd(:,i)),i) - mean(distStd(find(distStd(:,i)),i))) ./ std(distStd(find(distStd(:,i)),i));
end
% figure;
% for i = 1:length(incld2)
%     subplot(3,3,i); hist(distStdNorm(:,i),1000); ylim([0 500]); hold on;
% end
distStdNormSum = sum(distStdNorm,2);
distStdNormSumP = zeros([1 32492]);
distStdNormSumP(find(distStdNormSum)) = ...
    distStdNormSum(find(distStdNormSum)) - min(distStdNormSum);
figure; SurfStatView(distStdNormSumP, surfm);
SurfStatColormap([ 0 0 0; cbrewer('seq', 'YlGnBu', 100)]); SurfStatColLim([0 15]);

% %%
incld3 = [2 3 4 7 8 15 16 6 9 ];
distStdNorm = zeros([32492 1]);
for i = 1:length(incld3)
    distStdNorm(find(distStd(:,i)),i) = (distStd(find(distStd(:,i)),i) - mean(distStd(find(distStd(:,i)),i))) ./ std(distStd(find(distStd(:,i)),i));
end
% figure; hist(sum(distStdNorm,2), 1000);
distStdNormSum = sum(distStdNorm,2);
distStdNormSumP = zeros([1 32492]);
distStdNormSumP(find(distStdNormSum)) = ...
    distStdNormSum(find(distStdNormSum)) - min(distStdNormSum);
figure; SurfStatView(distStdNormSumP, surfm);
SurfStatColormap([ 0 0 0; cbrewer('seq', 'YlGnBu', 100)]); SurfStatColLim([0 15]);

% %%
primary = [1 10 11 14];
distStdNorm = zeros([32492 1]);
for i = 1:length(primary)
    distStdNorm(find(distStd(:,i)),i) = (distStd(find(distStd(:,i)),i) - mean(distStd(find(distStd(:,i)),i))) ./ std(distStd(find(distStd(:,i)),i));
end
% figure; hist(sum(distStdNorm,2), 1000);
distStdNormSum = sum(distStdNorm,2);
distStdNormSumP = zeros([1 32492]);
distStdNormSumP(find(distStdNormSum)) = ...
    distStdNormSum(find(distStdNormSum)) - min(distStdNormSum);
figure; SurfStatView(distStdNormSumP, surfm);
SurfStatColormap([ 0 0 0; cbrewer('seq', 'YlGnBu', 100)]); SurfStatColLim([0 15]);

% %%
% distStdNorm = zeros([32492 1]);
% for i = 1:17
%     distStdNorm(find(distStd(:,i)),i) = (distStd(find(distStd(:,i)),i) - mean(distStd(find(distStd(:,i)),i))) ./ std(distStd(find(distStd(:,i)),i));
% end
% % figure; hist(sum(distStdNorm,2), 1000);
% distStdNormSum = sum(distStdNorm,2);
% distStdNormSumP = zeros([1 32492]);
% distStdNormSumP(find(distStdNormSum)) = ...
%     distStdNormSum(find(distStdNormSum)) - min(distStdNormSum);
% figure; hist(sum(distStdNormSumP,2), 1000);
% figure; SurfStatView(distStdNormSumP, surfm);
% SurfStatColormap([ 0 0 0; cbrewer('seq', 'YlGnBu', 100)]); SurfStatColLim([0 15]);
% 
% %%
% figure; hist(sum(distStd(:,incld),2), 1000);
% figure; SurfStatView(sum(distStd(:,incld),2), surfm);
% 
% figure; hist(sum(distStd(:,primary),2), 1000);
% figure; SurfStatView(sum(distStd(:,primary),2), surfm);
% 
% figure; hist(sum(distStd(:,[incld primary]),2), 1000);
% figure; SurfStatView(sum(distStd(:,[incld primary]),2), surfm);
% 
% for i = 1:17
%     h = figure('visible','off'); SurfStatView(distStd(:,i), surfm);
%     saveas(h, ['RSN.yeo17.' num2str(i) '.distStd.png']);
%     close all
%     label = zeros([1 32492]);
%     label(find(vert(:,2) == i)) = 1;
%     h = figure('visible','off'); SurfStatView(label, surfm);
%     saveas(h, ['RSN.yeo17.' num2str(i) '.region.png']);
%     close all
% end
% 
% conn = load('/scr/murg2/HCP_new/HCP_Q1-Q6_GroupAvg_Related440_Unrelated100_v1/HCP_Q1-Q6_R468_rfMRI_groupAvg_left_corr_smoothed2_toR_nifti.mat');
% 
% cortex = find(data.dist(10015,:) == 0);
% conn(isnan(conn)) = 0;
% conn(noncortex, noncortex) = 0;
% [Yh,Ih] = sort(conn,1);
% [Yh1,Ih1]=sort(Ih,1);
% corrMath = zeros([32492 32492]);
% corrMath(find(Ih1 > (32492-(32492*.1)))) = 1;
% 
% distCorr = corrMath.* dist;

% Import dist

%% Calculate cluster-wise std where each is weighted individually
data.dist = load('/scr/litauen2/projects/nki_enhanced/dist.lh.mat');
cortex = find(cortex);%find(data.dist(10015,:));
distStd = zeros([32492 length(nonzeros(unique(vert(:,2))))]);
nets = nonzeros(unique(vert(:,2)));
for i = [1 2 3 4 6 7 8 9 10 11 13 14 15 16 17]
    for j = 1:length(nonzeros(unique(clus.label(find(ismember(clus.network, i))))))
        a = nonzeros(unique(clus.label(find(ismember(clus.network, i)))));
        meanDist(:,j) = mean(data.dist(cortex,find(clus.label == a(j))),2);
    end
    if j > 1
        % distStd(cortex,i) = mean(meanDist);
        distStd(cortex,i) = mean(abs(diff(meanDist,1,2)),2);
    else
        distStd(cortex,i) = meanDist;
    end
    clear meanDist
    disp(i);
end
