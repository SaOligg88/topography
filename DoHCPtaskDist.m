function[mo] = DoHCPtaskDist()

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

% Load DMN dist from yeo:
[distDMN, zonesDMN] = makeDMNdist_yeo();

%% Calc distance to tasks from DMN
percThr = 80;
% count = 1;
% for percThr = 0:10:100

for i = 1:length(b)
    disp(b{i});
end
for i = 1:length(b)
    data = tasks.(char(b{i}))(1:32492*2);
    dataAll(i,:) = data;
    dataL = data(1:32492);
    dataR = data(32493:32492*2);
    indFin = find(isfinite(dataL).*isfinite(dataR));
    ht(i) = corr(dataL(indFin),dataR(indFin));
    htsub(i,:) = dataL - dataR;
    
    vals = find(data > prctile(data, percThr));
    mo{i} = distDMN(vals);
    mom(i) = mean(distDMN(vals));
    momd(i) = median(distDMN(vals));
    %moc(i) = corr(data(isfinite(data)), distDMN(isfinite(data)), 'type', 'Kendall');
    
    valsNeg = find(data < prctile(data, percThr));
    data(valsNeg) = 0;
    htsubThr(i,:) = data(1:32492) - data(32493:32492*2);
    %disp(i);  
end
% indFin = (find(isfinite(mean(abs(htsubThr(ind(good),:))))));
% [r p ] = corrcoef(distDMNL(indFin), mean(abs(htsubThr(ind(good),indFin)))');
% rall(count) = r(1,2);
% pall(count) = p(1,2);
% disp(count);
% count = count + 1;
% end
ind = find(isfinite(mean(htsub(in,:))));
for i = 1:length(in)
    %ma(i) = mean(distDMNL(find(abs(htsubThr(in(i),:)) > 2.3)));
    %vals = find(abs(htsubThr(in(i),:)) > 3.1);
    %ma(i) = median(distDMNL(vals) .* zscore(abs(htsub(in(i),vals)))');
    %[r p ] = corrcoef(distDMNL(ind), (abs(htsubThr(in(i),ind)))');
    [r p ] = corrcoef(distDMNL(ind), (abs(htsub(in(i),ind)))');
    ra(i) = r(1,2);
    pa(i) = p(1,2);
end
ind = find(isfinite(mean(htsub(in,:))));
for i = 1:length(in)
    figure; SurfStatView(abs(htsub(in(i),:)),surfL); SurfStatColLim([-10 10]);
    saveas(gcf, ['hcp.task.' num2str(i) '.png']);
    close all
end

% for i = 1:32492
% %r = corrcoef([dataAll(:,i) dataAll(:,32492+i)]);
% r = corrcoef([dataAll(in,i) dataAll(in,32492+i)]);
% rv(i) = r(1,2);
% end

[s ind] = sort(momd);
[s indm] = sort(mom);
[s indc] = sort(moc);

c = zeros(length(b),1);
for i = 1:length(b)
    if length(findstr(char(b{ind(i)}),'neg')) ~= 0 
        c(i) = 1;
    end
    if length(findstr(char(b{ind(i)}),'wm_avg_')) ~= 0 
        c(i) = 1;
    end
    if length(findstr(char(b{ind(i)}),'_')) == 2 
        c(i) = 2;
    end
end

good = find(c == 2);
%good = find(c ~= -1);
clear p
for i = 1:length(good)
    disp(b{ind(good(i))});
    p(i) = b{ind(good(i))};
    p(i) = strrep(p(i), 'tfmri_', '');
    p(i) = strrep(p(i), '_', ' ');
end

c = zeros(length(good),1);
for i = 1:length(good)
    if size(findstr(char(b{ind(good(i))}),'motor') ~= 0 )
        c(i) = 1;
    end
    if size(findstr(char(b{ind(good(i))}),'_wm_') ~= 0 )
        c(i) = 2;
    end
        if size(findstr(char(b{ind(good(i))}),'_social_') ~= 0 )
        c(i) = 3;
        end
        if size(findstr(char(b{ind(good(i))}),'_emotion_') ~= 0 )
        c(i) = 4;
        end
        if size(findstr(char(b{ind(good(i))}),'_language_') ~= 0 )
        c(i) = 5;
        end
    if size(findstr(char(b{ind(good(i))}),'bk') ~= 0 )
        c(i) = 6;
    end
        if size(findstr(char(b{ind(good(i))}),'relation') ~= 0 )
        c(i) = 7;
        end
        if size(findstr(char(b{ind(good(i))}),'gambling') ~= 0 )
        c(i) = 8;
    end
end
d = zeros(length(good),3);
d(find(c == 1),1) = 1;
d(find(c == 0),3) = 1;
d(find(c == 2),2) = 1;
d(find(c == 6),3) = 0.5; d(find(c == 6),1) = 0.2;
d(find(c == 3),1) = 0.5; d(find(c == 3),2) = 0.5;
d(find(c == 4),2) = 0.5; d(find(c == 4),3) = 0.5;
d(find(c == 5),3) = 0.5; d(find(c == 5),1) = 0.5;
d(find(c == 7),2) = 0.5; d(find(c == 7),3) = 0.5; d(find(c == 7),1) = 0.5;
d(find(c == 8),3) = 0.15; d(find(c == 8),2) = 0.15; d(find(c == 8),1) = 0.15;

g = figure;
[h,L,MX,MED,bw]=violin(mo(ind(good)),'facecolor',d,'facealpha',0.5,'edgecolor','none', 'mc','k','medc','k.','bw', 0.5);
set(gca, 'xTick',1:length(p),'XTickLabel', p);
rotateticklabel();

%% Individual level analysis
subs = loadHCP();
taskTypes={'MOTOR', 'EMOTION', 'GAMBLING', 'LANGUAGE',  'RELATIONAL', 'SOCIAL', 'WM'};
%clear moc2 %names2 names2Ind
moc2s = zeros(491,142);
moc2p = zeros(491,142);
moc2mean = zeros(491,142);
moc2median = zeros(491,142);
for i = 1:length(subs)
    disp([subs{i} '  ' num2str(length(subs) - i + 1)]); 
    % Include calculation on individual level of DMN peaks:
    
    count = 1;
    for k = 1:length(taskTypes)
        if exist(['/a/documents/connectome/_all/' subs{i} '/MNINonLinear/Results/tfMRI_' taskTypes{k} '/tfMRI_' taskTypes{k} '_hp200_s12_level2.feat/' subs{i} '_tfMRI_' taskTypes{k} '_level2_hp200_s12.dscalar.nii'], 'file')
            cifti = ft_read_cifti(['/a/documents/connectome/_all/' subs{i} '/MNINonLinear/Results/tfMRI_' taskTypes{k} '/tfMRI_' taskTypes{k} '_hp200_s12_level2.feat/' subs{i} '_tfMRI_' taskTypes{k} '_level2_hp200_s12.dscalar.nii']);
            tasks = fieldnames(cifti);
            for j = 1:length(tasks)
                if ~isempty(findstr(char(tasks(j)),'_tfmri_')) && length(findstr(char(tasks{j}),'_')) == 6
                        data = cifti.(char(tasks{j}))(1:32492*2);
                        in = find(data > 3.1);
                        %inFin = find(isfinite(data));
                        if length(in) ~= 0
                            %moc2s(i,count) = corr(data(inFin), distDMN(inFin), 'Type', 'Spearman'); 
                            %moc2p(i,count) = corr(data(inFin), distDMN(inFin), 'Type', 'Pearson'); 
                            %moc2mean(i,count) = mean(distDMN(in));
                            moc2median(i,count) = median(distDMN(in));
                        end
                end                
                count = count + 1;
            end
        end
    end
end
valid = [2:8 36:37 50:51 64:65 78:79 92:93 114:115  120:123];
mogMed.valid = valid;
%mog2s = moc2s(find(mean(moc2s,2) ~= 0),valid);
%mog2p = moc2p(find(mean(moc2p,2) ~= 0),valid);
%mog2mean = moc2mean(find(mean(moc2mean,2) ~= 0),valid);
mog2median = moc2median(find(mean(moc2median,2) ~= 0),valid);

for i = 1:size(mog2s,2)
%    mogMed.s(i) = median(nonzeros(mog2s(:,i)));
%    mogMed.p(i) = median(nonzeros(mog2p(:,i)));
%    mogMed.mean(i) = median(nonzeros(mog2mean(:,i)));
    mogMed.median(i) = median(nonzeros(mog2median(:,i)));
%    mogMed.All.s{i} = nonzeros(mog2s(:,i));
%    mogMed.All.p{i} = nonzeros(mog2p(:,i));
%    mogMed.All.mean{i} = nonzeros(mog2mean(:,i));
    mogMed.All.median{i} = nonzeros(mog2median(:,i));
end    
%[s mogMed.ind.s] = sort(mogMed.s);
%[s mogMed.ind.p] = sort(mogMed.p);
%[s mogMed.ind.mean] = sort(mogMed.mean);
[s mogMed.ind.median] = sort(mogMed.median);
%mogMed.names = names;
% Get names:
% clear names
% i = 2
% count = 1;
% c2 = 1;
% for k = 1:length(taskTypes)
%     cifti = ft_read_cifti(['/a/documents/connectome/_all/' subs{i} '/MNINonLinear/Results/tfMRI_' taskTypes{k} '/tfMRI_' taskTypes{k} '_hp200_s12_level2.feat/' subs{i} '_tfMRI_' taskTypes{k} '_level2_hp200_s12.dscalar.nii']);
%     tasks = fieldnames(cifti);
%     for j = 1:length(tasks)
%         if findstr(char(tasks(j)),'_tfmri_') & length(findstr(char(tasks{j}),'_')) == 6                
%             nam = strrep(char(tasks{j}), ['x' subs{i} '_tfmri_'], '');
%             nam = strrep(nam, '_hp200_s12', '');
%             nam = strrep(nam, '_level2', '');
%             names{count} = nam;
%             count = count + 1;   
%             namesInd(c2) = 1;
%         end
%         c2 = c2 + 1;
%     end
% end
% Make figues:
% c = zeros(length(names),1);
% for i = 1:length(names)
%     if size(findstr(char(names{i}),'motor') ~= 0 )
%         c(i) = 1;
%     end
%     if size(findstr(char(names{i}),'wm') ~= 0 )
%         c(i) = 2;
%     end
%         if size(findstr(char(names{i}),'social') ~= 0 )
%         c(i) = 3;
%         end
%         if size(findstr(char(names{i}),'emotion') ~= 0 )
%         c(i) = 4;
%         end
%         if size(findstr(char(names{i}),'language') ~= 0 )
%         c(i) = 5;
%         end
%     if size(findstr(char(names{i}),'bk') ~= 0 )
%         c(i) = 6;
%     end
%         if size(findstr(char(names{i}),'relation') ~= 0 )
%         c(i) = 7;
%         end
%         if size(findstr(char(names{i}),'gambling') ~= 0 )
%         c(i) = 8;
%     end
% end
% 
% d = zeros(length(names),3);
% d(find(c == 1),1) = 1;
% d(find(c == 0),3) = 1;
% d(find(c == 2),2) = 1;
% d(find(c == 6),3) = 0.5; d(find(c == 6),1) = 0.2;
% d(find(c == 3),1) = 0.5; d(find(c == 3),2) = 0.5;
% d(find(c == 4),2) = 0.5; d(find(c == 4),3) = 0.5;
% d(find(c == 5),3) = 0.5; d(find(c == 5),1) = 0.5;
% d(find(c == 7),2) = 0.5; d(find(c == 7),3) = 0.5; d(find(c == 7),1) = 0.5;
% d(find(c == 8),3) = 0.15; d(find(c == 8),2) = 0.15; d(find(c == 8),1) = 0.15;
% 
% clear p
% for i = 1:length(names)
%     b(i) = names(i);
%     p(i) = strrep(b(i), '_', ' ');
% end
% mogMed.p = p;
% mogMed.d = d;

ind = mogMed.ind.median;
mog = mogMed.All.median(ind);
d = mogMed.d(ind,:);
p = mogMed.p(ind);
g = figure;
[h,L,MX,MED,bw]=violin(mog,'facecolor',d,'facealpha',0.5,'edgecolor','none', 'mc','k','medc','k.');
set(gca, 'xTick',1:length(mogMed.names),'XTickLabel', p);
rotateticklabel();
%title('Spearman');
%title('Pearson');
%title('Mean');
title('Median');
%%
% [s indc] = sort(moc);

covs = csvread('/scr/litauen1/solar/hcp/HCP_pheno_01.csv', 1, 0);
count = 1;
clear covariates
for i = 1:length(subs)
    disp(i);
    indcov = find(covs(:,1) == str2num(subs{i}));
    if indcov
        covariates(count,:) = [covs(indcov,1:3) moc(i,:)];
        count = count + 1;
    end
end
filename = '/scr/litauen1/solar/hcp/HCP_pheno_tasks_motor.csv';
savePhenoCSV(filename, covariates, ['ID','Sex','Age',names]);


%% Load surface area

subs = loadHCPall();
for i = 1:length(subs)
    disp(subs{i});
    unix(['wb_command -surface-vertex-areas /a/documents/connectome/_all/' subs{i} '/MNINonLinear/' subs{i} '.L.midthickness.164k_fs_LR.surf.gii area/' subs{i} '.mid.L.gii']);
    saL = gifti(['area/' subs{i} '.mid.L.gii']);
    surfAreaL(i) = sum(saL.cdata);
    unix(['wb_command -surface-vertex-areas /a/documents/connectome/_all/' subs{i} '/MNINonLinear/' subs{i} '.R.midthickness.164k_fs_LR.surf.gii area/' subs{i} '.mid.R.gii']);
    saR = gifti(['area/' subs{i} '.mid.L.gii']);
    surfAreaR(i) = sum(saR.cdata);
end

% Include covariates for total surface area
surfAreaTotal = surfAreaL + surfSareaR;

subs = loadHCPall();
for i = 1:length(subs)
    disp(subs{i});
    saL = gifti(['area/' subs{i} '.mid.L.gii']);
    saLall(:,i) = saL.cdata;
%    save_avw(['area/' subs{i} '.mid.L.nii'],saL.cdata,'f',[1 1 1 1]);
%    saR = gifti(['area/' subs{i} '.mid.L.gii']);
%    surfAreaR(i) = sum(saR.cdata);
end

surf_gii = gifti('/a/documents/connectome/_all/100307/MNINonLinear/100307.L.midthickness.164k_fs_LR.surf.gii');
surfm.coord = surf_gii.vertices'; surfm.tri = surf_gii.faces;
surf_gii = gifti('/a/documents/connectome/_all/100307/MNINonLinear/100307.L.very_inflated.164k_fs_LR.surf.gii');
surfi.coord = surf_gii.vertices'; surfi.tri = surf_gii.faces;
surfL.tri = surfm.tri; surfL.coord = (surfi.coord .* 0.8 + surfm.coord .* 0.2);

surf_gii = gifti('/a/documents/connectome/_all/100307/MNINonLinear/100307.R.midthickness.164k_fs_LR.surf.gii');
surfm.coord = surf_gii.vertices'; surfm.tri = surf_gii.faces;
surf_gii = gifti('/a/documents/connectome/_all/100307/MNINonLinear/100307.R.very_inflated.164k_fs_LR.surf.gii');
surfi.coord = surf_gii.vertices'; surfi.tri = surf_gii.faces;
surfR.tri = surfm.tri; surfR.coord = (surfi.coord .* 0.8 + surfm.coord .* 0.2);

surf.tri = [surfL.tri; surfR.tri+length(surfL.coord)];
surf.coord = [surfL.coord surfR.coord];


gifti('/a/documents/connectome/_all/100307/MNINonLinear/100307.L.midthickness.164k_fs_LR.surf.gii area/100307.mid.L.gii');


% Create masks gifti file
163842 - (55*55*55)
saLall_full = [saLall; zeros([2533,510])];
saLall_3d = reshape(saLall_full,[55 55 55 510]);
save_avw('/scr/litauen1/solar/hcp/mid.L.nii',saLall_3d,'f',[1 1 1 1]);
