addpath(genpath('./utils'));

%% variables: 
hemi = 'L'; % or: 'R'
% maskRSN = 1; noiseVals = [5]; % RSN7
maskRSN = 2; noiseVals = [5 12]; % RSN17

%% Load data:
dist = loadDist_group(hemi);
myelin = loadMyelin_group(hemi);
[surf, surfi, surfm, lab] = loadHCPsurf_group(hemi, ,maskRSN);
edg = SurfStatEdg(surf);


slm = struct();
slm.tri = surf.tri';

vert = zeros([length(surf.coord) 1]);
vals = unique(lab);
for i = 2:length(vals)   
    vert(find(lab == vals(i))) = i-1;
end
% remove noise areas:
for i = 1:length(noiseVals)
    vert(find(vert == noiseVals(i))) = 0;
end
vert(find(myelin == 0)) = 0; % remove myelin == 0 nodes

list = nonzeros(unique(vert))';
for i = list;
    ind = find(vert == i);
    d1 = dist(ind,ind);
    m1 = myelinMatL(ind,ind);
    indTril = find(tril(ones([length(ind) length(ind)]),-1));
    d2 = d1(indTril);
    m2 = m1(indTril);
    group = ones([length(indTril) 1]);  
    
    indOther = setdiff(find(vert), ind);
    d3 = [d2; reshape(dist(ind,indOther), prod(size(dist(ind,indOther))), 1)];
    m3 = [m2; reshape(myelinMatL(ind,indOther), prod(size(myelinMatL(ind,indOther))), 1)];
    m3 = sqrt(m3);
    group = [group; (ones([prod(size(myelinMatL(ind,indOther))) 1]) .* 2)];
    
    name = ['analysis ' num2str(i)];
%      [h,a,c,s] = aoctool(d3, m3, group, 0.05,...
%          'dist','myelin',name,'off','separate means');
%      [c,m(:,:,i),h,nms] = multcompare(s);
     
     [p,t,stats] = anova1(m3, group,'off');
     [c,m,h,nms] = multcompare(stats);
     disp(i);
     disp(p);
     disp([nms num2cell(m)]);
%     disp(c);
%     disp(s.pmm);
%     disp(s.pmmcov);
%    figure; boxplot(m3, group); title(name);
%    disp([name '   ' p]);
%mR(i) = mean(m3(group == 1));
%mY(i) = mean(m3(group == 2));

end

%% Do nodewise

a = gifti([dir 'HCP_Q1-Q6_R468_rfMRI_groupAvg_MGTR_left_corr_smoothed2_toR.gii']);
cortex2 = gifti([dir 'lh.cortex.gii']);
%cortex_lhh = read_label('fsaverage5', 'lh.cortex');
%cortex = [cortex(:,1)+1];
cortex3 = find(sum(cortex2.cdata,2)~=3);
%bad = [1 2 3 4 5 6 7 8 9 10 11 12];
%[cortex ind] = setdiff(cortex3, bad);
noncortex = ones([1 length(surf.coord)]);
noncortex(cortex3) = 0;
%a.cdata(find(noncortex), find(noncortex)) = 0;
connH = a.cdata(cortex3, cortex3);
distH = dist(cortex3, cortex3);
%distH(find(noncortex), find(noncortex)) = 0;
myelinMatLH = myelinMatL(cortex3, cortex3);
%myelinMatLH(find(noncortex), find(noncortex)) = 0;

R = zeros([1 length(surf.coord)]);
P = zeros([1 length(surf.coord)]);
tic
for i = 1:length(cortex3)
    [R(cortex3(i)) P(cortex3(i))] = partialcorr(connH(:,i), myelinMatLH(:,i), distH(:,i), 'type', 'Spearman'); 
    % disp(cortex3(i));
    if i == 15000
        disp(cortex3(i));
        toc
    end
end
% figure; 
% d = 1 ./ distH(:,1000);
% d(find(isinf(d))) = 1;
% subplot(1,3,1); hist(connH(:,1000), 1000); hold on;
% subplot(1,3,2); hist(myelinMatLH(:,1000), 1000); hold on;
% subplot(1,3,3); hist(d, 1000); hold on;
figure; SurfStatView(R, surf);

% Faster method:
tic
A = connH;
B = myelinMatLH;
An=bsxfun(@minus,A,mean(A,1)); %%% zero-mean
Bn=bsxfun(@minus,B,mean(B,1)); %%% zero-mean
An=bsxfun(@times,An,1./sqrt(sum(An.^2,1))); %% L2-normalization
Bn=bsxfun(@times,Bn,1./sqrt(sum(Bn.^2,1))); %% L2-normalization
C=sum(An.*Bn,1); %% correlation
R2 = zeros([1 length(surf.coord)]);
R2(cortex3) = C;
toc
figure; SurfStatView(R2, surf);


%% take similarity to only top connections with myelin, factoring
% out distance. 
% Test against relationship to to others.



%% Distance / connectivity correlation
R1 = zeros([1 length(surf.coord)]);
P1 = zeros([1 length(surf.coord)]);
tic
for i = 1:length(cortex3)
    [R1(cortex3(i)) P1(cortex3(i))] = corr(connH(:,i), distH(:,i), 'type', 'Spearman'); 
    if i == 15000
        disp(cortex3(i));
        toc
    end
end
figure; SurfStatView(R1, surf); 

% Faster method:
tic
A = connH;
B = distH;
An=bsxfun(@minus,A,mean(A,1)); %%% zero-mean
Bn=bsxfun(@minus,B,mean(B,1)); %%% zero-mean
An=bsxfun(@times,An,1./sqrt(sum(An.^2,1))); %% L2-normalization
Bn=bsxfun(@times,Bn,1./sqrt(sum(Bn.^2,1))); %% L2-normalization
C=sum(An.*Bn,1); %% correlation
R1 = zeros([1 length(surf.coord)]);
R1(cortex3) = C;
toc
figure; SurfStatView(R1, surf);
