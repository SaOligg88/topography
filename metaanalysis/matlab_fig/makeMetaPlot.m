
clear rVals tree D leafOrder
names = h5read('neurosynth_meta.mat','/names');
inputNorm = h5read('neurosynth_meta.mat','/insputNorm');

for i = 1:prod(size(inputNorm))
    rVals(i) = str2num(inputNorm{i});
end
rVals = reshape(rVals,[size(inputNorm)]);
% Threshold:
rVals(find(rVals < 3.1 & rVals > -3.1)) = 0;

tree = linkage(rVals,'average');
D = pdist(rVals);
leafOrder = optimalleaforder(tree,D);

%cluster
%Z = linkage(D);
%dendrogram(Z);

h = figure('Position', [100, 50, 900, 900]); 
imagesc(rVals(leafOrder,:));
set(gca, 'yTick',1:length(rVals),'YTickLabel', names(leafOrder));

if size(rVals,2) == 7
    xaxis = {'medial parietal','lateral parietal','ventrolateral pfc',...
        'medial temporal','ventromedial pfc','dorsolateral pfc','lateral temporal'};
else
    for i = 1:size(rVals,2)
        xaxis{i} = [num2str((i-1)*5) '-' num2str(((i-1)*5)+5)];
    end
end
set(gca, 'xTick',1:size(rVals,2),'XTickLabel', xaxis)
rotateXLabels(gca, 90);

colormap(flip(cbrewer('div','RdBu', 50))); 
%set(gca, 'CLim', [-0.15, 0.15]);
set(gca, 'CLim', [-10, 10]);
%set(gca, 'CLim', [-1*max(max(rVals)), max(max(rVals))]);
hcb = colorbar;
set(get(hcb,'YLabel'),'String','(Pearson''s r-value)','rotation',270)
%set(hcb,'YTick',[-0.15 0 0.15])
%set(hcb,'YTickLabel',{'-0.15','0','0.15'},...
%    'FontSize',10)

set(gca,'box','off','color','none')
set(h,'color','w');

%t = text(50, -5, {'X-axis' 'label'}, 'FontSize', 14);
%set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
%'Rotation',45);

xlabel('Distance from default-mode network peaks (mm)','FontSize',16,'FontWeight','bold');
ylabel('NeuroSynth Feature Term','FontSize',16,'FontWeight','bold');

print -depsc 'Fig.dist.pattern.dist.eps'
pause(10); close all

%% Make brain pics
a = zeros(32492*2,1);
a(find(distDMN < 32.5)) = 1;
a(find(distDMN >= 32.5 & distDMN < 50)) = 2;
a(find(distDMN >= 50)) = 3;

rsnL = gifti('../../Dropbox/misc/yeoTopo/RSN_L.gii');
rsnR = gifti('../../Dropbox/misc/yeoTopo/RSN_R.gii');
rsn = [rsnL.cdata(:,1); rsnR.cdata(:,1)];
a(find(rsn == 37)) = 0;

surf = GetHCPSurfs;
figure; SurfStatView(a,surf);
SurfStatColormap(([0 0 0; 0 139 159; 255 146 0; 47 165 80])./255);

DoFlatten(a(1:32492)', 'hcp.meta.dist.L', 'label_col_left', '32');
DoFlatten(a(32493:32492*2)', 'hcp.meta.dist.R', 'label_col_right', '32');
