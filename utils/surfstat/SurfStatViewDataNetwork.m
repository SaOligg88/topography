function [ a, cb ] = SurfStatViewDataNetwork( data, surf )
clf;
[ld wd] = size(data);
set(gcf, 'Position',[0 0 500 wd*200])
colormap(spectral(256));
background='white';

%%
t=size(surf.tri,1);
v=size(surf.coord,2);
tmax=max(surf.tri,[],2);
tmin=min(surf.tri,[],2);
for i=1:t-1
    tmax(i+1)=max(tmax(i+1),tmax(i));
    tmin(t-i)=min(tmin(t-i),tmin(t-i+1));
end
cut=min([find((tmin(2:t)-tmax(1:t-1))==1) t]);
cuv=tmax(cut);
tl=1:cut;
tr=(cut+1):t;
vl=1:cuv;
vr=(cuv+1):v;

h=1/(wd+1);
w=0.48;
h1 = 1;
%%
cortex = find(sum(data,2) ~= 0);
noncortex = setdiff([1:32492], cortex);

for i = 1:2:wd*2
    d = data(vl,(i/2)+0.5);
    d(noncortex) = (min(d(cortex))-1);
    a(i)=axes('position',[0.01 h1-(((i/2)+0.5)/wd) w h]);
    trisurf(surf.tri(tl,:),surf.coord(1,vl),surf.coord(2,vl),surf.coord(3,vl),...
        double(d),'EdgeColor','none');
    view(-90,0); 
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;
    lighting phong; material dull; shading flat;

    a(i+1)=axes('position',[0.505 h1-(((i/2)+0.5)/wd) w h]);
    trisurf(surf.tri(tl,:),surf.coord(1,vl),surf.coord(2,vl),surf.coord(3,vl),...
        double(d),'EdgeColor','none');
    view(90,0);
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;
    lighting phong; material dull; shading flat;
end
%%    

for i=1:wd%length(a)
    c = nonzeros(data(vl,i));
    [~,idx] = sort(c);
   % [sy,idy] = sort(idx);
    maxD=c(idx(floor(length(idx)-(length(idx)*0.001))));
    minD=c(idx(floor(length(idx)-(length(idx)*0.99))))-5;
    %maxD=c(idx(length(c)));
    %minD=c(idx(1))-5;
    
    %clim=[min(data(vl,i))-1,max(data(vl,i))];
    clim = [minD, maxD];
    set(a((i*2)-1),'CLim',clim); set(a(i*2),'CLim',clim);
    set(a((i*2)-1),'Tag',['SurfStatView ' num2str((i*2)-1) ' ' num2str(0)]);
    set(a((i*2)),'Tag',['SurfStatView ' num2str((i*2)) ' ' num2str(0)]);
end

% cb=colorbar('location','South');
% set(cb,'Position',[0.35 0.085 0.3 0.03]);
% set(cb,'XAxisLocation','bottom');
% h=get(cb,'Title');
% set(h,'String',title);

whitebg(gcf,background);
set(gcf,'Color',background,'InvertHardcopy','off');

dcm_obj=datacursormode(gcf);
set(dcm_obj,'UpdateFcn',@SurfStatDataCursor,'DisplayStyle','window');

%set(gcf,'PaperPosition',[0.25 2.5 6 4.5]);
set(gcf,'PaperPosition',[0 0 5 30]);

return
end
