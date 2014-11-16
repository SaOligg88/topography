function[dist] = DoDistMatHCP(sub)

addpath('/scr/litauen1/Dropbox/misc/surfstat/');
addpath('/scr/litauen1/Dropbox/misc/topology');

% import surface:
filename = ['/scr/dattel2/' num2str(sub) '/MNINonLinear/fsaverage_LR32k/' num2str(sub) ...
    '.L.midthickness.32k_fs_LR.surf.gii'];
surf_gii = gifti(filename);
surf.coord = surf_gii.vertices';
surf.tri = surf_gii.faces;
clear surf_gii;

aparc = gifti(['/scr/dattel2/' num2str(sub) '/MNINonLinear/fsaverage_LR32k/' num2str(sub) ...
    '.L.aparc.a2009s.32k_fs_LR.label.gii']);
aparc = aparc.cdata;

% Make noncortex mask
noncortex = zeros([1 length(surf.coord)]);
noncortex(find(aparc == 0)) = 1;
cortex = zeros([1 length(surf.coord)]);
cortex(find(noncortex == 0)) = 1;

dist = zeros([length(surf.coord) length(surf.coord)]);
cortVals = find(cortex);
[surf max_degree] = surfGetNeighbors(surf);
noncortexIDs = find(noncortex);
surf.nbr(:,noncortexIDs) = 0;
for i = 1:length(noncortexIDs)
    surf.nbr(find(surf.nbr == noncortexIDs(i))) = 0;
end

for p = 1:length(cortVals)
    tic;
    %    mask = zeros([1 length(surf.coord)]);
    %    mask(cortVals(i)) = 1;
    %    dist(cortVals(i),:) = surfGeoDist(surf, mask, noncortex);
    %    disp([num2str(i) ' is done running of ' num2str(length(cortVals))]);
    %end

    origin = cortVals(p);
    geoDepth = zeros(1,length(surf.coord));%dist(cortVals(p),:); %

    % indexPreCal(i)= 0 if the computation is not done yet on i-th vertex in
    % previous iteration, indexPreCal(i)= 1 if computation is done. 
    indexPreCal    = zeros(1,length(surf.coord));
    indexPreCal(find(noncortex)) = 100000;

    % start the firefront from these vertices 
    current_set    = origin;
    % index the initial set
    indexPreCal(current_set)=1;
    %indexPreCal(find(geoDepth))=1;
    
    k = 1;
    k_1=0;
    tmp=1;
    % iteration #
    l=2;
    current_cal=surf.nbr(:,current_set)';

    % stop iteration when no computation is done in the previous step, 
    % i.e., when k=k_1
    %disp('start iteration')
    while k > k_1
        k_1=k;

        for i=1:prod(size(current_cal))%length(current_cal)
            if (current_cal(i) ~= 0)
                t=current_cal(i);
                if (indexPreCal(t)==0)
                    tmp=1000000;
                    % if there is a neighbor with pre-computation

                    for j=1:max_degree

                       % for all existing neigbors (sometimes there is 0 vertex nbr) 
                       if surf.nbr(j,t) ~= 0
                           % compute depth if a current vertex has neighbor 
                           % vertices on which depth is computed in the 
                           % previous iteration 

                           if ( indexPreCal(surf.nbr(j,t))>0 && ...
                                indexPreCal(surf.nbr(j,t))<l )
                               nbrCoord=surf.coord(:,surf.nbr(j,t));
                               curCoord=surf.coord(:,t);
                               dist_cur=sqrt(sum((nbrCoord - curCoord).^2,1));
                               geoDepth(t)=min(dist_cur+geoDepth(surf.nbr(j,t)), tmp);
                               tmp=geoDepth(t);
                               indexPreCal(t)=l;
                               % inculde neiborhood of the current vertex in
                               % next iteration 
                               current_set=[current_set surf.nbr(:,t)'];

                           end
                       end   
                    end
                   k=k+1;
                end 
            end
        end

        % avoid double computation
        current_cal=unique(current_set);
        % remove 0 neighbor
        a=current_cal>0;
        current_cal=current_cal(a);
        current_set=[];
        l=l+1;
        %disp(l)
    end
    dist(cortVals(p),:) = geoDepth;
    %dist(p,:) = geoDepth;
    %dist(:,cortVals(p)) = geoDepth;
    disp([num2str(p) ' is done running of ' num2str(length(cortVals))]);
    toc;
end


% % load conn
dim = 32492;
filename = ['/scr/murg2/HCP_Q3_glyphsets_left-only/' num2str(sub) '/rfMRI_REST_left_corr_avg.gii.data'];
file = fopen(filename,'r');
conn = fread(file,[dim,dim],'float32');
conn(isnan(conn)) = 0;
conn(find(noncortex), find(noncortex)) = 0;
[Yh,Ih] = sort(conn,1);
[Yh1,Ih1]=sort(Ih,1);
corrMath = zeros([32492 32492]);
corrMath(find(Ih1 > (32492-(32492*.1)))) = 1;

distCorr = corrMath.* dist;
meandistCorr = zeros([1 32492]);
cortexN = find(cortex);
for i = 1:length(cortexN); meandistCorr(cortexN(i)) = mean(nonzeros(distCorr(:,cortexN(i)))); end
figure; SurfStatView(meandistCorr, surf);
% % save
save(['/scr/kaiser2/corticalDist/' num2str(sub) '.L.dist.32k.mat'], '-v7.3','dist');
% save(['/scr/kaiser2/corticalDist/' num2str(sub) '.L.meandistCorr.32k.mat'], '-v7.3','meandistCorr');
% 



    