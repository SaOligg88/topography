function[geoDistL32, geoDistPrimaryL32] = DoAllHCP()

List=dir('/a/documents/connectome/_all/');
for i = 1:length(List)-2
    ListSub(i,:) = List(i+2).name;
end

for i = 1:length(ListSub)
    disp(['Running Subject: ' num2str(ListSub(i,:))]);
    [geoDistL, geoDistPrimaryL] = DoDistRootsHCP(ListSub(i,:));
    geoDistL32(:,i) = ConvertTo32k(ListSub(i,:), 'geoDistL', geoDistL);
    geoDistPrimaryL32(:,i) = ConvertTo32k(ListSub(i,:), 'geoDistPrimaryL', geoDistPrimaryL);
    
%    geoDist32(:,i) = ConvertTo32k(ListSub(i,:), 'geoDist', geoDist);
%    geoDistPrimary32(:,i) = ConvertTo32k(ListSub(i,:), 'geoDistPrimary', geoDistPrimary);
end

% figure; SurfStatView(geoDist, surf);

