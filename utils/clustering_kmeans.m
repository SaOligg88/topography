function[clusters] = clustering_kmeans(k, x, y, surf)

clusters = kmeans([x y],k);

a = zeros(length(surf.coord),1);
for i = 1:k
    a(ismember(label, find(clusters == i))) = i;
end
figure; SurfStatView(a, surf);