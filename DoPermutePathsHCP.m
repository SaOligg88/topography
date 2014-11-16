function [randPath] = DoPermutePathsHCP(clus, n)

yeo17.s{1} = [10 14 11];
yeo17.s{2} = [6 3 16]; 
yeo17.s{3} = [4 15];
yeo17.s{4} = [7];
yeo17.s{5} = [8];
yeo17.s{6} = [2];

uniquelabels = nonzeros(unique(clus.label));   
networks = [ 2 3 4 6 7 8 10 11 14 15 16];

for k = 1:n+1
    clus.lab = clus.edge;
    c = clus.network;
    [b ind] = sort(rand(length(networks),1));
    labR = networks(ind);
    if (n ~= 0 & k ~= n+1);
        for i = 1:length(networks)
            c(find(clus.network == networks(i))) = labR(i);
        end
    end
    for i = 1:6
        yeo17.l{i} = nonzeros(unique(clus.label(find(ismember(c, yeo17.s{i})))));
    end
    a = zeros(length(uniquelabels));
    for i = 1:5
        a(yeo17.l{i},yeo17.l{i+1}) = clus.lab(yeo17.l{i},yeo17.l{i+1});
        a(yeo17.l{i+1},yeo17.l{i}) = clus.lab(yeo17.l{i+1},yeo17.l{i});
    end
    
    [D,P] = dijk(a, yeo17.l{1}, yeo17.l{6});
    randPath(k) = length(find(D == 5));
    disp(k);
end

if n ~= 0
    figure; hist(randPath(1:n),n/100); ylim([0 n/20]);  hold on
    plot([randPath(n+1) randPath(n+1)],[0 n/20],'r','LineWidth',2);
    title(['Parcellation ' 17 ' -- Num Permutations: ' num2str(n)]);
end

