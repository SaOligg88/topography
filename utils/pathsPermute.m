function [randPath] = pathsPermute(clus, n, steps, networks)

% config steps path:
%[steps, networks] = config_steps_path();
% my 18:
%[steps, networks] = config_my18_path();


uniquelabels = nonzeros(unique(clus.label));   
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
    for i = 1:length(steps.s)
        steps.l{i} = nonzeros(unique(clus.label(find(ismember(c, steps.s{i})))));
    end
    a = zeros(length(uniquelabels));
    for i = 1:length(steps.s)-1
        a(steps.l{i},steps.l{i+1}) = clus.lab(steps.l{i},steps.l{i+1});
        a(steps.l{i+1},steps.l{i}) = clus.lab(steps.l{i+1},steps.l{i});
    end
    
    [D,P] = dijk(a, steps.l{1}, steps.l{length(steps.s)});
    randPath(k) = length(find(D == length(steps.s)-1));
    disp(k);
end

% get p-value
[s,ind] = sort(randPath);
p = 1 - find(ind == length(randPath))/length(randPath);

if n ~= 0
    figure; hist(randPath(1:n),n/10); ylim([0 n/20]);  hold on
    plot([randPath(n+1) randPath(n+1)],[0 n/20],'r','LineWidth',2);
    title(['Parcellation ' num2str(length(networks)) ' -- Num Permutations: ' num2str(n) ' -- p-value: ' num2str(p)]);
end




