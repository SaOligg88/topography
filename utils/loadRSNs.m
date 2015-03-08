function[lab] = loadRSNs(hemi, maskRSN)

lab = gifti(['./data/RSN_' hemi '.gii']);
if hemi == 'L'
	nodes = 1:32492;
elseif hemi == 'R'
	nodes = 32493:32492*2;
end
lab = lab.cdata(nodes,maskRSN); % This is where the map gets set (there are 4)
