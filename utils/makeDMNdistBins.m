% make binned maps
a = zeros(length(distDMN),1);
for i = 85:-5:0
    a(find(distDMN <= i)) = i/5;
end
DoFlatten(a(1:32492)', 'hcp.dist.5.DMN.L', 'label_col_left', '32');
DoFlatten(a(32493:32492*2)', 'hcp.dist.5.DMN.R', 'label_col_right', '32');

