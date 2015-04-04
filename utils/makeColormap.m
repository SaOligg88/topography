function[colors] = makeColormap(filename, clus)

newCols = csvread(filename);
colors = [0 0 0];
for i = 1:length(unique(clus.edgeNet))
    addCol = newCols(clus.edgeNet == i,:);    
    colors = [colors; addCol(1,:)];
end