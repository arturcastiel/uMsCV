
Q = [centelem(:,1:2), kmap(elem(:,5),5)];
rep = 20;
nclust = 8;
deg = 3;



for ii = 1:rep
    
    F = clusterdata(Q,'Linkage','centroid','MaxClust',30);
    Q(:,end) = F;
    G = F;
    F = clusterdata(Q,'Linkage','centroid','MaxClust',deg);
    Q(:,end) = F;
end


% 
% 
% for ii = 30
%     G = clusterdata(Q,'Linkage','centroid','MaxClust',ii);
%     Q(:,end) = G;
% end
postprocessorTMS(full(F ),full(G),0,superFolder,'kmedia');