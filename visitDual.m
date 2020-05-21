elemdual = zeros(size(elem,1),1);
elemdual(GlobalBoundary) = 10;
elemdual(coarseElemCenter) = 20;
postprocessorTMS(full(elemdual), full(elemdual),0,superFolder,'elemdual')

lop = sparse(size(elem,1), max(edges_ordering));

for ii = 1:max(edges_ordering)
    ref = (edges_ordering == ii);
    lop(ref,ii) = ii;    
end
postprocessorOP(lop,1,  superFolder, 'lop')


