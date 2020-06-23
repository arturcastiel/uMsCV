function [centelem, elem, elemarea, esurn1, esurn2, inedge, bedge, elemloc,  coarseelem,r] = reordering(centelem, elem, elemarea, esurn1, esurn2, inedge, bedge, elemloc,  coarseelem, npar,coord)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
cmat = sparse(size(elem,1),size(elem,1));
cmat(sub2ind(size(cmat),inedge(:,3),inedge(:,4))) = true;
cmat(sub2ind(size(cmat),inedge(:,4),inedge(:,3))) = true;
nCon = [];
for node =  1:size(coord,1)
    point = esurn2(node)+1;
    length = esurn2(node+1) - esurn2(node);
    element = esurn1(point :(point+length-1));
    nCon = [nCon; combnk(element,2)];
end
nCon = unique(nCon, 'rows');
cmat(sub2ind(size(cmat),nCon(:,1),nCon(:,2))) = true;
cmat(sub2ind(size(cmat),nCon(:,2),nCon(:,1))) = true;

r = symrcm(cmat)';
% r = symrcm(cmat)';
% [~, rb] = sort(r);
%% reoordering
centelem = centelem(r,:); 
elem = elem(r,:); 
elemarea = elemarea(r,:);
elemloc = elemloc(r,:);

esurn1 = r(esurn1);
bedge(:,3) = r(bedge(:,3));
inedge(:,3) = r(inedge(:,3));
inedge(:,4) = r(inedge(:,4));

%wells(:,1) = r(wells(:,1));
for ii = 1:npar
   coarseelem{ii} =  r(coarseelem{ii})';
end
end

