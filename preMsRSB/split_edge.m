function [edges_split, edges_inv, internal_split] = split_edge(edges, perm_matrix)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 global dualRegion  elemloc coarseElemCenter intRegion npar boundRegion
[~, ref] = sort(coarseElemCenter);


internals = false(size(elemloc));
internals(edges) = true;
internals(coarseElemCenter) = true;
internals = ~internals;

internals_ref = (perm_matrix* internals == true);

internal_split = false(sum(internals),npar);

for ii = 1:npar
    aux = false(size(elemloc));
    aux(intRegion{ii}) = true;
    aux = (perm_matrix * aux);
    internal_split(:, ii) = aux(internals_ref); 
end

internal_split = internal_split(:, elemloc(coarseElemCenter(ref))');

edges_ref = (perm_matrix*edges == true);
edges_split = (perm_matrix * dualRegion);
edges_split = (edges_split(edges_ref,:));
edges_split = edges_split(:, elemloc(coarseElemCenter(ref))');
edges_inv = false(size(edges_split,1), size(edges_split,1));

auxmat = edges_split == 1;

for ii = 1:size(edges_split,1)
   mm  = find(auxmat(ii,:));
   for jj = mm
      edges_inv(ii, (auxmat(:,jj) == true)) = true;       
   end    
end


% for ii = 1:npar
%     
%     
% end
% edges_m = full(edges_m(bb,:)) == 1;
% edges_m = (edges_m == true); 

% for ii = 1:npar
%     auxvec = (perm_matrix * (dualRegion(:,ii) == 2)) == true;    
%     edges_split(:,ii) = auxvec(edges_ref);    
% % end
% 
%     edges_ref(edges_ref ==1) = edges_split(:,11) == 2;
% %     
%     bb = ((perm_matrix' * edges_ref) == true);
% %     1

end

