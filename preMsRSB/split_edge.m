function [edges_split, edges_inv] = split_edge(edges, perm_matrix)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 global dualRegion  elemloc coarseElemCenter
[~, ref] = sort(coarseElemCenter);


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
1
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

