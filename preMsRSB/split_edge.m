function [edges_split] = split_edge(edges)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 global dualRegion perm_matrix elemloc coarseElemCenter
[~, ref] = sort(coarseElemCenter);
 
edges_ref = (perm_matrix*edges == true);
edges_split = (perm_matrix * dualRegion);
edges_split = (edges_split(edges_ref,:));
edges_split = edges_split(:, elemloc(coarseElemCenter(ref))');
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

