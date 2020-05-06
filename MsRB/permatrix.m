function [perm_vec, perm_matrix, bmat] = permatrix()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global elemloc npar
n = size(elemloc,1);
[~, r] = sort(elemloc);
ord = @(i,j,A)( (j- ones(size(j)))*size(A,1) +i);
matrix = sparse(size(elemloc,1),size(elemloc,1));
ind = ord([1:size(elemloc,1)]', r, matrix);
matrix(ind) = true;
perm_matrix = matrix;
pelemloc = perm_matrix * elemloc;
block_matrix = false(size(elemloc,1),size(elemloc,1));
for ii = 1:npar
    ref = pelemloc == ii;
    block_matrix(ref,ref) = true; 
end
bmat = sparse(block_matrix);
perm_vec = r;
end

