function [bOp,sOp ] = velStruc()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global inedge bedge elemloc npar

mat = false(size(inedge,1), npar);
for ii = 1:size(inedge,1)
    mat(ii,elemloc(inedge(ii,3))) = true;
    mat(ii,elemloc(inedge(ii,4))) = true;
end

out = 1;

mat2 = mat((sum(mat,2) == 2), :);

sOp = mat2;
bOp = mat;
end

