function [outputArg1,outputArg2] = impFlux(TransF, F, pc, neumman_flux)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end


function [perm_coarse]= perm_coarse()
    global elemloc
    [~,r] = sort(elemloc);
     ord = @(i,j,A)( (j- ones(size(j)))*size(A,1) +i);
    matrix = sparse(size(r,1),size(r,1));
    ind = ord([1:size(r,1)]', r, matrix);
    matrix(ind) = true;
end