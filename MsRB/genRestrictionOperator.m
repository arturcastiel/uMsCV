function [ OR ] = genRestrictionOperator(flag)
%genProlongationOperator Generates Restriction Operator
%   INPUT:
%   nelem = number of elements in the mesh
%   npar = number of partions
%   OUTPUT:
%   OR = Restriction operator
    global coarseelem elemloc elem npar coarseElemCenter
    % flag 1 - olav, flag 0 ams
    flag = true;
    nelem = size(elem,1);
    OR = sparse(zeros([npar , nelem]));       
    refC = elemloc(sort(coarseElemCenter));
    for ii = 1:npar
        if flag
            ij = refC(ii); 
        else
            ij = ii;
        end
        OR(ii,coarseelem{ij}) = true;    
    end
    %OR = sparse(OR);
end

