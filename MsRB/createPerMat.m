function [matrix, r] = createPerMat(id_classify)
    [~, r] = sort(id_classify);
    ord = @(i,j,A)( (j- ones(size(j)))*size(A,1) +i);
    matrix = sparse(size(id_classify,1),size(id_classify,1));
    ind = ord([1:size(id_classify,1)]', r, matrix);
    matrix(ind) = true;
end


%% 
% 
%