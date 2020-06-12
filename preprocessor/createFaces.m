function [F, fElem, bfaces] = createFaces(elem)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
elem = elem(:, 1:4);
% F4 = [];
% F1 = elem(:,1:2);
% F2 = elem(:,2:3);
ref_trian = elem(:,end) == 0;

F_trian = [elem(ref_trian,1:2); elem(ref_trian,2:3); [elem(ref_trian,3) elem(ref_trian,1)]];
F_quad =  [elem(~ref_trian,1:2); elem(~ref_trian,2:3); elem(~ref_trian,3:4); [elem(~ref_trian,4) elem(~ref_trian,1)]];
% 
% 
% 
% if tFlag
%     F3 = [elem(:,3) elem(:,1)];
% else
%     F3 = elem(:,3:4);
%     F4 = [elem(:,4) elem(:,1)];
% end
% F = [F1; F2; F3; F4];
F = [F_trian; F_quad];
%F = sort([F1; F2; F3; F4],2);
F = sort(F,2);
F = unique(F,'rows');
flag = zeros(size(elem,1),1);
fElem = zeros(size(elem));
for index = 1:size(F,1)    
    ref = sum(ismember(elem,F(index,:)),2) == 2;
    flag(ref) = flag(ref) + 1;
    for ii = find(ref)'        
        fElem(ii, flag(ii)) = index;
    end
end

a = setdiff(unique(fElem),0);
tmpmat = fElem(:);
tmpmat = tmpmat(tmpmat~=0);
tmp = [a,histc(tmpmat,a)];
bfaces =  (tmp(:,2) == 1);
end

