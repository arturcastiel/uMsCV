function [out] = semiplan(coord_list, lpoint1, lpoint2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
R = [0 , -1; 1, 0];
center_line = 0.5*(lpoint1+lpoint2);
line_vec = lpoint2 - lpoint1;
normal = R * line_vec';
test_vec = coord_list - repelem(center_line, size(coord_list,1),1);
cnormal = repelem(normal', size(coord_list,1),1);
cross_product = sum((test_vec .* cnormal),2);
out = -1*ones(size(cross_product,1),1);
out(cross_product > 0) = 1;
end

