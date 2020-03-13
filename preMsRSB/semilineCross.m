function [out] = semilineCross(elemId, lpoint1, lpoint2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global elem coord
global_elem = elem(elemId,1:4);
if all(global_elem(:,4) == 0)
    global_elem = global_elem(:,1:3);
end
test_nodes = setdiff(unique(global_elem),0);
transvec = zeros(max(test_nodes),1);
transvec(test_nodes) = 1:size(test_nodes,1);
local_elem = transvec(global_elem);
ref = semiplan(coord(test_nodes, 1:2), lpoint1 , lpoint2);
res = sum(ref(local_elem),2);
ref_end = (res>= 3) | (res <= -3);
out = ~ref_end;
end

