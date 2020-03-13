function [out] = lineCross(elemId, lpoint1, lpoint2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global centelem
R = [0 , -1; 1, 0];
d = 10^10;
ref = semilineCross(elemId,lpoint1, lpoint2);
normal = (R*(lpoint2 - lpoint1)')';
normal = d * normal;
center = 0.5*(lpoint1 + lpoint2);
A = lpoint2 + normal;
B = lpoint1 + normal;
C = lpoint1 - normal;
D = lpoint2 - normal;
poly = [A;B;C;D];

flag = inpolygon(center(1), center(2), poly(:,1), poly(:,2));
last = inpolygon(centelem(elemId(ref),1) , centelem(elemId(ref),2) , poly(:,1), poly(:,2));
if flag == 1
    ref(ref) = last;
else
    ref(~ref) = ~last; 
end
out = ref;
end

