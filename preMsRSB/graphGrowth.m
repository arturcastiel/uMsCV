function [elemloc] = graphGrowth(elemloc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global inedge elem
npar = max(elemloc);
for ii = 1:npar
    ref1 = (elemloc(inedge(:,3)) == ii) & (elemloc(inedge(:,4)) == 0);
    ref2 = (elemloc(inedge(:,4)) == ii) & (elemloc(inedge(:,3)) == 0);
    ref = ref1 | ref2;
    while sum(ref) ~= 0
        elemloc(unique(inedge(ref,3:4))) = ii;
        ref1 = (elemloc(inedge(:,3)) == ii) & (elemloc(inedge(:,4)) == 0);
        ref2 = (elemloc(inedge(:,4)) == ii) & (elemloc(inedge(:,3)) == 0);
        ref = ref1 | ref2;   
    end
end
end

