function [I] = graphGrowthInt(coarseElemCenter, boundRegion)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global inedge elem npar  superFolder
I = false(size(elem,1),npar);

for ii = 1:npar
    seed = coarseElemCenter(ii);
    color = zeros(size(elem,1),1);
    color(coarseElemCenter) = 2;
    color(seed) = 1;
    
    color(boundRegion{ii}) = 2;
    color_old = -1;
    while ~all(color == color_old)
       color_old = color;
       ref =(color_old == 1);
       auxmat = ismember(inedge(:,3:4),find(ref));
       refF = xor(auxmat(:,1),auxmat(:,2));
       paint = unique(inedge(refF,3:4));
       refb = color_old(paint) == 0;
       paint = paint(refb);
       color(paint) = 1;
      
%         for ii = 1:npar
%             ref1 = (elemloc(inedge(:,3)) == ii) & (elemloc(inedge(:,4)) == 0);
%             ref2 = (elemloc(inedge(:,4)) == ii) & (elemloc(inedge(:,3)) == 0);
%             ref = ref1 | ref2;
%             if sum(ref) ~= 0
%                 elemloc(unique(inedge(ref,3:4))) = ii;
%             end
%         end

    end
    
   color(coarseElemCenter) = 0;
   color(boundRegion{ii}) = 0;
   I(:,ii) = color == 1;
end

end

