function  paintRec(ref, color)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%width = 100; % whatever
%height = 50; % whatever...
 maxl =  max(find(ref));
 minl =  min(find(ref));
 len = maxl - minl;
 xc = 0.5*(maxl + minl);

width = len;
height = len;
xCenter = xc; % Wherever...
yCenter = xc; % Wherever...
xLeft = xCenter - width/2;
yBottom = yCenter - height/2;
%color = [0,0,0];
rectangle('Position', [xLeft, yBottom, width, height],  'EdgeColor', color,'FaceColor', color, 'LineWidth', 1);
end

