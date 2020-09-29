function [outputArg1,outputArg2] = drawRectangle(width, height, x,y,color)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%width = 100; % whatever
%height = 50; % whatever...
xCenter = x; % Wherever...
yCenter = y; % Wherever...
xLeft = xCenter - width/2;
yBottom = yCenter - height/2;
%color = [0,0,0];
rectangle('Position', [xLeft, yBottom, width, height],  'EdgeColor', color,'FaceColor', color, 'LineWidth', 1);
end

