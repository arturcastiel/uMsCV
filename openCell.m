function [A] = openCell(fileName)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fid = fopen(fileName,'r');
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);
end

