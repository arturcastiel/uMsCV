function [A] = convertInput(A, folder)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global osMode
if folder == 0
    folder = pwd;
end
if strcmp(osMode,'windows')
    A{28}  = [folder, '\results'];
    A{286}  = [folder, '\grids'];
    A{294}  = [folder, '\geo'];
elseif strcmp(osMode,'windows')
    A{28}  = [folder, '/results'];
    A{286}  = [folder, '/grids'];
    A{294}  = [folder, '/geo'];   
end

end

