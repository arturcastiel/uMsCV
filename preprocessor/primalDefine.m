function [elemloc] = primalDefine(flag, gelemloc)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    global mesh nx ny coarsemesh
    if flag == 0
        elemloc = gelemloc;
    elseif flag == 1
        elemloc = partitionMesh(mesh, nx, ny);
    elseif flag == 2
        filename = coarsemesh;
        elemloc = gridpartition(filename);
    end
end

