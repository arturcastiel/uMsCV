function [forming_primal, primal,elemloc1,elemloc2] = primalDefine(flag)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    global mesh nx ny coarsemesh elemloc npar
    forming_primal = [];
    primal = [];
    dual = [];
    if flag == 0
        elemloc = gelemloc;
    elseif flag == 1
        elemloc = partitionMesh(mesh, nx, ny);
        npar = max(elemloc);
    elseif flag == 2
        filename = coarsemesh;
        [forming_primal, primal] = gridpartition(filename);
        elemloc = primal.elemloc;
        elemloc1 = elemloc;
        elemloc = graphGrowthA(elemloc);
        elemloc2 = elemloc;
       elemloc = graphIntegrity(elemloc);
       elemloc = graphIntegrity2(elemloc);
       elemloc = graphIntegrityC(elemloc);
        npar = max(elemloc);


    end
    npar = max(elemloc);

end

