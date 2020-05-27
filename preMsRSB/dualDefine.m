function [ coarseElemCenter, coarse_interface_center, coarse_strips, boundRegion, intRegion, GlobalBoundary,H, outSupport , ...
    refCenterInCoaseElem, dictionary,edgesCoarseDict,coarseDiricht, dualRegion, edges_ordering] =  dualDefine(dType, primal_forming, primal, npar, coarseneigh, centelem, exinterface, multiCC, splitFlag)
    out = 0;
    if dType == 1
        % Olav
        1
    elseif dType == 2
        % Artur 1
        2
    elseif dType == 3
        %modArt
        
     [coarseElemCenter, coarse_interface_center, coarse_strips, boundRegion, intRegion, GlobalBoundary, H, outSupport, refCenterInCoaseElem, ...
         dictionary,edgesCoarseDict,coarseDiricht, dualRegion, edges_ordering] = create_dualF(primal_forming, primal, coarseneigh, centelem, exinterface, multiCC, splitFlag)  
%         
%      [coarseElemCenter, coarse_interface_center, coarse_strips, boundRegion, intRegion, GlobalBoundary, H, outSupport, refCenterInCoaseElem, ...
%          dictionary,edgesCoarseDict,coarseDiricht, edges_ordering] = smartDual(npar,primal_forming, primal, coarseneigh, centelem, exinterface, multiCC, splitFlag)  

    end
end


function [out] = minDis(center, points)
    mcenter = repelem(center, size(points,1),1);
    dists = vecnorm(mcenter - points, 2,2);
    ref = dists == min(dists);
    if sum(ref) == 1
        out = ref;
    else        
        pref = false(size(ref,1),1);
        target = find(ref);
        pref(target(1)) = true;
        out = pref;
    end
end