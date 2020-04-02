function [ coarseElemCenter, coarse_interface_center, coarse_strips, boundRegion, intRegion, GlobalBoundary,H, outSupport , ...
    refCenterInCoaseElem, dictionary,edgesCoarseDict,coarseDiricht, coarse_element_bridge,  coarse_element_target] = alpreMsRB(npar,coarseelem, coarseneigh, centelem, exinterface, multiCC, splitFlag)    
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%finding neighbors 
global intCoord elem coord inedge bedge esurn1 esurn2 elemloc edgesOnCoarseBoundary intinterface polyReg cReg


bcflag = multiCC;
coarseElemCenter = zeros(npar,1);
%% recalculate centers
    for ii = 1:npar
        ref = elemloc == ii;
        center_point = weiszfeld(centelem(ref,1:2));
        ref(ref) = minDis(center_point , centelem(ref,1:2));
        coarseElemCenter(ii) = find(ref);
    end
    
%% finding edge that represents the interface between two coarse volumes
    ref = intCoord(:, 3) == 2;
    points = intCoord(ref, 1:2);
    coarse_interface_center = zeros(size(points, 1), 1);
    p1 = coord(inedge(edgesOnCoarseBoundary,1),1:2);
    p2 = coord(inedge(edgesOnCoarseBoundary,1),1:2);
    edge_center_point = (p1+p2)./2;
    for ii = 1:size(points, 1)
       edge_ref =  minDis(points(ii,:), edge_center_point);
       coarse_interface_center(ii) = edgesOnCoarseBoundary(edge_ref);
    end
    
%% creating the coarse path connecting two coarse volumes
    coarse_element_bridge = elemloc(inedge(coarse_interface_center,3:4));
    coarse_element_target = inedge(coarse_interface_center,3:4);
    
    
%% moving centers to comply with lorens bc    
    bcflag = true;
    
    if bcflag == 1
        % moving center of coarse volumes
        for ii = 1:npar
           int_sum = sum(coarseneigh(ii,end-3:end));
           if int_sum == 1
%                ref = (intCoord(:,3) == 3) & (intCoord(:,4) == ii);
%                point = intCoord(ref,1:2);
               cref = elemloc == ii;
               point = centerinterface(exinterface{ii} ,bedge, coord,elemloc);

               cref(cref) = minDis(point, centelem(cref,1:2));
               coarseElemCenter(ii) = find(cref);               
           
     
           elseif int_sum ==2
               ref = (intCoord(:,3) == 4) & (intCoord(:,4) == ii);
               point = intCoord(ref,1:2);
               cref = elemloc == ii;
               cref(cref) = minDis(point, centelem(cref,1:2));
               coarseElemCenter(ii) = find(cref);               
           end  
        end
        % moving center of interfces
        all_boundary_nodes = unique(bedge(:,1:2));
        bcoarse = any(coarseneigh(:,end-3:end),2);
        bEdges = any(ismember(inedge(edgesOnCoarseBoundary,1:2),all_boundary_nodes),2);
        for ii = 1:size(coarse_element_bridge,1)
           flag = all(bcoarse(coarse_element_bridge(ii,:)));
           if flag == 1
              left =  coarse_element_bridge(ii,1);
              right =  coarse_element_bridge(ii,2);
              %ref1 = all(ismember(elemloc(inedge(edgesOnCoarseBoundary,3:4)),[left,right]),2);
              aux = elemloc(inedge(edgesOnCoarseBoundary,3:4));
              ref1 = (aux(:,1) == left) & (aux(:,2) == right);
              ref2 = (aux(:,2) == left) & (aux(:,1) == right);
              ref = (ref1 | ref2) &  bEdges;
              if sum(ref) ~= 0
                  if all(elemloc(inedge(edgesOnCoarseBoundary(ref),3:4)) == coarse_element_bridge(ii,:)) == 1
                      coarse_element_target(ii,:) = inedge(edgesOnCoarseBoundary(ref),3:4);
                  else
                      coarse_element_target(ii,1) = inedge(edgesOnCoarseBoundary(ref),4);
                      coarse_element_target(ii,2) = inedge(edgesOnCoarseBoundary(ref),3);
                  end
              end
              
           end
         
            
        end
        
        
    end    

%% creat coarse path center-nodes boundary for those coarse cells located on the boundaries of the physical domain
    boundaryTarget = zeros(npar,1);
    numBoundary = sum(coarseneigh(:,end-3:end), 2);
    bedge_vol = unique(bedge(:,3));
    for ii = 1:npar
        if numBoundary(ii) == 1
            edges_ref = exinterface{ii};
            edge_nodes = bedge(edges_ref,1:2);
            center_bedges =  0.5*(coord(edge_nodes(:,1),1:2) + coord(edge_nodes(:,2),1:2));
            average_point = mean(center_bedges);
            target_edge = minDis(average_point, center_bedges);
            boundaryTarget(ii) = bedge(edges_ref(target_edge),3);     
        elseif numBoundary(ii) == 2
            ref =  (intCoord(:,3) == 4) & (intCoord(:,4) == ii);
            center = intCoord(ref,1:2);
            center_bedges = centelem(bedge_vol,1:2);            
            target_vol= minDis(center, center_bedges);
            boundaryTarget(ii,:) = bedge_vol(target_vol);  
        end
    end

    

%% creating the distances between adjacencies elements    
    P1 = centelem(inedge(:, 3),1:2);
    P2 = centelem(inedge(:, 4),1:2);
    dist = vecnorm(P2-P1, 2, 2) ;


%% increasing weights for strips on the edgesOnCoarseBoundary
    vol_cb = unique(inedge(edgesOnCoarseBoundary,3:4));
    bad_ref = any(ismember((inedge(:,3:4)), vol_cb),2);
    %
    
%% creating the distances
    dist(:) = 1;
    dist_weight = dist;
    
    dist_weight(bad_ref) = 3000*dist(bad_ref);
%     G = graph(inedge(:,3),inedge(:,4), dist_weight);
    %G = graph(inedge(:,3),inedge(:,4));

    
    
    
%% finding the shortest path
    coarse_strips = cell(size(coarse_element_bridge,1),1);
    half_strips = cell(size(coarse_element_bridge,1),2);
    
    for ii = 1:npar
        center = coarseElemCenter(ii);
        p1 = centelem(center,1:2);
        lref = (elemloc == ii);
        transvec = find(lref);
        transback = zeros(size(elem,1),1);
        transback(transvec) = 1:size(transvec,1);
        edge_ref = all(ismember( inedge(:,3:4), transvec),2);
        auxmat = transback(inedge(edge_ref,3:4));
        for jj = 1:size(coarse_strips,1)
            flag = ismember(coarse_element_bridge(jj,:), ii);
            if (flag(1) == 1)  | (flag(2) == 1)
                flag2 = (elemloc(coarse_element_target(jj,:))' == ii);
                target = coarse_element_target(jj,flag2);
                p2 = centelem(target,1:2);                
                easy_dist = all(ismember( inedge(:,3:4), find(lineCross(1:size(elem,1), p1,p2))),2);
                ldist = dist;
                ldist(easy_dist) = ldist(easy_dist) * 0.001;   
                cnflag1 = sum(coarseneigh(coarse_element_bridge(jj,1),npar+1));
                cnflag2 = sum(coarseneigh(coarse_element_bridge(jj,2),npar+1));
                if (splitFlag ~=1) && (cnflag1 == 1) && (cnflag2 == 1)
                    easy_dist = all(ismember(inedge(:,3:4), find(any(ismember(elem(:,1:4), unique(bedge(:,1:2))),2))),2);
                    ldist(easy_dist) = ldist(easy_dist) * 0.00001;  
                end
                G = graph(auxmat(:,1), auxmat(:,2), ldist(edge_ref));
                path = shortestpath(G, transback(center), transback(target));
                path = transvec(path);
                half_strips{jj, flag} = path;
            end         
        end
    end
    
    for ii = 1:size(coarse_strips,1)
        coarse_strips{ii} = vertcat(half_strips{ii,:})';
    end
    

%% finding the shortest path for boundary elements

bcoarse_strips = cell(npar,1);
if bcflag == 0
    for ii = 1:npar
        if numBoundary(ii) > 0
            center = coarseElemCenter(ii);
            target = boundaryTarget(ii);
            p1 = centelem(center,1:2);
            p2 = centelem(target,1:2);  
            
            easy_dist = all(ismember( inedge(:,3:4), find(lineCross(1:size(elem,1), p1,p2))),2);
            ldist = dist;
            ldist(easy_dist) = ldist(easy_dist) * 0.001;       
            lref = (elemloc == ii);
            transvec = find(lref);
            transback = zeros(size(elem,1),1);
            transback(transvec) = 1:size(transvec,1);
            edge_ref = all(ismember( inedge(:,3:4), transvec),2);
            auxmat = transback(inedge(edge_ref,3:4));
            G = graph(auxmat(:,1), auxmat(:,2), ldist(edge_ref));
            path = shortestpath(G, transback(center), transback(target));
            bcoarse_strips{ii}= transvec(path)';
        end
    end
end

 %% finding boundary regions
    boundRegion = cell(npar,1);
    coarse_brige = sort(coarse_element_bridge,2);
    for ii=1:npar
        celsur = find(coarseneigh(ii,1:npar));
        cnpar = [];
        for jj = celsur
            celsur_local = intersect(celsur, find(coarseneigh(jj,1:npar)));
            
            celsur_local_boundary =  coarseneigh(celsur_local, end-3:end);
            cf =  any(coarseneigh(ii, end-3:end));
            all_par = sort([jj * ones(size(celsur_local))' , celsur_local'], 2);
            cnpar = unique([cnpar; all_par], 'rows');
            for ff = 1:size(celsur_local_boundary,1)
                if (sum(celsur_local_boundary(ff,:)) > 0) & (cf == 1)
                     boundRegion{ii} =  unique([boundRegion{ii}, bcoarse_strips{celsur_local(ff),1}]);
                end
            end            
        end            
        for jj = 1:size(cnpar,1)
            ref = ismember(coarse_brige, cnpar(jj,:), 'rows');
            boundRegion{ii} = unique([boundRegion{ii}, coarse_strips{ref}]); 
        end
    end
    
    

%% Creating interaction Region

intRegion = cell(npar,1);

polyReg = cell(npar,1);
cReg = cell(npar,1);


for ii = 1:npar
    all_el = union(find(coarseneigh(ii, 1:npar)), ii);
    cand_el = ismember(elemloc,all_el);
    cand_el(boundRegion{ii}) = false;
    ref = all(ismember(inedge(:,3:4), find(cand_el)),2);
    auxmat = inedge(ref,3:4);
    transvec = unique(auxmat);
    transback = zeros(size(elem,1),1);
    transback(transvec) = 1:size(transvec,1);
    auxmat = transback(auxmat);
    G = graph(auxmat(:,1), auxmat(:,2));
    con_comp = conncomp(G);
    ref_el = (con_comp == con_comp(transback(coarseElemCenter(ii))));
    intRegion{ii} = transvec((ref_el));
end


for ii = 1:npar
    boundRegion{ii} = setdiff(boundRegion{ii}, coarseElemCenter);
end

GlobalBoundary = false(size(elem,1),1);
ref = unique(horzcat(boundRegion{:}));
GlobalBoundary(ref) = true;

for ii = 1:npar
    intRegion{ii} = setdiff(intRegion{ii}, coarseElemCenter);
end

H = cell(npar,1);
for ii =1:npar
    H{ii} = intersect(intRegion{ii}, find(GlobalBoundary)) ;
end

%% generate out of support
outSupport = false(size(elem,1),npar);
size(outSupport)
wholeSet = 1:size(elem,1);
for jj = 1 : npar
    %ismember(wholeSet,intRegion{jj})
    outSupport(:,jj) = ~ismember(wholeSet,intRegion{jj})';
    
    %outSupport{jj} = logical(setdiff(wholeSet, intRegion{jj}));
    %outSupport{jj} = setdiff(outSupport{jj},boundRegion{jj}); 
end
refGlobal2Local

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