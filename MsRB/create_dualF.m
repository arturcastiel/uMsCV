%function  [ coarseElemCenter, coarse_interface_center, coarse_strips, boundRegion, intRegion, GlobalBoundary,H, outSupport , ...
 %   refCenterInCoaseElem, dictionary,edgesCoarseDict,coarseDiricht, edges_ordering]  = 


function  [ coarseElemCenter, coarse_interface_center, coarse_strips, boundRegion, intRegion, GlobalBoundary,H, outSupport , ...
   refCenterInCoaseElem, dictionary,edgesCoarseDict,coarseDiricht, dualRegion, edges_ordering]  = create_dualF(primal_forming, primal, coarseneigh, centelem, exinterface, multiCC, splitFlag)    
%   Detailed explanation goes here
    global elem coord inedge bedge elemloc edgesOnCoarseBoundary

    %[npar,primal_forming, primal, coarseneigh, centelem, exinterface, multiCC, splitFlag]
    
    if multiCC == 2
        bcflag = true;
    else
        bcflag = false;
    end
    onflag = false;
    bcflag = true;
    mdnode = true;
    correctionweightflag = false;
    shortestflag = true;
    icbflag  = false;
    
    
    coarseElemCenter = primal.coarseElemCenter;
    npar = size(coarseElemCenter,1);    

     %% moving centers to comply with lorens bc        
    if bcflag
        for ii = 1:npar
            lfaces = setdiff(primal_forming.elfaces(ii,:),0);
            lbfaces = primal_forming.bflag(lfaces);
            
            
            if any(lbfaces)
                if sum(lbfaces) == 2
                    auxmat = primal_forming.faces(lfaces(lbfaces),:);
                    intersect_node = intersect(auxmat(1,:), auxmat(2,:));
                    point = primal_forming.coord(intersect_node,:);
                elseif sum(lbfaces) == 1
                    nodes = primal_forming.faces(lfaces(lbfaces),:);
                    p1 = primal_forming.coord(nodes(1),:);
                    p2 = primal_forming.coord(nodes(2),:);
                    point = (p1 + p2)*0.5;
                end
                
                refA = elemloc(bedge(:,3)) == ii;               
                if sum(refA) ~= 0
                    ref = false(size(elemloc));
                    ref(bedge(refA,3)) = true;
                else
                    ref = elemloc == ii;
                end
                ref(ref) = minDis(point , centelem(ref,1:2));
                coarseElemCenter(ii) = find(ref);
           
            
            elseif (sum(primal_forming.bnodes(setdiff(primal_forming.elem(ii,:),0))) == 1) && onflag
                    nodes = setdiff(primal_forming.elem(ii,:),0);
                    ref = primal_forming.bnodes(nodes);
                    point = primal_forming.coord(nodes(ref),:);
                    refA = elemloc(bedge(:,3)) == ii;               
                    if sum(refA) ~= 0
                        ref = false(size(elemloc));
                        ref(bedge(refA,3)) = true;
                    else
                        ref = elemloc == ii;
                    end
                    ref(ref) = minDis(point , centelem(ref,1:2));
                    coarseElemCenter(ii) = find(ref);
            end
        end
    end
    
    %% finding edge that represents the interface between two coarse volumes
    n = size(primal_forming.nface,1);
    coarse_interface_center = zeros(n,1);
    for ii = 1:n
        c1 = primal_forming.faces(ii,1);
        c2 = primal_forming.faces(ii,2);
        fc1 = primal_forming.bnodes(c1);
        fc2 = primal_forming.bnodes(c2);
        cflag = primal_forming.nface(ii,2) ~= 0;
        left = primal_forming.nface(ii,1);
        right = primal_forming.nface(ii,2);
        if right == 0
            refB =  elemloc(bedge(:,3)) == left;
            p1 = coord(bedge(refB,1),1:2);
            p2 = coord(bedge(refB,2),1:2);
            edge_center_point = (p1+p2)*.5;
            cpoint = primal_forming.coord(primal_forming.faces(ii,:),:);
            cpoint = mean(cpoint,1);
            edge_ref =  minDis(cpoint, edge_center_point);
            refB(refB == true) = edge_ref;
            coarse_interface_center(ii) = find(refB);         
        else
            refB = all(ismember(elemloc(inedge(edgesOnCoarseBoundary,3:4)),primal_forming.nface(ii,:)),2);
            p1 = coord(inedge(edgesOnCoarseBoundary(refB),1),1:2);
            p2 = coord(inedge(edgesOnCoarseBoundary(refB),2),1:2);
            edge_center_point = (p1+p2)*.5;            
            if (fc1 | fc2) & mdnode
                if fc1
                    cpoint = primal_forming.coord(c1,:);
                else
                    cpoint = primal_forming.coord(c2,:);
                end
            else
                cpoint = primal_forming.coord(primal_forming.faces(ii,:),:);
                cpoint = mean(cpoint,1);
            end
            edge_ref =  minDis(cpoint, edge_center_point);
            refB(refB == true) = edge_ref;
            coarse_interface_center(ii) = edgesOnCoarseBoundary(refB);
        end
    end
    
    bcoarse = coarse_interface_center(primal_forming.bflag);
    coarse_interface_center = coarse_interface_center(~primal_forming.bflag);
    
    %% creating the coarse path connecting two coarse volumes
    coarse_element_bridge = elemloc(inedge(coarse_interface_center,3:4));
    coarse_element_target = inedge(coarse_interface_center,3:4);
    
%% creat coarse path center-nodes boundary for those coarse cells located on the boundaries of the physical domain
boundaryTarget = [];
if ~bcflag
    boundaryTarget = bedge(bcoarse,3);
end

%% creating the distances between adjacencies elements    

    dist = ones(size(inedge,1),1);
    if ~shortestflag
        P1 = centelem(inedge(:, 3),1:2);
        P2 = centelem(inedge(:, 4),1:2);
        dist = vecnorm(P2-P1, 2, 2) ;
    end
   % come back
    %dist_weight = dist;    

% %% increasing weights for strips on the edgesOnCoarseBoundary
    if icbflag 
        vol_cb = unique(inedge(edgesOnCoarseBoundary,3:4));
        bad_ref = any(ismember((inedge(:,3:4)), vol_cb),2);        
        dist(bad_ref) = 3000*dist(bad_ref);
    end
 %% finding the shortest path
    coarse_strips = cell(size(coarse_element_bridge,1),1);
    half_strips = cell(size(coarse_element_bridge,1),2);    
    if ~correctionweightflag
        correction_constant = 1;
    else
        correction_constant = 0.001;
    end
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
                % come back
                ldist = dist;
                easy_dist = all(ismember( inedge(:,3:4), find(lineCross(1:size(elem,1), p1,p2))),2);               
                ldist(easy_dist) = ldist(easy_dist) * correction_constant;   
%                 cnflag1 = sum(coarseneigh(coarse_element_bridge(jj,1),npar+1));
%                 cnflag2 = sum(coarseneigh(coarse_element_bridge(jj,2),npar+1));
%                 if (splitFlag ~=1) && (cnflag1 == 1) && (cnflag2 == 1)
%                     easy_dist = all(ismember(inedge(:,3:4), find(any(ismember(elem(:,1:4), unique(bedge(:,1:2))),2))),2);
%                     ldist(easy_dist) = ldist(easy_dist) * 0.00001;  
%                 end
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
    
    
    %%find shortest path for boundaries
    n = size(boundaryTarget,1);
    bcoarse_strips = cell(n,1);
    if ~bcflag
        for ii = 1:n
                ccell = elemloc(boundaryTarget(ii));
                center = coarseElemCenter(ccell);
                target = boundaryTarget(ii);
                p1 = centelem(center,1:2);
                p2 = centelem(target,1:2);
                easy_dist = all(ismember( inedge(:,3:4), find(lineCross(1:size(elem,1), p1,p2))),2);
                ldist = dist;
                ldist(easy_dist) = ldist(easy_dist) * 0.001;
                lref = (elemloc == ccell);
                transvec = find(lref);
                transback = zeros(size(elem,1),1);
                transback(transvec) = 1:size(transvec,1);
                edge_ref = all(ismember( inedge(:,3:4), transvec),2);
                auxmat = transback(inedge(edge_ref,3:4));
                G = graph(auxmat(:,1), auxmat(:,2), ldist(edge_ref));
                path = shortestpath(G, transback(center), transback(target));
                bcoarse_strips{ii} = transvec(path)';
        end
    end
   

 
  %% finding boundary regions
  boundRegion = cell(npar,1);
  
  for ii =1:npar
%       nodes = setdiff(setdiff(primal_forming.elem(ii,:),0),ii);
%       neigh_vol = setdiff(find(any(ismember(primal_forming.elem,nodes),2)),0);
      nodes = setdiff(primal_forming.elem(ii,:),0);
      neigh_vol = setdiff(find(any(ismember(primal_forming.elem,nodes),2)),ii);
      
      ref1 = all(ismember(coarse_element_bridge, neigh_vol),2);
      ref2 = any(coarse_element_bridge == ii,2);
      ref = (ref1 & ~ref2);
      boundRegion{ii} = unique(horzcat(coarse_strips{ref}));
  end
 
 
  
 if ~bcflag
     
    for ii = 1:npar
      nodes = setdiff(primal_forming.elem(ii,:),0);
      neigh_vol = setdiff(setdiff(find(any(ismember(primal_forming.elem,nodes),2)),0),ii);       
      refB = ismember(elemloc(boundaryTarget),neigh_vol);
      
      ref = primal_forming.bflag;
      ref(ref == true) = refB;
      nodeC = primal_forming.faces(ref,:);
      
      refC = any(ismember(nodeC, nodes),2);
      refB(refB == true) = refC;
      
      
      boundRegion{ii} = setdiff(union(boundRegion{ii},   horzcat(bcoarse_strips{refB})' ),coarseElemCenter);      
      boundRegion{ii} = reshape(boundRegion{ii},1,[]);
    end
 end
%  
for ii = 1:npar
    boundRegion{ii} = setdiff(boundRegion{ii}, coarseElemCenter);
end
GlobalBoundary = false(size(elem,1),1);
ref = unique(horzcat(boundRegion{:}));
GlobalBoundary(ref) = true;
%  
%  
 %% Creating interaction Region
 
 intRegion = cell(npar,1);
 
 polyReg = cell(npar,1);
 cReg = cell(npar,1);
 for ii = 1:npar
     %nodes = setdiff(setdiff(primal_forming.elem(ii,:),0),ii);
     nodes = setdiff(primal_forming.elem(ii,:),0);   
     all_el = setdiff(find(any(ismember(primal_forming.elem,nodes),2)),0);
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
    outSupport(:,jj) = ~ismember(wholeSet,intRegion{jj})';
end
 refGlobal2Local
 
 
 %% creating dual

dualRegion = sparse(size(elem,1),npar);
for ii = 1:npar
   dualRegion(H{ii},ii) = 1; 
   dualRegion(boundRegion{ii},ii) = 2;    
end


    
     %% creating edges ordering
 
 edges_ordering = zeros(size(elemloc));
%  for ii = 1:size(coarse_strips,1)
%      edges_ordering(setdiff(boundRegion{ii}, coarseElemCenter)) = ii;
%  end

 
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
