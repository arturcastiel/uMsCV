function [forming_primal, primal] = gridpartition(filename)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
    [forming_primal] = create_forming_primal(filename);
    [primal] = create_primal(forming_primal);
    %aa = createfdual(forming_primal);
%     F = forming_primal.faces;
%     tcoord = forming_primal.coord;
    % F = coarse.faces;
    % tcoord = coarse.coord;
%     for ii = 1:size(F)
%        drawLine(F(ii,1), F(ii,2), tcoord)
%     end
%     %[forming_primal, primal] = createfdual(forming_primal);

end

function [rprimal] = create_primal(primal)
    rprimal.elemloc = inVol(primal.coord,primal.elem);
    n = size(primal.elem,1);
    rprimal.coarseElemCenter = zeros(n,1);
    for ii = 1:n
        center =   primal.centelem(ii, :);
        ref = (rprimal.elemloc == ii);
        rprimal.coarseElemCenter(ii) = find(minDis(center, ref));        
    end
%     rprimal.coarse_element_bridge = elemloc(inedge(coarse_interface_center,3:4));
%     coarse_element_target = inedge(coarse_interface_center,3:4);

end



function [coarse] = create_forming_primal(filename)
    [tcoord, tnode] = getcoord(filename);
    tcoord = tcoord(:,1:2);
    telem = getelemf(filename, tnode);
    % trian_flag = all(telem(:,4) == 0);
    [F, tfaces,bflag] = createFaces(telem);   
    npar = size(telem,1);
    [tcentelem]  = findcentelem(tcoord, telem);
    [nface] = create_faceneigh(tfaces);
    bnodes = ismember(1:size(tcoord,1), unique(F(bflag,:)))';
    
    %elemloc = inVol(primal.coord,primal.elem);
    coarse = struct('coord', tcoord, 'elem',telem, 'faces', F, 'elfaces', tfaces, 'bflag', bflag, 'npar', npar, 'centelem', tcentelem, 'nface', nface, 'bnodes',bnodes); %, 'elemloc', elemloc);    
end


function [dual] = createfdual(coarse)
    mid_faces = 0.5.*(coarse.coord(coarse.faces(:,1),:) + coarse.coord(coarse.faces(:,2),:));
    refN = 1:size(coarse.coord,1);
    refC = [1:size(coarse.centelem,1)] + size(refN,2) ;
    refM = [1:size(mid_faces,1)] + (size(refN,2) + size(refC,2));
    nodes  = [coarse.coord; coarse.centelem; mid_faces];
    elem = []';
    for ii = 1:size(coarse.coord,1)
        tvolumes = find(any(coarse.elem == ii,2))';
        tfaces = find(any(coarse.faces == ii,2))';
        for vol = tvolumes
            lface = intersect(setdiff(coarse.elfaces(vol,:),0), tfaces);
            quad = [ii, refM(lface(1))  refC(vol), refM(lface(2))];
            elem = [elem; quad];
        end        
    end
    dual.coord = nodes;
    dual.elem = elem;
    dual.npar = size(elem,1);
    [F, tfaces,bflag] = createFaces(elem);  
    dual.faces = F;
    dual.elfaces = tfaces;
    dual.bflag = bflag;    
    dual.centelem = findcentelem(nodes, elem);
    %dual.elemloc = inVol(primal.coord,primal.elem);

end


function [tcentelem]  = findcentelem(tcoord, telem)
    t1 = tcoord(telem(:,1),:);
    t2 = tcoord(telem(:,2),:);
    t3 = tcoord(telem(:,3),:); 
    t4 = zeros(size(t3));
    ref = (telem(:,end) ~= 0);
    t4(ref,:) = tcoord(telem(ref,4),:);
    tcentelem = (t1+t2+t3+t4);
    tcentelem(~ref,:) = (1/3) * (tcentelem(~ref,:));
    tcentelem(ref,:) = (1/4) * (tcentelem(ref,:));    
end

function [elemloc] = inVol(tcoord, telem)
    global elem centelem
    npar = size(telem,1);
    elemloc = zeros(size(elem,1),1);
    for ii =1:npar
        tnodes = telem(ii,:);
        ref = (tnodes == 0);
        tnodes(ref) = [];
        center = mean(tcoord(tnodes,:));
        tnodes = [tnodes, tnodes(1)];
        xv = tcoord(tnodes,1);
        yv = tcoord(tnodes,2);
        cflag = inpolygon(center(1), center(2), xv, yv);
        cref = inpolygon(centelem(:,1), centelem(:,2), xv, yv);
        if cflag == 1
            elemloc(cref) = ii * 1;
        else
            elemloc(~cref) = ii * 1; 
        end
    end
end



function [face] = create_faceneigh(F)
    n = max(max(F));
    face = zeros(n,2);
    for ii = 1:n
       el = find(any(ismember(F,ii),2));
       m = max(size(el));
       if m == 1
           face(ii, 1) = el;
       elseif m == 2
           face(ii, :) = el;
       end       
    end
end


function [nref] = minDis(center, nref)
    global centelem
    points = centelem(nref,1:2);
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
    nref(nref==1) = out;
end
