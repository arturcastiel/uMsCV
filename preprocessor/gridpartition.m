function [primal, dual, elemloc] = gridpartition(filename)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
    [primal] = createfprimal(filename);
    
    F = primal.faces;
    tcoord = primal.coord;
    % F = coarse.faces;
    % tcoord = coarse.coord;
    for ii = 1:size(F)
       drawLine(F(ii,1), F(ii,2), tcoord)
    end
    [dual] = createfdual(primal);
    [mdual] = createfdual(dual);
    elemloc = inVol(mdual.coord,mdual.elem);

end


function [coarse] = createfprimal(filename)
    [tcoord, tnode] = getcoord(filename);
    tcoord = tcoord(:,1:2);
    telem = getelemf(filename, tnode);
    % trian_flag = all(telem(:,4) == 0);
    [F, tfaces,bflag] = createFaces(telem);   
    npar = size(telem,1);
    [tcentelem]  = findcentelem(tcoord, telem);
%     t1 = tcoord(telem(:,1),:);
%     t2 = tcoord(telem(:,2),:);
%     t3 = tcoord(telem(:,1),:);
%     t4 = zeros(size(t3));
%     ref = (telem(:,end) ~= 0);
%     t4(ref,:) = tcoord(ref, :);       
%     tcentelem = (t1+t2+t3+t4).*(1 ./ (3*ones(size(t1,1),1) + ref));
    coarse = struct('coord', tcoord, 'elem',telem, 'faces', F, 'elfaces', tfaces, 'bflag', bflag, 'npar', npar, 'centelem', tcentelem);    
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
end


function [tcentelem]  = findcentelem(tcoord, telem)
    

    t1 = tcoord(telem(:,1),:);
    t2 = tcoord(telem(:,2),:);
    t3 = tcoord(telem(:,3),:);
   
    t4 = zeros(size(t3));
    ref = (telem(:,end) ~= 0);
    if any(ref)
        1
    end
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



