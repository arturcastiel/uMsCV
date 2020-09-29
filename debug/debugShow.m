
hh = 9;


   flag = (~(primal_forming.bflag));
   
   el = primal_forming.nface(flag,:)
    %coarse_interface_center = coarse_interface_center(~primal_forming.bflag);
    cneigh = find(coarseneigh(hh,1:end-4));
    
ref = ismember(el,hh);
ref = any(ref,2);
intern = unique(horzcat(coarse_strips{find(ref)}));



figHandle = figure;
figure(figHandle);

lb = [1,2,3,4];
lb = [158,159,160,161];

bad = ismember(elemloc,lb);


% bad = elemloc < 5;
% bad = elemloc > (max(elemloc) - 4);


xm = max(coord(:,1)) - min(coord(:,1));
ym = max(coord(:,2)) - min(coord(:,2));
len = 1000;
set(figHandle, 'Position', [0 0 len len*(ym/xm)])
set(figHandle,'color','w')
axis off

n = 1;
Ap = [0,0];
Bp = [0, n];
Cp = [n, n];
Dp = [n, 0];

points = [Ap;Bp;Cp;Dp; Ap];

hold on
%plot(points(:,1), points(:,2), 'LineWidth',2, 'color', [0 0 0],'LineStyle','--')

color1 =[1,0, 0] 
%fill(x,y,'r')

% 
all_el = union(find(GlobalBoundary), coarseElemCenter);
ref = ismember(inedge(:,3:4),all_el);
ref = all(ref,2);
ref = ~ref;
%ref = ~all(ismember(inedge(:,3:4), all_el),2);
aux = ismember(elemloc(inedge(:,3:4)), lb);
raux = all(aux,2);
ref(raux) = false;

%for ii = find(ref)'
for ii = find(~raux)'
    nodes = inedge(ii,1:2)';
    el_face = coord(nodes,:);
    plot(el_face(:,1), el_face(:,2), 'LineWidth',0.5, 'color', [0.85 0.85 0.85],'LineStyle','-')
end
% 


%internos

color1 =[1,0, 0]
intern = H{hh}';
edges_vol = setdiff(intern, coarseElemCenter);
mop = setdiff(intern, coarseElemCenter) ;
%mop = setdiff(mop, find(bad));
for el = mop
    if elem(el,4) == 0
        nodes = [elem(el,1:3), elem(el,1)]';
    else
        nodes = [elem(el,1:4), elem(el,1)]';
    end
    elpoint = coord(nodes,:);
    fill(elpoint(:,1),elpoint(:,2),color1)
    hold on
end
% 



color1 =[255, 255, 0] /255
%edges_vol = setdiff(intern, coarseElemCenter);
mop = coarseElemCenter(hh) ;
%mop = setdiff(mop, find(bad));
for el = mop
    if elem(el,4) == 0
        nodes = [elem(el,1:3), elem(el,1)]';
    else
        nodes = [elem(el,1:4), elem(el,1)]';
    end
    elpoint = coord(nodes,:);
    fill(elpoint(:,1),elpoint(:,2),color1)
    hold on
end

color1 =[255, 255, 0] /255
%edges_vol = setdiff(intern, coarseElemCenter);
mop = coarseElemCenter' ;
%mop = setdiff(mop, find(bad));
for el = mop
    if elem(el,4) == 0
        nodes = [elem(el,1:3), elem(el,1)]';
    else
        nodes = [elem(el,1:4), elem(el,1)]';
    end
    elpoint = coord(nodes,:);
    fill(elpoint(:,1),elpoint(:,2),color1)
    hold on
end







% 
color1 =[160,255, 255]/255 
%edges_vol = setdiff(intern, coarseElemCenter);
mop = boundRegion{hh} ;
%mop = setdiff(mop, find(bad));
for el = mop
    if elem(el,4) == 0
        nodes = [elem(el,1:3), elem(el,1)]';
    else
        nodes = [elem(el,1:4), elem(el,1)]';
    end
    elpoint = coord(nodes,:);
    fill(elpoint(:,1),elpoint(:,2),color1)
    hold on
end

%
% color2 = [1, 1, 0] 
% 
% for el = setdiff(coarseElemCenter',find(bad))
%     if elem(el,4) == 0
%         nodes = [elem(el,1:3), elem(el,1)]';
%     else
%         nodes = [elem(el,1:4), elem(el,1)]';
%     end
%     elpoint = coord(nodes,:);
%     fill(elpoint(:,1),elpoint(:,2),color2)
%     hold on
% end

refb = all(ismember(elemloc(inedge(edgesOnCoarseBoundary,3:4)),lb),2);

for ii = edgesOnCoarseBoundary(~refb)'
    nodes = inedge(ii,1:2)';
    el_face = coord(nodes,:);
    plot(el_face(:,1), el_face(:,2), 'LineWidth',2, 'color', [0 0 0],'LineStyle','-')
end
% 
% 
for ii = 1:size(bedge,1)'
    
    nodes = bedge(ii,1:2)';
    el_face = coord(nodes,:);
    plot(el_face(:,1), el_face(:,2), 'LineWidth',2, 'color', [0 0 0],'LineStyle','-')
end