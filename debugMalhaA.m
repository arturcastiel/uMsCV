figHandle = figure;
figure(figHandle);


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


all_el = union(find(GlobalBoundary), coarseElemCenter);

ref = ismember(inedge(:,3:4),all_el);
ref = all(ref,2);
ref = ~ref;

%ref = ~all(ismember(inedge(:,3:4), all_el),2);

for ii = find(ref)'
    nodes = inedge(ii,1:2)';
    el_face = coord(nodes,:);
    plot(el_face(:,1), el_face(:,2), 'LineWidth',1, 'color', [0.85 0.85 0.85],'LineStyle','-')
end

%edges_vol = find(GlobalBoundary');
for el = find(GlobalBoundary')
    if elem(el,4) == 0
        nodes = [elem(el,1:3), elem(el,1)]';
    else
        nodes = [elem(el,1:4), elem(el,1)]';
    end
    elpoint = coord(nodes,:);
    fill(elpoint(:,1),elpoint(:,2),color1)
    hold on
end

color2 = [1, 1, 0] 

for el = coarseElemCenter'
    if elem(el,4) == 0
        nodes = [elem(el,1:3), elem(el,1)]';
    else
        nodes = [elem(el,1:4), elem(el,1)]';
    end
    elpoint = coord(nodes,:);
    fill(elpoint(:,1),elpoint(:,2),color2)
    hold on
end



for ii = edgesOnCoarseBoundary'
    nodes = inedge(ii,1:2)';
    el_face = coord(nodes,:);
    plot(el_face(:,1), el_face(:,2), 'LineWidth',2, 'color', [0 0 0],'LineStyle','-')
end


for ii = 1:size(bedge,1)'
    
    nodes = bedge(ii,1:2)';
    el_face = coord(nodes,:);
    plot(el_face(:,1), el_face(:,2), 'LineWidth',2, 'color', [0 0 0],'LineStyle','-')
end

% 
% points = [Ap;Bp;Cp;Dp; Ap];
% 
% hold on
% plot(points(:,1), points(:,2), 'LineWidth',2, 'color', [0 0 0],'LineStyle','--')
%