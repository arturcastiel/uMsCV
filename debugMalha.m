figHandle = figure;
figure(figHandle);
set(figHandle, 'Position', [0 0 700 700])
set(figHandle,'color','w')
axis off

n = 1;
Ap = [0,0];
Bp = [0, n];
Cp = [n, n];
Dp = [n, 0];

points = [Ap;Bp;Cp;Dp; Ap];

hold on
plot(points(:,1), points(:,2), 'LineWidth',2, 'color', [0 0 0],'LineStyle','--')

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
    nodes = [elem(el,1:3), elem(el,1)]';
    elpoint = coord(nodes,:);
    fill(elpoint(:,1),elpoint(:,2),color1)
    hold on
end

color2 = [1, 1, 0] 

for el = coarseElemCenter'
    nodes = [elem(el,1:3), elem(el,1)]';
    elpoint = coord(nodes,:);
    fill(elpoint(:,1),elpoint(:,2),color2)
    hold on
end




points = [Ap;Bp;Cp;Dp; Ap];

hold on
plot(points(:,1), points(:,2), 'LineWidth',2, 'color', [0 0 0],'LineStyle','--')

