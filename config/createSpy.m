%set(gca,'XColor', 'none','YColor','none')
function createSpy(aMat, edges,flag)
colormat = load('color.dat');
figHandle = figure;
figure(figHandle);
set(figHandle, 'Position', [0 0 700 700])
set(figHandle,'color','w')
axis off
%hold on
%whitebg('white')
middle_line  = flag;

Q  = aMat;
n = size(Q,1);

num = 0.015;

m =  num*n;


Ap = [0,0]  + [1-m, 1-m];
Bp = [0, n] + [1-m,  m];
Cp = [n, n] + [ m,  m];
Dp = [n, 0] + [m, 1-m];

points = [Ap;Bp;Cp;Dp; Ap];



hold on
plot(points(:,1), points(:,2), 'LineWidth',1, 'color', [0.75, 0.75, 0.75],'LineStyle','--')

hold on

   
tmp = find(edges);
ledges = min(tmp) - 0.5;
medges = max(tmp) + 0.5;



if middle_line 
    L1 = [Ap(:,1), ledges; Dp(:,1), ledges]
    L2 = [Ap(:,1), medges; Dp(:,1), medges]

    plot(L1(:,1), L1(:,2), 'LineWidth',1, 'color', [0.75, 0.75, 0.75],'LineStyle','--')
    plot(L2(:,1), L2(:,2), 'LineWidth',1, 'color', [0.75, 0.75, 0.75],'LineStyle','--')

    U1 = [ledges, Ap(:,2) ; ledges, Cp(:,2)];
    U2 = [medges, Ap(:,2) ; medges, Cp(:,2)];

    plot(U1(:,1), U1(:,2), 'LineWidth',1, 'color', [0.75, 0.75, 0.75],'LineStyle','--')
    plot(U2(:,1), U2(:,2), 'LineWidth',1, 'color', [0.75, 0.75, 0.75],'LineStyle','--')
end

%U2 = [Ap(:,1), medges; Dp(:,1), medges]


spy(Q,1)

 xlim([min(points(:,1)) max(points(:,1))])
    ylim([min(points(:,2)) max(points(:,2))])
%hold off
end