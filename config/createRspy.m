%set(gca,'XColor', 'none','YColor','none')
function createRspy(aMat, edges,flag, id_type)
global perm_matrix
id_type = perm_matrix * id_type;
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

num = 0.015/2;

m =  num*n;



ll = max(id_type)

Ap = [0,0]  + [1-m, 1-m];
Bp = [0, n] + [1-m,  m];
Cp = [n, n] + [ m,  m];
Dp = [n, 0] + [m, 1-m];

points = [Ap;Bp;Cp;Dp; Ap];

for ii = 1:ll
   ref =  id_type == ii;
   color = [0.4,0.6,0.7]
   %color = [215/255, 229/255, 234/255];
   %color = [239/255, 206/255, 128/255];
   color = [255/255, 100/255, 100/255];

   
   paintRec(ref, color);
   
%    maxl =  max(find(ref));
%    minl =  min(find(ref));
%    len = maxl - minl;
%    xc = 0.5*(maxl + minl);
%    drawRectangle(len,len,xc,xc, [1, 1, 1/ii]);
end











hold on
plot(points(:,1), points(:,2), 'LineWidth',1, 'color', [0.75, 0.75, 0.75],'LineStyle','--')

hold on

   
tmp = find(edges);
ledges = min(tmp) - 0.5;
medges = max(tmp) + 0.5;



if middle_line 
    L1 = [Ap(:,1), ledges; Dp(:,1), ledges];
    L2 = [Ap(:,1), medges; Dp(:,1), medges];

    plot(L1(:,1), L1(:,2), 'LineWidth',1, 'color', [0.75, 0.75, 0.75],'LineStyle','--')
    plot(L2(:,1), L2(:,2), 'LineWidth',1, 'color', [0.75, 0.75, 0.75],'LineStyle','--')

    U1 = [ledges, Ap(:,2) ; ledges, Cp(:,2)];
    U2 = [medges, Ap(:,2) ; medges, Cp(:,2)];

    plot(U1(:,1), U1(:,2), 'LineWidth',1, 'color', [0.75, 0.75, 0.75],'LineStyle','--')
    plot(U2(:,1), U2(:,2), 'LineWidth',1, 'color', [0.75, 0.75, 0.75],'LineStyle','--')
end

%U2 = [Ap(:,1), medges; Dp(:,1), medges]


spy(Q,'k',1)

 xlim([min(points(:,1)) max(points(:,1))])
    ylim([min(points(:,2)) max(points(:,2))])
%hold off
end