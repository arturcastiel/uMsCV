%%creates a matrix with all edges coarseedge
global inedge bedge coord coarseedge
colormat = load('color.dat');

disp('Plotando Coarse Mesh')
figHandle = figure;
figure(figHandle);
%spy(TransF, 'r')
spy(perm_matrix*TransF*perm_matrix', 'b')
xSize = max(coord(:,1));
ySize = max(coord(:,2));

mxSize = min(coord(:,1));
mySize = min(coord(:,2));
xp = xSize - mxSize;
yp = ySize - mySize;

ratio = xp/yp;
set(figHandle, 'Position', [0 0 700 700/ratio])
set(gca, 'XTick',[])
set(gca, 'xticklabel',[])
set(gca, 'YTick',[])
set(gca, 'yticklabel',[])
text = num2str(npar);
%set(figHandle, 'name', ['Coarse Mesh: ', text ,' coarse cells'],'NumberTitle','off');
%title(['\fontsize{16}black {\color{magenta}magenta '...
%'\color[rgb]{0 .5 .5}teal \color{red}red} black again'])



    %title(['Coarse Mesh: ', text ,' coarse cells'],'FontSize', 14)
    whitebg('white')
    xlim([mxSize xSize])
    ylim([mySize ySize])
    
    hold on