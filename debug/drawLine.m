function [] = drawLine(p1, p2,coord)

x1 = coord(p1,1);
y1 = coord(p1,2);
x2 = coord(p2,1);
y2 = coord(p2,2);

plot([x1,x2],[y1,y2], 'LineWidth' , 1.3, 'color', [0,.8,0])


%plot([x1,x2],[y1,y2],'LineWidth' , 1.3)
hold on
end