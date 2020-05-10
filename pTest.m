

F = mdual.faces;
tcoord = mdual.coord;
% F = coarse.faces;
% tcoord = coarse.coord;
for ii = 1:size(F)
   drawLine(F(ii,1), F(ii,2), tcoord)
end