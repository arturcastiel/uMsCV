

F = primal_forming.faces;
%F = forming_primal.faces;

%F = aa.faces;

tcoord = primal_forming.coord;
%tcoord = forming_primal.coord;
%tcoord = aa.coord;

%bflag = primal_forming.bflag;
% F = coarse.faces;
% tcoord = coarse.coord;
for ii = 1:size(F)
    %if bflag(ii)
        drawLine(F(ii,1), F(ii,2), tcoord)
    %end
end

pcoord = primal_forming.coord(primal_forming.bnodes,:);
pcoord = primal_forming.centelem;

scatter(pcoord(:,1), pcoord(:,2));