%ref dos elementos com termo fonte
global sourceTerm actTerm bnodes refDir


ref1 =    (centelem(:, 1) <  5/8) & (centelem(:, 1) >  3/8);
ref2 =    (centelem(:, 2) <  5/8) & (centelem(:, 2) >  3/8);
ref = ref1 & ref2;

p1 = coord(elem(:,1),:);
p2 = coord(elem(:,2),:);
p3 = coord(elem(:,3),:);
%ref(:) = 1;
A = p1(ref,1).* p2(ref,2) +  p2(ref,1).* p3(ref,2) +  p3(ref,1).* p1(ref,2);
B = p1(ref,2).* p2(ref,1) +  p2(ref,2).* p3(ref,1) +  p3(ref,2).* p1(ref,1);
S = (A-B)./2;
actTerm = ref;
sourceTerm = S;



ee = 5*(10^-2);
kxx = (centelem(:,2).^2) + (ee * (centelem(:,1).^2));
kxy = -(1- ee)*(centelem(:,1) .* centelem(:,2));
kyy= (ee * (centelem(:,2).^2)) +  (centelem(:,1).^2);

kmap = [ [1:size(elem,1)]' , kxx, kxy, kxy, kyy];
kmap(:,2) = 1;
kmap(:,3) = 0;
kmap(:,4) = 0;
kmap(:,5) = 1 ;



wells(:,4) = 0;
wells(:,end) = 0;
wells = [];
bnodes = unique(bedge(:,1:2));

refDir = (any(ismember(elem(:,1:4),[1,2,3,4]),2));