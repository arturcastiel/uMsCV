function [TransF, Q] = removeDirichlet(TransF,F)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global wells
dvec = sparse(size(F,1),size(F,2) );
doc = sparse(size(F,1),size(F,2) );

if size(wells,2) == 6
    ref = wells(:,5) > 400;
else
    ref = false(0,1);
end
wellelem = wells(ref,1);
dvec(wellelem) =   wells(ref,end);
doc(wellelem) =   1;
Q = F -  TransF*dvec;
TransF(wellelem, :) = 0;
TransF(:,wellelem) = 0;
Q(wellelem) = 0;
TransF  = TransF + diag(doc);
Q = Q + dvec;

end

