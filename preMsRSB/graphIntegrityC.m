function [nelemloc] = graphIntegrityC(elemloc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global inedge elem
npar = max(elemloc);



nelemloc = elemloc;
bad_elements = [];

for ii = 1:npar
    edges_boundary = nelemloc(inedge(:,3:4));
    ref = edges_boundary(:,1) ~= edges_boundary(:,2);
    lol = any(elemloc(inedge(ref,3:4)) == ii,2);
    ref(ref == true) = lol;
    
   
    all_bel = unique(inedge(ref,3:4));
    all_bel = all_bel(elemloc(all_bel) == ii);
    for jj = all_bel'
       refb = any(inedge(:,3:4) == jj,2);
       lmat = inedge(refb,3:4);
       if ~any(nelemloc(lmat(lmat ~= jj)) == ii)
           lep = nelemloc(lmat(lmat ~= jj));
           uv = unique(lep);
           n  = histc(lep,uv);
           refll = n == max(n);
           nelemloc(jj) = uv(refll);
       end
    end
   %belem = intersect((unique(elemloc(inedge(ref,3:4))),  find(elemloc == ii)));
end

end

