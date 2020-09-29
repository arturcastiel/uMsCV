function [nelemloc] = graphIntegrity2(elemloc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global inedge elem
npar = max(elemloc);



nelemloc = elemloc;
bad_elements = [];
for ii = 1:npar
    lref = (nelemloc == ii);
    transvec = find(lref);
    transback = zeros(size(elem,1),1);
    transback(transvec) = 1:size(transvec,1);
    edge_ref = all(ismember( inedge(:,3:4), transvec),2);
    auxmat = transback(inedge(edge_ref,3:4));
    G = graph(auxmat(:,1), auxmat(:,2));
    color_coarse =  conncomp(G);
    
    if max(color_coarse) ~= 1
        1+1;
    end
    if max(color_coarse) ~= 1
        [a, b] = hist(color_coarse, unique(color_coarse));
        unconnected_pieces = b(~(a == max(a)));
        ref = ismember(color_coarse, unconnected_pieces);
        bad_elements = [bad_elements; transvec(find(ref))];
    end
    nelemloc(bad_elements) = NaN;
end
for ii = bad_elements'
    auxmat = inedge(any(inedge(:,3:4) == ii,2),3:4);
    auxvec =  reshape(auxmat.',1,[]);
    auxvec =  nelemloc(auxvec(auxvec ~= ii));
    [a, b] = hist(auxvec, unique(auxvec));
    ref = a == max(a);
    if all(ref)
        ref(:) = false;
        ref(1) = true;
    end
    color = b(ref);
    nelemloc(ii) = color;
end
end

