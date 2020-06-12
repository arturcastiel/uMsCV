function [perm_matrix, nodes, edges, id_classify] = genPermMatrix()
%genProlongationOperator Generates Prolongation Operator
    [id_classify, nodes, edges, internals] = createClassVec();
    [perm_matrix, per_vec] = createPerMat(id_classify);
end


function [id_classify, nodes, edges, internals] = createClassVec()
    global coarseElemCenter elemloc GlobalBoundary  inedge
    id_classify = zeros(size(elemloc));
    edges = false(size(elemloc));
    nodes = false(size(elemloc));
    edges(GlobalBoundary) = true;
    nodes(coarseElemCenter) = true;
    internals = ~(edges | nodes);
    ref = all(ismember(inedge(:,3:4), find(internals)),2);
    transvec = find(internals);
    transback = zeros(size(elemloc,1),1);
    transback(transvec) = 1:size(transvec,1);
    auxmat = transback(inedge(ref,3:4));
    G = graph(auxmat(:,1) , auxmat(:,2));
    id_classify(transvec) = conncomp(G);
    max_id = max(id_classify);
    id_classify(edges) = max_id + 1;% + edges_ordering(edges);
    %max_id = max(id_classify);
    id_classify(nodes) = max_id + 2;
end

function [matrix, r] = createPerMat(id_classify)
    [~, r] = sort(id_classify);
    ord = @(i,j,A)( (j- ones(size(j)))*size(A,1) +i);
    matrix = sparse(size(id_classify,1),size(id_classify,1));
    ind = ord([1:size(id_classify,1)]', r, matrix);
    matrix(ind) = true;
end



function [id_classify, nodes, edges, internals] = createClassVec2()
    global coarseElemCenter elemloc GlobalBoundary  inedge piq
    id_classify = zeros(size(elemloc));
    nodes = false(size(elemloc));
    edges = GlobalBoundary;
    nodes(coarseElemCenter) = true;
    internals = ~(edges | nodes);
    
    %%classifying internals
    ref = all(ismember(inedge(:,3:4), find(internals)),2);
    transvec = find(internals);
    transback = zeros(size(elemloc,1),1);
    transback(transvec) = 1:size(transvec,1);
    auxmat = transback(inedge(ref,3:4));
    G = graph(auxmat(:,1) , auxmat(:,2));
    id_classify(transvec) = conncomp(G);
    
    
    max_id = max(id_classify);
   
    ref = all(ismember(inedge(:,3:4), find(edges)),2);
    transvec = find(edges);
    transback = zeros(size(elemloc,1),1);
    transback(transvec) = 1:size(transvec,1);
    auxmat = transback(inedge(ref,3:4));
    G = graph(auxmat(:,1) , auxmat(:,2));
    id_classify(transvec) = conncomp(G) + max_id;    
    piq = max_id;
    
     max_id = max(id_classify);
    
    
    %id_classify(edges) = max_id + 1;% + edges_ordering(edges);
    %max_id = max(id_classify);
    id_classify(nodes) = max_id + 1;
end