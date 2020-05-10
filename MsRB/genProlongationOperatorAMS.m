function [ OP, CT] = genProlongationOperatorAMS(TransF, F)
%genProlongationOperator Generates Prolongation Operator
%   INPUT:
%   nelem = number of elements in the mesh
%   npar = number of partions
%   maxint = max number of iteration
    global npar superFolder
    
    [perm_matrix, per_vec, id_classify] = genPermMatrix();
   
    pTransF = perm_matrix*TransF*perm_matrix';
    pF = perm_matrix*F;
    type_vec = perm_matrix*id_classify;
    nodes = type_vec == max(type_vec);
    edges = (type_vec == (max(type_vec) -1));
    internals = ~(nodes | edges); 
    
    Ass = pTransF(internals, internals);
    Ase = pTransF(internals, edges);
    Asn = pTransF(internals, nodes);
    
    Aes = pTransF(edges, internals);
    Aee = pTransF(edges, edges);
    Aen = pTransF(edges, nodes);
    
    Ans = pTransF(nodes, internals);
    Ane = pTransF(nodes, edges);
    Ann = pTransF(nodes, nodes); 
    
    Qs = pF(internals);
    Qe = pF(edges);
    Qn = pF(nodes);
    
    Mee = Aee + diag(sum(Aes,2));
%     Mnn = Ann + diag(sum(Ans,2)) + diag(sum(Ane,2));
%     Aes(:) = 0;
%     Ans(:) = 0;
%     Ane(:) = 0;
    %total = [Ass Ase Asn; Aes Mee Aen; Ans Ane Mnn];
    Mtp = (Mee\Aen);
    B = [Ass\(Ase*Mtp - Asn) ; -Mtp; eye(npar,npar)];
    OP = perm_matrix'*B;
    CT = [Ass\Qs - Ass\(Ase* (Mee\Qe)); Mee\Qe; sparse(npar,1) ];
    CT = perm_matrix'*CT;
    %postprocessorTMS(full(sum(OP,2)),full(perm_matrix'*CT),0,superFolder,'Soma-da-galinha11');
end