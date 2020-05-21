function [ OP, CT] = genProlongationOperatorAMS(TransF, F)
%genProlongationOperator Generates Prolongation Operator
%   INPUT:
%   nelem = number of elements in the mesh
%   npar = number of partions
%   maxint = max number of iteration
    global npar superFolder edges_ordering
    
    [perm_matrix, per_vec, id_classify] = genPermMatrix();
   
    pTransF = perm_matrix*TransF*perm_matrix';
    pF = perm_matrix*F;
    type_vec = perm_matrix*id_classify;
    nodes = type_vec == max(type_vec);  
    trange = max(edges_ordering);
    tmax = max(id_classify) - 1;
    tmin = (tmax - trange) + 1;
    
    %edges = (type_vec == (max(type_vec) -1));
    
    edges = (type_vec <= tmax) & (type_vec >= tmin);
    
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
    dog = blockDiag(perm_matrix * edges_ordering);
    
%     Nee =  (Mee.*dog) + (diag(sum(Mee.* ~dog,2)));
%     Mee = Nee;
    Mtp = (Mee\Aen);
    B = [Ass\(Ase*Mtp - Asn) ; -Mtp; eye(npar,npar)];
    OP = perm_matrix'*B;
    CT = [Ass\Qs - Ass\(Ase* (Mee\Qe)); Mee\Qe; sparse(npar,1) ];
    CT = perm_matrix'*CT;
    
    %total = [Ass Ase Asn; zeros(size(Aes)) Mee Aen; zeros(size(Ans)) zeros(size(Ane)) eye(size(Ann))];
    
    %postprocessorTMS(full(sum(OP,2)),full(perm_matrix'*CT),0,superFolder,'Soma-da-galinha11');
    endc


function [blocDiag] = blockDiag(vec)
    lap = (vec ~= 0);
    only_edge = vec(lap);
    blocDiag = sparse(size(only_edge,1), size(only_edge,1));    
    for  ii = 1:max(only_edge)
       ref = (only_edge == ii);
       blocDiag(ref,ref) = 1;        
    end
end
