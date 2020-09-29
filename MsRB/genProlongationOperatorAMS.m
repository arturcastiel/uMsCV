function [ OP, CT] = genProlongationOperatorAMS(TransF, F)
%genProlongationOperator Generates Prolongation Operator
%   INPUT:
%   nelem = number of elements in the mesh
%   npar = number of partions
%   maxint = max number of iteration
    type_flag = 1;

    global npar    edges_ordering perm_matrix id_classify local_edges  I internal_split; 
    qq = perm_matrix*I;
    
    pTransF = perm_matrix*TransF*perm_matrix';
    pF = perm_matrix*F;
    type_vec = perm_matrix*id_classify;
  
    
    nodes = type_vec == max(type_vec);  
    %trange = max(edges_ordering);
%     tmax = max(id_classify) - 1;
%     tmin = (tmax - trange) + 1;
    
    %edges = (type_vec == (max(type_vec) -1));
    edges = type_vec == (max(type_vec) - 1);
%     global piq;
%     edges = (type_vec > piq) & (type_vec < max(type_vec));
    
    %edges = (type_vec <= tmax) & (type_vec >= tmin);
    
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
    %dog = blockDiag(perm_matrix * edges_ordering);
    
%     Nee =  (Mee.*dog) + (diag(sum(Mee.* ~dog,2)));
%     Mee = Nee;
%     Mtp = solveMeeAen(Mee,Aen);
%     Btp = -Mtp;
%     Btp = bsxfun(@rdivide, Btp(:,:) ,sum( Btp(:,:) ,2));

    %Mtp = (Mee\Aen);    
    %Btp = -Mtp;
    if type_flag == 0
        Btp = -(Mee\Aen);
    elseif type_flag == 1
        Mtp = solveMeeAen(Mee,Aen);
        Btp = -Mtp;
        Btp = bsxfun(@rdivide, Btp(:,:) ,sum( Btp(:,:) ,2)); 
        %Btp = normNegativo(Btp, local_edges);
        %Btp = normB(Btp, local_edges);
        %Btp = bsxfun(@rdivide, Btp(:,:) ,sum( Btp(:,:) ,2)); 
        %ref = true(size(Btp,1), 1);
        %Btp(ref,:) = bsxfun(@rdivide, Btp(ref,:) ,sum( Btp(ref,:) ,2));     
    elseif type_flag == 2
        [iMee] = finvMee(Mee);
        %[Aen] = fAee(Aen);
        Btp = -iMee*Aen;
        ref = true(size(Btp,1), 1);
        %Btp(ref,:) = bsxfun(@rdivide, Btp(ref,:) ,sum( Btp(ref,:) ,2));     
        Btp = bsxfun(@rdivide, Btp(:,:) ,sum( Btp(:,:) ,2));  
    end
 
    %Mee = mMee;
    
    
    T = (Ase*(-Btp) - Asn);
    %[T] = fixT(T);
    
    B = [Ass\T ; Btp; eye(npar,npar)];
    %B = [Ass\(Ase*(-Btp) - Asn) ; Btp; eye(npar,npar)];

    
   % B = [Ass\(Ase*Mtp - Asn) ; Btp; eye(npar,npar)];
   %B = bsxfun(@rdivide, B(:,:) ,sum( B(:,:) ,2));        

   OP = perm_matrix'*B;
   PQ = perm_matrix'*qq; 
    
    
    %CT = [Ass\Qs - Ass\(Ase* (Mee\Qe)); Mee\Qe; sparse(npar,1) ];
    if type_flag == 0
        CT = [Ass\Qs - Ass\(Ase* (Mee\Qe)); Mee\Qe; sparse(npar,1) ];
    else
        mMee = -(Btp / Aen);
        CT = [Ass\Qs - Ass\(Ase* (mMee*Qe)); mMee*Qe; sparse(npar,1) ];
    end
    CT = perm_matrix'*CT;
 
    %Assi = eye(size(Ass))\Ass;
    %C = [Assi -Assi*Ase*mMee  ; mMee; sparse(npar,1) ];

    
    %total = [Ass Ase Asn; zeros(size(Aes)) Mee Aen; zeros(size(Ans)) zeros(size(Ane)) eye(size(Ann))];
    
    %postprocessorTMS(full(sum(OP,2)),full(perm_matrix'*CT),0,superFolder,'Soma-da-galinha11');
end

% 
function [Tm] = fixT(T)
    global internal_split
    inf_outsupport = sum((~internal_split) .* T,2);
    inf_insupport = internal_split .*T;
    sum_support = sum(inf_insupport,2);
    
    ref = sum_support ~= 0;
    inf_insupport(ref,:) = bsxfun(@rdivide, inf_insupport(ref,:), sum_support(ref));     
    
    Tm = inf_insupport;
end

function [Mtp] = fAee(Aen)
  global npar local_edges perm_matrix elem inedge bedge edges coord
    %edge_size = size(Mee,1);
    %Mtp = zeros(size(Aen));    
    Mtp = Aen;
    for ii = 1:npar
        ref = (local_edges(:,ii) ~= 1) == true;
%         for jj = find(ref)'
%             auxA(:,jj) = 0;
%             auxA(jj,jj) = 1;
%             auxB(jj) = 0;
%         end
        Mtp(ref,ii) = 0;      
    end

end

function [Mtp] = finvMee(Mee)
    global edges_inv
    edge_size = size(Mee,1);
    Mtp = zeros( edge_size, edge_size);
    auxB = false(edge_size,1);
    for ii = 1:edge_size
        ref = (edges_inv(:,ii) == 1)   == true;
      
        auxA = Mee(ref,ref);
        auxC = auxB;
        auxC(ii) = true;
        auxC = auxC(ref);
%         for jj = find(ref)'
%             auxA(:,jj) = 0;
%             auxA(jj,jj) = 1;
%             %auxB(jj) = 0;
%         end
        Mtp(ref,ii) = auxA\auxC;   
    end

end
% 
% 
% function [Mtp] = finvMee(Mee)
%     global edges_inv
%     edge_size = size(Mee,1);
%     %Mtp = sparse( edge_size, edge_size);    
%     Mtp = zeros( edge_size, edge_size);    
% 
%     for ii = 1:edge_size
%         auxA = Mee;
%         auxB = false(edge_size,1);
%         auxB(ii)= true;
%         ref = (edges_inv(:,ii) ~= 1) == true;
%         auxA(:,ref) = 0;
%         auxA = auxA + diag(ref);
% 
%         %       auxC = auxB;
%         
% %         for jj = find(ref)'
% %             auxA(:,jj) = 0;
% %             auxA(jj,jj) = 1;
% %             %auxB(jj) = 0;
% %         end
%         Mtp(:,ii) = auxA\auxB;   
%     end
% end


function [Mtp] = solveMeeAen1(Mee,Aen)
    global npar local_edges perm_matrix elem inedge bedge edges coord
    edge_size = size(Mee,1);
    Mtp = zeros( edge_size, npar);    
    parfor ii = 1:npar
        auxA = Mee;
        auxB = Aen(:,ii);
        ref = (local_edges(:,ii) ~= 1) == true;
        for jj = find(ref)'
            auxA(:,jj) = 0;
            auxA(jj,jj) = 1;
            auxB(jj) = 0;
        end
        Mtp(:,ii) = auxA\auxB;      
    end

end


function [Mtp] = solveMeeAen(Mee,Aen)
    global npar local_edges perm_matrix elem inedge bedge edges coord
    edge_size = size(Mee,1);
    Mtp = zeros( edge_size, npar);    
    parfor ii = 1:npar
        auxA = Mee;
        auxB = Aen(:,ii);
        auxvec = zeros(size(auxB));
        ref = (local_edges(:,ii) == 0) == true;
        
        ref2 = (local_edges(:,ii) ~= 1) == true;        
        %linha de producao
        %ref = ref2;
        for jj = find(ref)'
            auxA(:,jj) = 0;
            auxA(jj,jj) = 1;
            auxB(jj) = 0;
        end
        
%         aa = auxA(ref,ref);
%         bb = auxB(ref);
%         
%         auxvec(ref) = aa\bb;
        auxvec = auxA\auxB;
        %auxvec = -normalizeNegative(-auxvec);
        auxvec(ref2) = 0;
        %auxvec(~ref2) = -normalizeNegative(-auxvec(~ref2));
        %auxvec(~ref2) = -normalizeNegative(-auxvec(~ref2))
        %Mtp(:,ii) = normalizeNegative(auxvec);
        Mtp(:,ii) = auxvec;      
    end

end

function [Mtp] = solveMeeAen2(Mee,Aen)
    global npar local_edges perm_matrix elem inedge bedge edges coord
    edge_size = size(Mee,1);
    Mtp = zeros( edge_size, npar);    
    parfor ii = 1:npar
        auxA = Mee;
        auxB = Aen(:,ii);
        ref = (local_edges(:,ii) == 0) == true;
        1
        for jj = find(ref)'
            auxA(:,jj) = 0;
            auxA(jj,jj) = 1;
            auxB(jj) = 0;
        end
        
        
        Mtp(:,ii) = auxA\auxB;      
    end

end

% function [Mtp] = solveMeeAen(Mee,Aen)
%     global npar local_edges perm_matrix elem inedge bedge edges coord
%     edge_size = size(Mee,1);
%     Mtp = sparse( edge_size, npar);    
%     for ii = 1:npar
%         auxA = Mee;
%         auxB = Aen(:,ii);
%         ref = (local_edges(:,ii) ~= 1) == true;
%         for jj = find(ref)'
%             auxA(:,jj) = 0;
%             auxA(jj,jj) = 1;
%             auxB(jj) = 0;
%         end
%         Mtp(:,ii) = auxA\auxB;      
%     end
% 
% end

function [C] = normNegativo(B, local_edges)
    ref = true(size(B,1), 1);
        %Btp(ref,:) = bsxfun(@rdivide, Btp(ref,:) ,sum( Btp(ref,:) ,2));     
    C = B;
    top = max(B, [], 2);
    bot = min(B, [], 2);
    dtp = top - bot;  
    %ref = (local_edges == 0) == true;
    C = bsxfun(@minus,C, bot);
   
    C((local_edges == 0) == true) = 0;

    C = bsxfun(@rdivide, C, dtp);
    dtp = sum(C,2);  
    %C((local_edges == 2) == true) = 0;
    
    C = bsxfun(@rdivide, C, dtp);
    %C((local_edges == 2) == true) = 0;

    B = C;
    %C  = F;
    %C = bsxfun(@minus, C(ref,:) ,bot); 
    %C = bsxfun(@rdivide, C(ref,:) , dtp); 
   1 
end


function [bnorm] = normB(B, local_edges)
global npar
bnorm = zeros(size(B));
for ii = 1:npar
    reff2 = abs(B(:,ii)) > 10^-8;
    reff1 = local_edges(:,ii) == 1;
    reff = reff1 & reff2;
    bnorm(reff,ii) = normalizeNegative(B(reff,ii));
end
1
end

function [dnorm] = normalizeNegative(auxvec)
     dmin = min(auxvec);
     dt = max(auxvec) - min(auxvec);
     dnorm = (auxvec - dmin)/dt;
       
end



function [blocDiag] = blockDiag(vec)
    lap = (vec ~= 0);
    only_edge = vec(lap);
    blocDiag = sparse(size(only_edge,1), size(only_edge,1));    
    for  ii = 1:max(only_edge)
       ref = (only_edge == ii);
       blocDiag(ref,ref) = 1;        
    end
end
