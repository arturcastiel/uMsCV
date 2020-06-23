function [ OP, CT] = genProlongationOperatorAMS(TransF, F)
%genProlongationOperator Generates Prolongation Operator
%   INPUT:
%   nelem = number of elements in the mesh
%   npar = number of partions
%   maxint = max number of iteration
    type_flag = 1;

    global npar    edges_ordering perm_matrix id_classify local_edges; 
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
        %Btp = bsxfun(@rdivide, Btp(:,:) ,sum( Btp(:,:) ,2));        
    elseif type_flag == 2
        [iMee] = finvMee(Mee);
        %[Aen] = fAee(Aen);
        Btp = -iMee*Aen;
        Btp = bsxfun(@rdivide, Btp(:,:) ,sum( Btp(:,:) ,2));  
    end
 
    %Mee = mMee;
    
    B = [Ass\(Ase*(-Btp) - Asn) ; Btp; eye(npar,npar)];

   % B = [Ass\(Ase*Mtp - Asn) ; Btp; eye(npar,npar)];
   %B = bsxfun(@rdivide, B(:,:) ,sum( B(:,:) ,2));        

   OP = perm_matrix'*B;
    
    
    
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
        ref = (edges_inv(:,ii) == 1) == true;
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


function [Mtp] = solveMeeAen(Mee,Aen)
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




function [blocDiag] = blockDiag(vec)
    lap = (vec ~= 0);
    only_edge = vec(lap);
    blocDiag = sparse(size(only_edge,1), size(only_edge,1));    
    for  ii = 1:max(only_edge)
       ref = (only_edge == ii);
       blocDiag(ref,ref) = 1;        
    end
end
