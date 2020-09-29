global internal_split perm_matrix elemloc edges_inv local_edges edges


teste = zeros(size(OP_old));
map = zeros(size(OP_old));
lap = zeros(size(OP_old));
edgesp = perm_matrix*edges > 0;
for ii = 1:npar
    auxloc = zeros(size(elemloc));
    auxloc2 = zeros(size(elemloc));
        auxloc3 = zeros(size(elemloc));

    auxloc(internal_split(:,ii)) = 1;
    auxloc2(edgesp) = (local_edges(:,ii) == 1);
    auxloc3(edgesp) = (local_edges(:,ii) == 2);

    teste(:,ii) = perm_matrix'*auxloc;
    map(:,ii) = perm_matrix'*auxloc2;
    lap(:,ii) = perm_matrix'*auxloc3;

    
end
postprocessorOP(OP_old,1,  superFolder, 'TESTE_01')
postprocessorOP(teste,1,  superFolder, 'TESTE_02')

postprocessorOP(map,1,  superFolder, 'TESTE_03')
postprocessorOP(lap,1,  superFolder, 'TESTE_04')


