debugTest3;
global npar edges local_edges perm_matrix
ii = 14;

ref = local_edges(:,ii) == 1;

bef = (perm_matrix*edges == 1);
bef (bef == 1) = local_edges(:,ii);

bb = perm_matrix' * bef;
bb = full((bb ==1));

meshplot(elem(bb,1:4), coord(:,1), coord(:,2),'color',[0.4, 0.55 ,0.75], 'LineWidth' , 1.7);

debugInt