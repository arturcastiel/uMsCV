debugTest3
ii = 10

ref = (cop(:,1) == ii) | (cop(:,2) == ii);
m = coarse_strips(find(ref));
m =  unique(horzcat(m{:}));
m = m(elemloc(m) == ii);
meshplot(elem(m,1:4), coord(:,1), coord(:,2),'color',colormat(6,:), 'LineWidth' , 1.7);

meshplot(elem(coarseElemCenter,1:4), coord(:,1), coord(:,2),'color',colormat(90,:), 'LineWidth' , 2.1);
leb = unique(cqb);
meshplot(elem(leb,1:4), coord(:,1), coord(:,2),'color',colormat(3,:), 'LineWidth' , 2.1);

%meshplot(elem(boundRegion{13},1:4), coord(:,1), coord(:,2),'color',colormat(6,:), 'LineWidth' , 1.7);
