global centelem
filename = 'cemfra2.msh';
file = [pwd, '/grids/coarse/', filename];
[tcoord,tnode] = getcoord(file);
tcoord = tcoord(:,1:2);
telem = getelem(file,tnode,0,0);
trian_flag = all(telem(:,4) == 0);
if trian_flag == 1
    telem = telem(:,1:3);
else
    telem = telem(:,1:4);
end
npar = size(telem,1);
elemloc = zeros(size(elem,1),1);

for ii =1:npar
    tnodes = telem(ii,:);
    center = mean(tcoord(tnodes,:));
    tnodes = [tnodes, tnodes(1)];
    xv = tcoord(tnodes,1);
    yv = tcoord(tnodes,2);
    cflag = inpolygon(center(1), center(2), xv, yv);
    cref = inpolygon(centelem(:,1), centelem(:,2), xv, yv);
    if cflag == 1
        elemloc(cref) = ii * 1;
    else
        elemloc(~cref) = ii * 1; 
    end
end