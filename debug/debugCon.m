
%% create connectivity matrix
cmat = sparse(size(elem,1),size(elem,1));
cmat(sub2ind(size(cmat),inedge(:,3),inedge(:,4))) = true;
cmat(sub2ind(size(cmat),inedge(:,4),inedge(:,3))) = true;
nCon = [];
for node =  1:size(coord,1)
    point = esurn2(node)+1;
    length = esurn2(node+1) - esurn2(node);
    element = esurn1(point :(point+length-1));
    nCon = [nCon; combnk(element,2)];
end
nCon = unique(nCon, 'rows');
cmat(sub2ind(size(cmat),nCon(:,1),nCon(:,2))) = true;
cmat(sub2ind(size(cmat),nCon(:,2),nCon(:,1))) = true;

r = symrcm(cmat)';
[~, rb] = sort(r);

% mm = zeros(size(r));
% for ii = 1:size(mm,1)
%    ref = (r == ii);
%    mm(ii) = find(ref);
% end

%% reoordering
centelem = centelem(r,:); 
elem = elem(r,:); 
elemarea = elemarea(r,:);
elemloc = elemloc(r,:);

esurn1 = rb(esurn1);
%transvec = 1:
bedge(:,3) = rb(bedge(:,3));

inedge(:,3:4) = rb(inedge(:,3:4));
for ii = 1:npar
   coarseelem{ii} =  rb(coarseelem{ii});
end

save('ord.mat','r','rb')
clear cmat
% for ii = 1:npar
%    coarseelem{ii} =  r(coarseelem{ii})';
% end
% % function [coord,centelem,elem,esurn1,esurn2,nsurn1,nsurn2,bedge,inedge,...
%     normals,esureface1,esureface2,esurefull1,esurefull2,elemarea,dens,visc,...
%     satlimit,pormap,bcflag,courant,totaltime,numcase,phasekey,kmap,wells,elemloc,npar,coarseelem ,ghostelem,...
%     multiCC, pmethod,mshfile, keymsfv,interptype,keypath1,smethod,xyrz,r0,symaxe,coarseratio,auxcvfactor,...
%     nonlinparam,multdopt,goefreeopt,order,timeorder,recovtype,lsneightype,...
%     lsexp,keygravity,g,keycapil,ncaplcorey,mainpath,resfolder,benchkey,...
%     klb,limiterflag,rowposit] = preprocessor 
%nCon = unique(nCon, 'rows');

%sub2ind(size(A),3,2)