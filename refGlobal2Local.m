%this script finds the reference of the centers in the coarseelem

refCenterInCoaseElem = zeros(npar,1);
for ii = 1:npar
    set_elem = find(elemloc == ii);
    refCenterInCoaseElem(ii) = find(set_elem == coarseElemCenter(ii));    
end


 for ii = 1:npar
     set_elem = find(elemloc == ii);
     dictionary(ii) = dict();
     dictionary(ii).map = containers.Map(set_elem,[1:size(set_elem,1)] );
 end
 

 edgesCoarseDict = containers.Map(edgesOnCoarseBoundary',[1:size(edgesOnCoarseBoundary,1)] );
 
 
%  %% Checking if are there any coarse elements centroids on a dirichlet boundary
% 
% frontElem = find(sum(coarseneigh(:,npar+1:end),2));
% 
% centFront = coarseElemCenter(frontElem);
% 
% 
% refBed = (bedge(:,4) < 200 | bedge(:,5) < 200 );
% 
% centDircht = intersect( bedge(refBed,3),centFront);
% 
% coarseDiricht = find(ismember(coarseElemCenter,centDircht));
coarseDiricht =0