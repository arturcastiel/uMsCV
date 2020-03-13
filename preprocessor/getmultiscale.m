%--------------------------------------------------------------------------
% GETMULTISCALE
%--------------------------------------------------------------------------
% Script responsable for generating all multiscale properties
%--------------------------------------------------------------------------
% OUTPUT:  cell : coarseedge{}
%                 intinterface{}
%                 exinterface{}
%                 exinterfaceaxes{}
%                 interfacecenter{}
%                 coarseblockcenter()
%                 numinterface()
%		Matrix:   coarseneigh()
%                 intCoord()
%                 edgesOnCoarseBoundary
%                 intDist 
%       Struct:   coarseningRatio.min
%                 coarseningRatio.max
%                 coarseningRatio.avg
%--------------------------------------------------------------------------
% INPUT: Scripts don't take input. This script uses all variables stated 
% inside preprocessor.m
%--------------------------------------------------------------------------
%% GENERATING COARSE CELL EDGES AND COARSE CELL NEIGHBORS
%--------------------------------------------------------------------------
% coarseedge{npar} - ALL EDGES OF A COARSE CELL
% Cell Array containing a matrix with the edges of each partion
% ex: coarseedge{3} - first collumn is a reference to either bedge or
% inedge and the second column is flag: 1 - look bedge | 0 - look inedge
% edge)
%   Reference Flag
%     3        0
%     4        0
%         ...
%     23       1
%--------------------------------------------------------------------------
% intinterface - coarse edge neighbors - EACH OTHER NEIGHBORS
% intinterface{m,n} gives the refrence of the edges that belong to h M and N
% intinterface(3,4) =
%   Reference 
%     3    
%     4   
%    ...
%     23 
%--------------------
%--------------------------------------------------------------------------
% exinterface- coarse edge on the boundary
% exinterface{m} gives the edges of coarse cell M that belong  to the
% boundary
%   Reference 
%     3      
%     4     
%    ...
%     23     
%--------------------------------------------------------------------------
% exinterfaceaxes- coarse edge on the boundary x
% exinterface{m} gives the edges of coarse cell M that belong  to the
% boundary
%      
%    ___2___       1 - Right
%   |       |      2  - Top
% 3 |       | 1    3 - Left
%   |_______|      4 - Bottom
%       4
%   Conection between external interface and axes 
%   Coarse Cell  Right Top Left Bottom
%                 1     2    3    4
%     1    
%     2
%    ...
%    number
%   of coarse
%     cells
% ex: exinterfaceaxes{1,3} - interface of partion 3 located on the left
%
%--------------------------------------------------------------------------
% interfacecenter{} - coordinates of each interface of all coarse cells
% ex: interfacecenter{2}  returns the center coordinate of the interface
% boundary from coarse block 2
%      
%    ___2___       1 -  [1 0.5]
%   |       |      2  - [0.5 1]
% 3 |       | 1    3 -  [0 0.5]
%   |_______|      4 -  [0.5 0]
%       4
%   The number of interfaces may very for each coarse cell
%   The return value is a n x 2 matrix, where n is the number of interfaces.
%
%--------------------------------------------------------------------------
% coarseblockcenter() - coordinates of the center of a coarse block
% ex: coarseblockcenter(3) returns the center of coarse block 3
%      
%    ___2___       
%   |       |      O - [0.25 0.25]
% 3 |   O   | 1   
%   |_______|      
%       4
%   coarseblockcenter is a matrix n x 2, where n is the number of 
%   coarseblocks
%   ex: coarseblockcenter(3) = [0.33 0.45]
%
%--------------------------------------------------------------------------
% numinterface() - return the number of interfaces (sides) of coarseblock
%    Interface 3:  
%    ___2___       
%   |       |      O - [0.25 0.25]
% 3 |   O   | 1   
%   |_______|      
%       4
%	 ex: numinterface(3)  returns 4 (Coarse block 3 has 4 interfaces)
%--------------------------------------------------------------------------
% coarseneigh(i,j) - matrix that returns 1 if i is neighbor to jj
% PS: coarseneigh is [npar npar+4]. 
%From npar+1 to npar+4 means 
%   Coarse Cell  Right      Top      Left    Bottom
%                npar+1    npar+2   npar+3   npar+4
%      
%    ___2___       1 - Right
%   |       |      2  - Top
% 3 |       | 1    3 - Left
%   |_______|      4 - Bottom
%       4

% A last collumn is added to tell if
% a coarse element is neighbor to the boundary.
%					
%   Ex: npar = 20  
%	 ex: coarseneigh(3,18) returns 1 if 3 and 18 are neighbors;0 otherwise
%		 coarseneigh(14,21) returns 1 if 14 is on the boundary top

%--------------------------------------------------------------------------
% intDist(i,j) - matrix the distance of interface of element I to element J
% PS: intDist(3,4) represents the distance of the interface of element 3 and 4. 
%--------------------------------------------------------------------------
% intCoord - general matrix (x 6)
% [ coordX  coord Y TypeOfNode Flag1 Flag2 Flag 3]
%
%  Type of Node:                    Flag 1          Flag 2       Flag 3
%   1 - Coarse Center Node       #Coarse Cell         0            0
%   2 - Interface betweeen      #Coarse Cell 1   #Coarse Cell 2    0
%          Coarse Cells
%   3 - Interface Coarse Cell   #Coarse Cell     #Boundary Type    0
%            Boundary                            (Egg: 1 - Right)
%   4 - Boundary Edges Point    #Coarse Cell    #BoundaryType #BoundaryType
%                                               (Egg: 1 - Right  x 2 - Top)
% starts only if working in a multiscale mesh

global flagboundcoarse

flagboundcoarse = zeros(size(coord,1),1);

if meshtype == true
    %checking bedges for edges on the border
    %edges from bedges come already oriented
    for i=1:npar
        tmp = find(elemloc( bedge(:,3)) ==i);
        s_tmp = size(tmp);
        exinterface{i} = tmp;
        coarseedge{i} = [tmp, ones(s_tmp(1),1)];
    end
 clear tmp s_tmp
 
    %checking inedges for edges on coarse boundary partion
    for ii = 1:npar
        for jj = 1:npar
            if ( ii ~= jj)
                intinterface{ii,jj} = [find(elemloc(inedge(:,3)) == ii & elemloc(inedge(:,4)) == jj) ; find(elemloc(inedge(:,3)) == jj & elemloc(inedge(:,4)) == ii)];
            end
        end
    end
    
    %adding internal edges in intiinterface to coarseedges
    
    for ii = 1:npar
        auxmat = [];
        for jj = 1:npar
            if ~isempty(intinterface{ii,jj})
                auxmat = [auxmat ;intinterface{ii,jj}];
            end
        end
        coarseedge{ii} = [coarseedge{ii} ;auxmat , zeros(size(auxmat))];
        
    end
    clear auxmat
    %% splitting edges on the boarder 
    for ii=1:npar
        t = size(exinterface{ii});
        for jj = 1:t(1)
            ref = exinterface{ii}(jj);
            p1 = coord(bedge(ref,1),1:2);
            p2 = coord(bedge(ref,2),1:2);
            dx = p1(1) - p2(1);
            dy = p1(2) - p2(2); 
            if (dy == 0) & (p1(2) == 0) & (p2(2) == 0) 
            % case edge lays on the bottom
                fp = 4;
            elseif (dy == 0) & (p1(2) ~= 0) & (p2(2) ~= 0) 
            % case edge lays on top
                fp = 2;
            elseif (dx == 0) & (p1(1) == 0) & (p2(1) == 0) 
            % case edge lays on the left
                fp = 3;
            else
            % case edge lays on the right
                fp = 1;                               
            end
      
            exinterfaceaxes{ii,fp} = [exinterfaceaxes{ii,fp} ; ref];
            exinterfaceaxes{ii,fp} = unique(exinterfaceaxes{ii,fp});
        end
        
    end
    
    
    %% collecting data on the number of interfaces and calculating the centers and its distances
%     numinterface = zeros([npar 1]);
%     interfacecenter = cell([npar,1]);
    %internal interfaces
    
    %allocating intCoord
    
    indexCoord = 1;
    for ii = 1 : npar
            
            acum = 0;
        for jj = ii+1 : npar
            if ~isempty(intinterface{ii,jj})
                intDist(ii,jj) = sum(disCalc(inedge(intinterface{ii,jj},1:2),coord  ));
                        
                point = centerinterface(intinterface{ii,jj},inedge, coord,elemloc);
                intCoord = [intCoord ; point , 2 , ii, jj,0];
                interfacecenter{ii} = [interfacecenter{ii};indexCoord];
                interfacecenter{jj} = [interfacecenter{jj};indexCoord];
                indexCoord = indexCoord + 1;
                acum = acum + 1;
            end
        numinterface(ii) = acum;
        
        end
    end
    

    for ii= 1:npar
        acum = 0;
        for jj = 1:4
            if ~isempty(exinterfaceaxes{ii,jj})
                point = centerinterface(exinterfaceaxes{ii,jj},bedge, coord,elemloc);
                intCoord = [intCoord ; point , 3 , ii, jj,0];

                interfacecenter{ii} = [interfacecenter{ii};indexCoord];
                indexCoord = indexCoord + 1;

                acum = acum + 1;
            end
        end
        numinterface(ii) = numinterface(ii) + acum;
    end
    

    
    
    %% Matrix that says which element is neighbor to which
%     
    for ii=1:npar
        for jj=1:4
            if ~isempty(exinterfaceaxes{ii,jj})
                 coarseneigh(ii,npar+jj) = coarseneigh(ii,npar+jj) +1;
            end
        end
    end
    
    coarseneigh = any(coarseneigh,3);
    %%checking if neighbors share at least a single point
    %substituir por coarsedge
    for ii = 1:npar
        comp1 = unique(inedge(coarseedge{ii}(find(coarseedge{ii}(:,2) == 0),1),1:2));       
            for jj = 1:npar
                comp2 = unique(inedge(coarseedge{jj}(find(coarseedge{jj}(:,2) == 0),1),1:2));
                check = intersect(comp1,comp2);                
                if ~isempty(check) & (ii ~= jj)
                   coarseneigh(ii,jj) = coarseneigh(ii,jj) +1;
                end
            end
        
    end
    
    %edges of the boundary
    
    for index = 1:npar
       for ii = 1:4
           for jj = (ii+1):4 
              if ~isempty(exinterfaceaxes{index,ii}) & ~isempty(exinterfaceaxes{index,jj});
                  comp1 = unique(bedge(exinterfaceaxes{index,ii},1:2));
                  comp2 = unique(bedge(exinterfaceaxes{index,jj},1:2));
                  tmp = intersect(comp1,comp2);
                  if ~isempty(tmp)
                      
                      intCoord = [intCoord ; coord(tmp,1:2) , 4 , index,ii, jj];
                      indexCoord = indexCoord + 1;                     
                  end                  
              end
           end
       end
    end
    
%% Calculating the center of internal elements
%     coarseblockcenter = zeros([npar 2]);
%     change here

disp('Calculating Multiscale Settings')
if multiCC == 1

    for ii = 1 : npar
        %point = weiszfeld(intCoord(interfacecenter{ii},1:2));
        
        if partionCenter == 1;
            %%TESTE FINAL
            %MT SCALE SPE
            point = centConstr(ii,2:3);
%             point(:,2) = point(:,2) * max(coord(:,1));
%             point(:,3) = point(:,3) * max(coord(:,2));
        else
            point = weiszfeld(intCoord(interfacecenter{ii},1:2));
        end
                       
        intCoord = [intCoord ; point , 1 , ii, 0,0] ;
        coarseblockcenter(ii,1) = indexCoord;
        indexCoord = indexCoord + 1;

    end
    elseif multiCC == 2
        for ii = 1 : npar
            if sum(coarseneigh(ii,end-3:end)) == 0               
                
                if partionCenter == 1;
                    point = centConstr(ii,2:3);
                else
                    %% testes
                    refTmp = [intCoord(interfacecenter{ii},4),intCoord(interfacecenter{ii},5)];
                    for lll = 1:size(refTmp,1)
                       weight(lll) =  intDist(refTmp(lll,1) ,refTmp(lll,2));
                    end
                    %% testes
                    point = weiszfeld(intCoord(interfacecenter{ii},1:2));
                    
                    %% points = weiszfeldMod(intCoord(interfacecenter{ii},1:2),  weight);
                end
            else
                ymax = max(coord(:,2));
                ymin = min(coord(:,2));
                xmax = max(coord(:,1));
                xmin = min(coord(:,1));
                if (coarseneigh(ii,end-3) == 1) & (coarseneigh(ii,end-2) == 1)
                    %RIGHT TOP
                    point = [xmax,ymax];
                elseif (coarseneigh(ii,end-2) == 1) & (coarseneigh(ii,end-1) == 1)
                    %TOP LEFT
                    %disp('topleft')
                    point = [xmin,ymax];
                elseif (coarseneigh(ii,end-1) == 1) & (coarseneigh(ii,end) == 1)
                    %LEFT BOTTOM
                    %disp('leftbottom')
                    point = [xmin,ymin];
                elseif (coarseneigh(ii,end) == 1) & (coarseneigh(ii,end-3) == 1)
                    %BOTTOM RIGHT
                    %disp('BOTTOM RIGHT')
                    point = [xmax,ymin];

                else
                    %disp('entrou aqui');

                    point = centerinterface(exinterface{ii},bedge, coord,elemloc);

                    
                end
               
            
            
            
            
            
            end
            
            intCoord = [intCoord ; point , 1 , ii, 0,0] ;
            coarseblockcenter(ii,1) = indexCoord;
            indexCoord = indexCoord + 1;
            
        end
        
    end
    
%% Calculating Coarsening Ratio

nElem = zeros(npar,1);


for index = 1:npar
    nElem(index) = size(coarseelem{index},2);    
end

coarseningRatio = struct('avg', size(elem,1)/npar, 'max',  max(nElem), ...
    'min', min(nElem),'std', std(nElem),'coarselem',npar, 'fineelem',size(elem,1),'cr', size(elem,1)/npar );
        
  %% Creating edgesOnCoarseBoundary
  
[a b] = size( intinterface);
edgesOnCoarseBoundary = [];
for ii = 1: a
   for jj = 1:b
       if  ~isempty(intinterface{ii,jj})
            edgesOnCoarseBoundary = union(edgesOnCoarseBoundary,intinterface{ii,jj});
       end
   end
end

%% Multiscale for MPFA-D
%oneNodetEdges is an edge that shares one node with an
%edgesOnCoarseBoundary and do not touch a dirichlet boundary

auxmat = unique([bedge(:,1) bedge(:,4); bedge(:,2) bedge(:,5)],'rows'); 
dirichtPoint = auxmat((auxmat(:,2) <200),:) ;

%dirichtPoint = []
tmp = ismember(inedge(:,1:2),setdiff(unique(inedge(edgesOnCoarseBoundary,1:2)),unique(dirichtPoint(:,1)) ));
Both = find(tmp(:,1) & tmp(:,2));
point1 = setdiff(find(tmp(:,1)), Both);
point2 = setdiff(find(tmp(:,2)), Both);


point1 = setdiff(point1, edgesOnCoarseBoundary);
point2 = setdiff(point2, edgesOnCoarseBoundary);
point3 = setdiff( Both, edgesOnCoarseBoundary);
% 
% 
% point1 = setdiff(point1, dirichtPoint);
% point2 = setdiff(point2, dirichtPoint);





oneNodeEdges = unique([point1 ones(size(point1,1),1);point2 2*ones(size(point2,1),1); point3 3*ones(size(point3,1),1)  ],'rows');




%pointloc returns all primal coarse cell a point belongs to
pointloc = cell(size(coord,1),1);

auxcell = cell(size(npar,1),1);

for ii = 1:npar
    
    auxcell{ii} = setdiff(unique(elem(coarseelem{ii},1:4)),0);
end

for ii = 1:size(coord,1)
    
   
    
    for jj = 1:npar
        
       if  sum(ismember(auxcell{jj},ii)) > 0
            pointloc{ii} = [ pointloc{ii}, jj]; 
       
       end
    end
    
    
end
clear auxcell

%%creating a structure that inputs the region and outputs the points to
%%recalculate the weights


tmp = inedge(oneNodeEdges(:,1),1:2);
%allPoints = unique(tmp(oneNodeEdges(:,2),1:2));

ref1 =  (oneNodeEdges(:,2)) == 1;
ref2 =  (oneNodeEdges(:,2)) == 2;
%%
ref3 =  (oneNodeEdges(:,2)) == 3;




%pointWeight = input coarse cell - n
%              output points in the coarse cell n that the weights must be recalculated 

%alteracao 07 04
allPoints  = unique( [tmp(ref1,1) ;tmp(ref2,2) ;tmp(ref3,1);tmp(ref3,2)   ] );
%allPoints = unique(tmp);

pointWeight = cell(npar,1);
for ii = 1: npar
  
   tmp =  inedge(edgesOnCoarseBoundary,3:4);
   ref = elemloc(tmp(:,1))== ii | elemloc(tmp(:,2)) == ii;

   localpoints = unique(inedge(edgesOnCoarseBoundary(ref), 1:2));
   
   
   acum = [];

   ref = elemloc(inedge(oneNodeEdges(:,1),3)) == ii | elemloc(inedge(oneNodeEdges(:,1),4)) == ii;
   acum = unique(inedge(oneNodeEdges(ref,1),1:2));
   
   
   if ~isempty(localpoints)
        
        pointWeight{ii} = localpoints;
        pointWeight{ii} = intersect(allPoints,localpoints);
        pointWeight{ii} = union(acum,pointWeight{ii});
        
   end
   
      
   flagboundcoarse(localpoints) = ones(size(localpoints));
   
   
end

% dynamic W wieght and dynamic S weight

for ii = 1:npar
    pointWeight{ii} = unique(elem(elemloc == ii,1:4));
end



regularEdges = unique( setdiff( 1:(size(inedge,1)), [oneNodeEdges(:,1)  ; edgesOnCoarseBoundary ]))';

bedgePoint = intersect( unique(inedge(edgesOnCoarseBoundary,1:2)), unique(bedge(:,1:2)));

auxmat = ismember(bedge(:,1:2), bedgePoint); 
bedgeNode = find(auxmat(:,1) | auxmat(:,2));

pq = setdiff(unique(bedge(bedgeNode,1:2)'),bedgePoint');
%

%% alteracoes para testar semiEdges

ref1 = oneNodeEdges(:,2) == 1;
ref2 = oneNodeEdges(:,2) == 2;
ref3 = oneNodeEdges(:,2) == 3;
point1 = unique(inedge(oneNodeEdges(ref1,1),2));
point2 = unique(inedge(oneNodeEdges(ref2,1),1));
point3 = unique(inedge(oneNodeEdges(ref3,1),1));
point4 = unique(inedge(oneNodeEdges(ref3,1),2));
allPoints = unique([point1; point2; point3; point4;pq ]);
pointBoundary = unique(bedge(:,1:2));

%check taking it off
%allPoints = setdiff(allPoints, pointBoundary );
auxmat = ismember(inedge(regularEdges,1:2), allPoints);



ref1 = auxmat(:,1) == 1;
ref2 = auxmat(:,2) == 1;
ref3 = sum(auxmat,2) == 2;

ref1 = logical(ref1 .*~ref3);
ref2 = logical(ref2 .*~ref3);

point1 = regularEdges(ref1); % auxmat2(ref1);
point2 = regularEdges(ref2); % auxmat2(ref2);
point3 = regularEdges(ref3); %auxmat2(ref3);
%semiEdges = unique([point1 ones(size(point1,1),1)

semiEdges = unique([point1 ones(size(point1,1),1);point2 2*ones(size(point2,1),1); point3 3*ones(size(point3,1),1)  ],'rows');

regularEdges = setdiff(regularEdges, semiEdges(:,1));



%% pointWeight bedge points


for ii = 1:npar
   ref = elemloc(bedge(bedgeNode,3)) == ii; 
   points = bedge(bedgeNode(ref),1:2);
   pointWeight{ii} =  setdiff(union(pointWeight{ii},points),0);
end

    %%
    %semiEdges = 0
    
    disp('Multiscale properties successfully generated');
else
    semiEdges = [];
    bedgeNode = [];
    
    coarseedge = {};
    intinterface = {};
    exinterface = {};
    exinterfaceaxes = {};
    interfacecenter = {};
    coarseblockcenter = [];
    numinterface = [];
    coarseneigh = [];
    intCoord = [];
    edgesOnCoarseBoundary = [];
    intDist = [];
    coarseningRatio = [];
    oneNodeEdges =[];
    pointloc = [];
    pointWeight = [];
    regularEdges = [];
    
%       Struct:   coarseningRatio.min
%                 coarseningRatio.max
%                 coarseningRatio.avg
    
    
end
