function [elem,nbe,nelem,nodelim,intnode,flaglim,inboundedge] = ...
    getelem(verifymshfile,nnode,numwell,well)
%Initialize "inboundedge"
inboundedge = 0;
%Identify on the *.msh file a like-element type. Actualy, the *.msh
%gives entities which hold points, edges and elements (like-elements).
%This entities will all of them be reads.
%Open the *.msh file.
readmsh = fopen(verifymshfile);
        
%"nent" is the number of all entities understud as element in *.msh file
getmshdata = textscan(readmsh,'%u',1,'HeaderLines',7 + nnode);
%Attribute the data to "nent" (number of entity, that is: points, lines 
%and elements)
nent = getmshdata{1};
        
%Once we have "nent", both like-element and external edges entities may be 
%located. "meshtype" is the type of element used in domain discretization.
%2 ==> triangles; 3 ==> quadrangles. But that will be used also to define
%how many points entities there on the *.msh file

%There is an INTERNAL LINE
%debug

if numwell > 0 && any(well(:,9) == 3) || numwell > 0 && any(well(:,9) == 0)
    %Swept the "elements" of "*.msh" in order to fill the parameters above. 
    getmshdata = textscan(readmsh,'%*n%n%*n%n%*n%n%n%*[^\n]',nent);
    %Attribute the data to "auxmat"
    auxmat = cell2mat(getmshdata);
    %Gives the contribution to "entitype"
    entitype = auxmat(:,1);

    %Get "inboundedge" from "auxmat"
    inboundedge = ...
        auxmat(logical(auxmat(:,1) ~= 15 & auxmat(:,2) > 1000),3:4);
    %Orient in encreasing order the vertices which define each boundary
    %edge (ex.: 5, 1 ==> 1, 5).
    inboundedge = sort(inboundedge,2);

    %"nbe" verifies how many entities are edges which constitute the 
    %boundary (scalar)
    nbeaux = sum(entitype == 1);
    nbe = nbeaux - size(inboundedge,1);
    %"nodelim". These nodes constitute the limits of geometry. 
    %Get the amount of internal node (scalar)
    intnode = sum(logical(auxmat(:,1) == 15 & auxmat(:,2) > 1000));
    %We exclude any node inside the domain.    
    nodelim = sum(entitype == 15) - intnode;
%There is no well
else
    %Swept the "elements" of "*.msh" in order to fill the parameters above. 
    getmshdata = textscan(readmsh,'%*n %n %*n %n %*[^\n]',nent);
    %Attribute the data to "auxmat"
    auxmat = cell2mat(getmshdata);
    %Fill "entitype"
    entitype = auxmat(:,1);
    %"nbe" verifies how many entities are edges which constitute the 
    %boundary
    nbe = sum(entitype == 1);
    nbeaux = nbe;
    %"nodelimit". These nodes constitute the limits of geometry
    nodelim = sum(entitype == 15);
    intnode = 0;
end  %End of IF

%Define the amount of each entity:
%Verify how many "entities" are nodes. This value is attributed to 
%"ntri" verifies how many triangles there are in the domain
ntri = sum(entitype == 2); 
%"nquad" verifies how many quadrangles there are in the domain
nquad = sum(entitype == 3);
%"nelem" is the number of elements in the domain (triangles and/or 
%quadrangles)
nelem = ntri + nquad;

%Get "flaglim" from "auxmat". It is a vector with the flag associated to
%external nodes.
flaglim = auxmat(logical(auxmat(:,1) == 15 & auxmat(:,2) < 1000),2);

%Close the *.msh file
fclose(readmsh);

%---------------------------
%Define the matrix ("elem"):

%Open again the *.msh file
readmsh = fopen(verifymshfile);

%"elem" is a matrix which contain the nodes which constitute it
%Create and initialize the parameter "elem"
elem = zeros(nelem,5);  %nelem rows and 5 columns

%Fill the matrix "elem" with 4 nodes for a quadrangle element and 3 nodes 
%for a triangle element. In this case the fourth column receives null
%value. The fiveth column receive material propertie flag.
%There are tringles:
if ntri > 0  
    %Triangle
    getmshdata = textscan(readmsh,'%*u %*u %*u %u %*u %u %u %u',ntri,...
        'HeaderLines',8 + nnode + nodelim + intnode + nbeaux);
    elemaux = cell2mat(getmshdata);
    %Attribute to "elem" (fifth column)
    elem(1:ntri,5) = elemaux(:,1); 
    %Attribute to "elem" (another column)
    elem(1:ntri,1:3) = elemaux(:,2:4); 

    %Close the *.msh file
    fclose(readmsh);
end  %End of first IF

%There are quadrangles:
if nquad > 0  
    %Open again the *.msh file
    readmsh = fopen(verifymshfile);
    
    %Quadrangles
    getmshdata = textscan(readmsh,'%*u %*u %*u %u %*u %u %u %u %u',nquad,...
        'HeaderLines',8 + nnode + nodelim + intnode + nbeaux + ntri);
    elemaux = cell2mat(getmshdata);
    %Attribute to "elem" (fifth column)
    elem(ntri + 1:ntri + nquad,5) = elemaux(:,1); 
    %Attribute to "elem" (another column)
    elem(ntri + 1:ntri + nquad,1:4) = elemaux(:,2:5); 

    %Close the *.msh file
    fclose(readmsh);
end  %End of second IF

%Clear matrices
clear auxmat getmshdata elemaux; 

%It avoids that the fifth column of "elem" get null values
if any(logical(elem(:,5) == 0))
    %Fill the fifth column of "elem" with "1". Thus the material can be
    %attributed.
    elem(:,5) = 1;
end  %End of IF
end