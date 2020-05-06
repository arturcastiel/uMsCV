function [coord,nnode] = getcoord(verifymshfile)

%Open the *.msh file
readmsh = fopen(verifymshfile);
%"nnode" is the number of nodes in the discrete domain
getmshdata = textscan(readmsh,'%u',1,'HeaderLines',4);
%Attribute the data to "nnode"
nnode = getmshdata{1};

%"coord" is a matrix which contain the coordinate of each point of domain
getmshdata = textscan(readmsh,'%*u %f64 %f64 %f64',nnode,'HeaderLines',1);
%Fill the matrix "coord" with the x, y and z coordinates.
coord = cell2mat(getmshdata);

%Clear "getmshdata"
clear getmshdata;

%Close the *.msh file
fclose(readmsh);

end