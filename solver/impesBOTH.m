%% VPI Snapshots
vpi_vecF = [ 0.02 :0.05:1];

vpi_vecF = [ 0.02 :0.05:1];
% %%
% %plot perm
% postprocessorPerm(elem(:,end),0,superFolder)
%Write header
header
%Info on the mesh
meshProp

%Write Info on Coarse Mesh
coarseVTK

%% copying start.dat and msh to the result folder
copyfile('start.dat',superFolder)
strr = strcat(keypath1{1},'/',mshfile{1});
copyfile(strr,superFolder)
%% copying iterativeRoutine and msh to the result folder
copyfile('iterative/iterativeRoutine.m',superFolder)

% % % %%  ms MPFAD
  %  vpi_vec = vpi_vecF;
%  impesMPFAD2


%% starting simulation
% % regular MPFAD
%  vpi_vec = vpi_vecF;
 impes

%% starting simulation
% % regular TPFA
  %vpi_vec = vpi_vecF;
  %impesTPFA

%% 
readHeader
%% 
