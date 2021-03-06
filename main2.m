%FV-MsRB - MultiScale Restriction smoothed Basis method
%Two Phase Flow
%2d -Impl�cita na Press�o Expl�cita na Satura��o (IMPES)
% 21/09/2016 por Artur Castiel Reis de Souza

%% Main File Start Up
%  Declare all global variables to be used in the code
global tol_c coord centelem elem esurn1 esurn2 nsurn1 nsurn2 bedge inedge ...
    normals esureface1 esureface2 esurefull1 esurefull2 elemarea dens visc...
    satlimit pormap bcflag courant totaltime numcase oneNodeEdges keypath2 foldername...
    elemloc npar coarseelem  exinterface exinterfaceaxes ghostelem pointWeight bedgeNode ...
    numinterface interfacecenter  coarseblockcenter coarseneigh intCoord flagboundcoarse ...
    intRegion   boundRegion GlobalBoundary H outSupport coarseElemCenter ...
    coarseningRatio wells mshfile edgesOnCoarseBoundary refCenterInCoaseElem ...
    dictionary edgesCoarseDict coarseDiricht intinterface pointloc regularEdges semiEdges ...
    coarseedge ordem splitFag bold nc pt flagSuavizador suavizador MetodoSuavizador
suavizador = 'SOR';
MetodoSuavizador =  'S_Multiescala';
flagSuavizador = 'on';
nameFile = 'start_fivespot.dat';
%nameFile = 'start_linear.dat';
%nameFile = 'start_ameba.dat';
%nameFile = 'start_serra.dat';
pMethod = 0;
splitFlag = 0;
prepareStart
multiscale = 'on';
smetodo = 'FOU'
pt = 3;
nc = 12;
%Globals2D_CPR;
global osMode 
osMode = 'linux';


%% Tolerance Settings
tol_c = 0.00001; flagboundcoarse = 0;

%% Adds Preprocessor folder to the path
% Including other folders to the path
path(path,'preprocessor');

%% Including the folder of the preprocessor of the method MsRB
path (path,'preMsRSB');

%% Including the folder of the postprocessor of the method MsRB
path (path,'postprocessor');

%% Including the folder of the multiscale library
path (path,'MsRB');
path (path,'preMsMPFAD');

%% Including Non-Linear Folder
% path(path,'nlfv_lpew');
% path(path,'nonlinear');

%% Including the folder of the solver and it subfolders
path (path,'solver');
if strcmp(osMode,'windows')
    path (path,'solver\TPFA');
    path (path,'solver\MPFAD\onephase');    
    path (path,'solver\MPFAD\onephase\preMPFAD');
elseif strcmp(osMode,'linux')
    path (path,'solver/TPFA');
    path (path,'solver/MPFAD/onephase');    
    path (path,'solver/MPFAD/onephase/preMPFAD');
end    
    %preMsMPFAD
%% Adds the path of the fuctions used for debbuging
path(path,'debug');

%% Adds the path of the fuctions used for creating complex permeability fields
path(path,'permeability');
path(path,'permfield');

%% Adds the path of class 
path(path,'class');

%% Adds the path of config options
path(path,'config');
colormat = load('color.dat');

%% Adds the path of SPE options
path(path,'SPE');

%% Adds the path Post Multiscale treatments
path(path,'postMs');

%% Adds the path of the Iterative Treatment
path(path,'iterative');


%% ====================== Preprocesador Global ============================
    [coord,centelem,elem,esurn1,esurn2,nsurn1,nsurn2,bedge,inedge,...
    normals,esureface1,esureface2,esurefull1,esurefull2,elemarea,dens,visc,...
    satlimit,pormap,bcflag,courant,totaltime,numcase,phasekey,kmap,...
    wells,elemloc,npar,coarseelem, ghostelem, multiCC,pmethod,mshfile, keymsfv,interptype,keypath1 ] = preprocessor;
% %% primal creation
% Q = [centelem(:,1:2), ones(size(elem,1),1)];
% 'cosine' 'sqeuclidean' 'correlation' 'cityblock'
% n = 35;
% mat = rand(n,2);
% mat = [mat, ones(size(mat,1),1)];
% seed = randi(size(centelem,1), n,1)

%seed = [1, 10, 900, 100, 5];

value = Q(seed,:);

F = kmeans(Q,n) %,'Distance','cosine', 'Start',Q(seed,:));
% 
% L =  kmeans(Q(F==1,:) ,4, 'Distance','cityblock');
% 
% F(F==1) = (L);
% %F = clusterdata(Q,'Linkage','centroid','MaxClust',n);
% 
% npar = max(F);
elemloc = F;
%% ===================== Preprocessador Multiescala =======================
%debugTest3;
%debugPoint;

multiCC = 1;
splitFlag = 1;


 [ elemloc, npar,coarseelem, coarseedge,intinterface,exinterface,exinterfaceaxes,...
    numinterface,interfacecenter, coarseblockcenter, coarseneigh, intCoord, multiCC, ...
    coarseningRatio, semiEdges,bedgeNode, edgesOnCoarseBoundary,oneNodeEdges, ...
    pointloc,pointWeight,regularEdges] = prems(true, splitFlag, pMethod, true, multiCC);
    if strcmp(multiscale, 'off')
        semiEdges = [];
    end
%% ================ Multiscale Preprocessador for MsRSB ===================
if ~isempty(semiEdges)
bold = 1;    
[ intRegion ,  boundRegion, GlobalBoundary, H, outSupport, ...
    coarseElemCenter,refCenterInCoaseElem, dictionary,edgesCoarseDict,coarseDiricht]   = preMsRB(npar,coarseneigh, centelem,coarseelem, ...
    coarseblockcenter,exinterface,multiCC);
%% Alternative Multiscale Preprocessor 
% tic;
bold = 2;
 [coarseElemCenter, coarse_interface_center, coarse_strips, boundRegion, intRegion, GlobalBoundary, H, outSupport, refCenterInCoaseElem, ...
     dictionary,edgesCoarseDict,coarseDiricht, cop, cqb] = alpreMsRB(npar,coarseelem, coarseneigh, centelem, exinterface, multiCC, splitFlag);

 % toc;

end
%% ======================== CPR Preprocessor ==============================
% % if strcmp(smetodo, 'cpr')
% % P=ordem-1; pre_CPR = meshQuad(P,elem,coord);
% % kxi=pre_CPR.kxi; eta=pre_CPR.eta;%nodes for global velocity field
% % EToE=pre_CPR.EtoE; EToV=pre_CPR.EtoV; EToF=pre_CPR.EtoF;
% % x=pre_CPR.x; y=pre_CPR.y;
% % dx = max(max(x(:,2:end)-x(:,1:end-1))); dy = max(max(y(:,2:end)-y(:,1:end-1)));
% % nE = pre_CPR.nE; npf = pre_CPR.eNodesPerFace*pre_CPR.eFaces;
% % nN = pre_CPR.eNodes; Nfp=pre_CPR.eNodesPerFace; Fmask = pre_CPR.FaceMask;
% % if nE<=12, figure;meshQuad.plotVertices(pre_CPR.EtoV,pre_CPR.VX,pre_CPR.VY);title('Mesh');axis square; end 
% % Create comunication and boundary maps
% % map = meshQuad.BuildMaps2Dv2(pre_CPR);
% % 
% % Load DG tools
% % NDG = DGtoolsQuad(pre_CPR);
% %     Vnd=NDG.V2D; J=NDG.J;
% %     Dr = NDG.D_kxi;   Ds = NDG.D_eta;
% %     rx = NDG.kxi_x;   sx = NDG.eta_x;
% %     ry = NDG.kxi_y;   sy = NDG.eta_y;
% % Build surface normals
% % [nx,ny,sJ] = DGtoolsQuad.Normals(NDG); 
% % Fscale = sJ./J(Fmask(:),:); 
% % disp('CPR preprocessing successfuly executed!')
% % end
%==========================================================================
%% superfolder maker
itOn = 0;
superFolderMK
% malha 32x32

%% distorcao de malhas estruturadas
%[auxcoord]=distortedramd;
%% s� em malha com furo
% x=bedge(129:144,1);
% y=bedge(129:144,2);
% bedge(129:144,1)=y;
% bedge(129:144,2)=x;
%% s� em malha "Tipo1malha1"
%   x=bedge(:,1);
%   y=bedge(:,2);
%   bedge(:,1)=y;
%   bedge(:,2)=x;
%   x1=elem(:,1);
%   x2=elem(:,3);
%   elem(:,1)=x2;
%   elem(:,3)=x1;
%% S� malha "drainoblique1" no entanto revise outras malhas deste tipo
%    x=bedge(:,1);
%    y=bedge(:,2);
%    bedge(:,1)=y;
%    bedge(:,2)=x;
%% use se deseja resolver o problema de Queiroz et al 2014
% unstructured mesh com furo reordenando o sentido da fronterira no
% contorno interior
%  x=bedge(135:150,1);
%  y=bedge(135:150,2);
%  bedge(135:150,1)=y;
%  bedge(135:150,2)=x;
%  bedge(135:150,4:5)=102; % 36x36 cuidado
%  bcflag(2,1)=102;
%  bcflag(2,2)=2;
%-----------------------------
% unstructured mesh com furo reordenando o sentido da fronterira no
% contorno interior
%  x=bedge(71:78,1);
%  y=bedge(71:78,2);
%  bedge(71:78,1)=y;
%  bedge(71:78,2)=x;
%  bedge(71:78,4:5)=102; % benchmark23_3_GROSSO 18x18
%  bcflag(2,1)=102;
%  bcflag(2,2)=2;

% unstructured mesh com furo reordenando o sentido da fronterira no
% contorno interior
% x=bedge(289:320,1);
% y=bedge(289:320,2);
% bedge(289:320,1)=y;
% bedge(289:320,2)=x;
% bedge(289:320,4:5)=102; % benchmark23_3 72x72
% bcflag(2,1)=102;
% bcflag(2,2)=2;

% unstructured mesh com furo reordenando o sentido da fronterira no
% contorno interior
%  x=bedge(577:640,1);
%  y=bedge(577:640,2);
%  bedge(577:640,1)=y;
%  bedge(577:640,2)=x;
%  bedge(577:640,4:5)=102; % benchmark23_3 144x144
%  bcflag(2,1)=102;
%  bcflag(2,2)=2;
%% Modifica��o Malha Kershaw
%bedge(:,4:5)=101;
%----------------------------
%% tratamento malha Hermeline
%bedge(:,4:5)=101;
% malha 16x16
%x=bedge(16:24,1);
%y=bedge(16:24,2);
%bedge(16:24,1)=y;
%bedge(16:24,2)=x;
% malha 32x32
% x=bedge(33:48,1);
% y=bedge(33:48,2);
% bedge(33:48,1)=y;
% bedge(33:48,2)=x;
% malha 64x64
% x=bedge(64:96,1);
% y=bedge(64:96,2);
% bedge(64:96,1)=y;
% bedge(64:96,2)=x;
% malha 128x128
%x=bedge(128:192,1);
%y=bedge(128:192,2);
%bedge(128:192,1)=y;
%bedge(128:192,2)=x;
%%
tic
%% escolha o m�todo de interpola��o
% 1-->LPEW1 Gao e Wu 2010
% 2-->LPEW2 Gao e Wu 2010
interptype=2;
%% Digite
% monofasico ---> quando deseje rodar um problema de escoamento monof�sico ou
% bifasico   ---> quando deseja rodar um problema de ecoamento bif�sico ;
% quando rodar bif�sico verifique linha 46-50 do Ks_Interp_LPEW2.m
simu='monofasico';
%% escolha o tipo de erro discreto que deseja usar
% erromethod1 ---> erro utilizado por Gao e Wu 2010
% erromethod2 --->  ''     ''     por Lipnikov et al 2010
% erromethod3 --->  ''     ''     por Eigestad et al 2005
% erromethod4 --->  ''     ''     por Shen e Yuan 2015
erromethod='erromethod1';
%% Defina o tipo de solver de press�o
% tpfa      --> m�todo Linear dos volumes finito TPFA
% mpfad     --> (MPFA-D)--> Alvaro
% lfvLPEW   --> m�todo linear basedo no m�todo n�o linear usando LPEW (MPFA-HD), ter cuidado linha 52 e 54 do preNLFV
% lfvHP     --> (MPFA-H)--> Allan
% lfvEB     --> m�todo completamente baseado na face (MPFA-BE).
% nlfvLPEW  --> (NLFV-PP)--> Abdiel %ESTAVA ESSE
% nlfvDMPSY --> m�todo n�o linear que preserva DMP baseado no artigo (Gao e Wu, 2013) e (Sheng e Yuan, 20...)
% nlfvIntFree --> M�todo de volume finito livre de interpla��o
pmetodo='mpfad';
%% metodo de intera��o: iterpicard, iternewton, iterbroyden, itersecant,
% m�todo de itera��o pr�prio de m�todos n�oo lineares
%iterfreejacobian,iterdiscretnewton, JFNK
%iteration='iterdiscretnewton';
%iteration='iterbroyden';
iteration='iterpicard';
%iteration='iterhybrid';

%% digite segundo o benchmark
benchmark='gaowu9';


%% adequacao de flag
%adeSPE
%verficar a alteracao em wells no preprocessador
%wells = [6490, 1,  301,  1,  0, 1;wells]
%adeqPerm
%% Multiscale Solver

meshAnalysis
%permChannel
%changepermeability 
%datageneration
%solver;