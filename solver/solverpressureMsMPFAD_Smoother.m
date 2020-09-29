function [pressure,errorelativo,flowrate,flowresult,OP_old,b2,b3,tempo,v]=solverpressureMsMPFAD_Smoother(kmap,fonte,...
    mobility,wells,S_old,V,nw,no,N,auxflag,Hesq, Kde, Kn, Kt, Ded,nflag,OP_old,S_cont)
global w s flagboundcoarse superFolder pt nc npar coarseElemCenter refDir
% function [pressure,errorelativo,flowrate,flowresult]=solverpressure(kmap,nflagface,nflagno,fonte,...
%     tol, nit,p_old,mobility,gamma,wells,S_old,V,nw,no,N,parameter,metodoP,...
%     auxflag,interptype,Heseyh v, Kde, Kn, Kt, Ded,weightDMP,auxface,benchmark,iteration,nflag)
% TO DO
% adaptar calculo PRE_LPEW_2M para contar com mobilidade
% adaptar calculo dos fluxos globais e fluxos do inside coarse mesh para
% contar com mobilidade
global pointWeight elem coarseelem  edgesOnCoarseBoundary npar bedge oneNodeEdges flagSuavizador ...
       suavizador MetodoSuavizador Wc Wf bold
errorelativo=0;

wsdynamic = dynamic(pointWeight);


[w,s]=Pre_LPEW_2T(kmap,mobility,V,S_old,nw,no,N);
mobRegion = mobilityfaceRegion(S_old,nw,no,auxflag,S_cont,mobility);
%assembly da matriz
%[ TransF, F] = globalmatrixmpfadn( w,s, Kde, Ded, Kn, Kt, nflag, Hesq,wells,mobility,fonte);
[ TransF, F] = globalmatrixmpfadnobc( w,s, Kde, Ded, Kn, Kt, nflag, Hesq,wells,mobility,fonte);

%[perm_matrix, per_vec] = genPermMatrix();
%%MultiScale
%pre condition matrix 
% TransFc = TransF; %original
% n = size(TransFc,2);
% TransFc(1:size(TransFc,2)+1:end) = diag(TransFc) - sum(TransFc,2);
auxflag = 0;
% %OP_old
%OR  = genRestrictionOperator(size(elem,1), npar);

OR  = genRestrictionOperator();

% if OP_old == -1;
%     %OP =  genProlongationOperator(OR', TransFc, 2/3, 1800);
%     [OP,CT] = genProlongationOperatorAMS(TransF, F);
% 
% else
%     %OP =  genProlongationOperator(OP_old, TransFc, 2/3,300); 
%     [OP,CT] = genProlongationOperatorAMS(TransF, F);
% end
%
 if size(wells,2) > 1
     ref = wells(:,5) > 400;
end

%load('file.mat','AA','BB')
% disp('OP generation')
% tic
% toc
tic
[ TransF, F] = globalmatrixmpfad_bc(TransF, F);
%qq = TransF\F;
[ TransF, F] = removeDirichlet(TransF,F);
po = TransF\F;

%TransF(wells(ref,1),:) = 0;

[OP,CT] = genProlongationOperatorAMS(TransF, F);



%OP = bsxfun(@rdivide, OP(reff,:) ,sum( OP(reff,:) ,2));        





toc;

%
%OR(:,wells(ref,1)) = 0;


disp('OP generation')
tic



%OP(wells(ref,1),:) = 0;
toc
% for ii = find(refDir)'
%     Trans(ii,:) = 0;
%     Trans(ii, ii) = 1; 
%     F(ii) = 1;    
% end
po = TransF\F;

%postprocessorOP(OP,0,superFolder,'Operadores');
OP_old = OP;

% 
% %% wells treatment
% if size(wells,2) > 1
%     ref = wells(:,5) > 400;
%     % %testes
%      OP(wells(find(ref),1),:) = 0 ;
% end

%OP = bsxfun(@rdivide, OP(:,:) ,sum( OP(:,:) ,2)); 

 %% Precondicionador para resolver o sistema A*x=b
ac = OR * TransF * OP; 
bc = (OR * F) - OR*TransF *CT;
%% 


% % % % 
% % % if size(wells,2) > 1
% % %     ref = wells(:,5) > 400;
% % %     % %testes
% % %     OP(wells(find(ref),1),:) = 0 ;
% % % end

% for jj = wells(ref,1)'
%     refC = (coarseElemCenter == jj);
%     if sum(refC) ~= 0
%         refB = wells(:,1) == jj;
%         value = wells(refB,end);
%         pc(refC) = value;    
%     end
% end
%% Solu��o multiescala
pc = ac\bc;

% for jj = wells(ref,1)'
%     refC = (coarseElemCenter == jj);
%     if sum(refC) ~= 0
%         refB = wells(:,1) == jj;
%         value = wells(refB,end);
%         pc(refC) = value;    
%     end
% end
pd = OP*pc + CT;
%pd = OP*pc;
%pd = OP*pc;
if size(wells,2) > 1
    ref = wells(:,5) > 400;
    pd(wells(ref,1)) = wells(ref,end);
end


%[ pm] = multiscaleScheme(An,Bn,TransF, F)

%% Suavizador
disp('Iterativo')
tic 
pv = pd;

pd = iterativeMs(TransF, F, ac,bc, OP, OR, pd);

postprocessorTMS(full((pd)), full((pv)),0,superFolder,'pd - pv')
%pd = po
toc
postprocessorTMS(full((po)), full((po)),0,superFolder,'original')

postprocessorTMS(full((OP*pc)), full((CT)),0,superFolder,'OPpc - CT')
%pd = OP*pc;

disp([normError(pd, po ,2) , normError(pd, po ,inf)])
disp([normError(pv, po ,2) , normError(pv, po ,inf)])

%C = CT\F;
tempo = 0;
v = 0;
%% Calculates the flow on the Edges of the Coarse Boundary
%flowPd = flowrateMsMPFAD(edgesOnCoarseBoundary,pd,w,s,Kde,Ded,Kn,Kt,Hesq,nflagno,auxflag);
[flowPd, flowresultPd,velocityPd]=calflowrateMPFADn(pd,w,s,Kde,Ded,Kn,Kt,Hesq,nflag,auxflag,mobility);
flowPd2 = flowPd( edgesOnCoarseBoundary + size(bedge,1));
%[flowPd3, flowresultPd3,velocityPd3]=flowrateMPFAD(p,w,s,Kde,Ded,Kn,Kt,Hesq,nflagno,auxflag);
%poq = pd;
%iterativeRoutine
%% Neuman
%recalculating weights
%tests replacing flowPd by velocityPd
%alterar velocityPd


[wsdynamic] = Pre_LPEW_2MS(kmap, wsdynamic,velocityPd,mobRegion,S_old,nw,no);    


%excluding influence of edges and nodes that connect different coarse regions
[TransFn, Fn] = minusMPFAD(TransF,F, edgesOnCoarseBoundary, w,s, Kde, Ded, Kn, Kt, nflag, Hesq,mobility);
%debugUncoupling %tests TransFn to check if the uncoupling is done alright
%adding the influence back with weights that consider only the coarse
%volumes inside
%alterar flowPd2
%nneuman(TransF,F, pc, flowPd2)
[TransFn, Fn] = addMPFAD(TransFn,Fn, edgesOnCoarseBoundary, oneNodeEdges,flowPd2,wsdynamic,Kde, Ded, Kn, Kt, nflag, Hesq ,mobRegion ); % VERIFICAR DEPOIS!!!!
%debugUncoupling

% flowPd(abs(flowPd) < 0.00000001) = 0;
% [~,ref] = sort(coarseElemCenter)
%pc = pc(ref);



pp =  neumanmMPFAD(TransFn,Fn, coarseelem , edgesOnCoarseBoundary, flowPd,pd );
%pp =  neumanmMPFAD2(TransFn,Fn,pd );

[flowPp, flowresult,velocity]=flowratePPMPFAD(pp,w,s,wsdynamic, Kde,Ded,Kn,Kt,Hesq,nflag,auxflag,mobility,mobRegion);
% erro = max(abs(pd-pp)./max(pd))
% pause(1.5)

%superFolder = 'C:\Users\admin\Google Drive\Programa��o\FV-MsRB\';
%postprocessorName(full(pd),full(pp),superFolder, 'PP-Teste');

flowPms = compoundFlow( flowPp,flowPd2, size(bedge,1) );
%flowPms(abs(flowPms) < 0.00000001) = 0;
%postMPFA

flowrate = flowPms;
pressure = pd;

% p = TransF\F;
% [flowPn, flowresultPn,velocityPn]=calflowrateMPFAD(p,w,s,Kde,Ded,Kn,Kt,Hesq,nflag,auxflag,mobility);
%flowrate = flowPn;
%pressure = p;
% flowresult = flowresultPn;
% ooo = find(abs(flowresultPn) > 0.0000001);
flowresult = fluxSummation(flowPms);
consTestMpfa
% postprocessorTMS(full(pd),full(pp),0,superFolder,'Neummna');
% if bold == 1
%     postprocessorTMS(full(pd),full(po),0,superFolder,'Primeiro');
% elseif bold == 2
%     postprocessorTMS(full(pd),full(po),0,superFolder,'Segundo');
%end
disp([normError(pd, po ,2) , normError(pd, po ,inf)])
disp([normError(pv, po ,2) , normError(pv, po ,inf)])
ll = [npar pt nc bold kmap(end,end) normError(pd, po ,2) normError(pd, po ,inf)];
dlmwrite('simuSingle.txt', full(ll),'-append')

OP_old = OP;
end