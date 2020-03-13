%% IMPES
% Implicit Pressure , Explicit Saturation


%% Saturation Preprocessor 
[N,Fo,V,S_old,S_cont]=presaturation(wells);


%%Calculating Fixed Parameters
[Hesq, Kde, Kn, Kt, Ded] = Kde_Ded_Kt_Kn(kmap);
fonte=0;
nflagno= contflagnoN(bedge);
[w,s]=Pre_LPEW_2(kmap,N);





%Initializing IMPES variables
CFL = courant;
cont = 1;
vpi_old = 0;
VPI(1) = 0;
cumulateoil(1) = 0;
oilrecovery(1) = 1;
watercut(1) = 0;
vpi = totaltime(2);
t_old = 0;
step=0;
time2=0;

%% exponente das permeabilidades relativas cuidado...! 
nw=2;
no=2;

%% este flag s� use quando o problema � Buckley-Leverett com fluxo imposto
% na face com flag 202, por�m a satura��o ser� imposto na face
auxflag=[0];

while vpi_old < vpi %true
    
    if cont == 180
        break
    end
    [mobility] = mobilityface(S_old,nw,no,auxflag,S_cont);
    [p,errorelativo,flowrate,flowresult]=solverpressure(kmap,fonte,...
        mobility,wells,S_old,V,nw,no,N,auxflag,Hesq, Kde, Kn, Kt, Ded,nflagno);
    
    
    
    %% calculo do fluxo fracional em cada elemento da malha
    f_elem = fractionalflow(S_old,nw,no);
    
    %% caculo do passo de tempo
    % order sempre = 1
    d_t=steptime(f_elem,S_old,flowrate,CFL,S_cont,auxflag,nw,no); %,order);
    
    %% calculo da satura��o explicito
    %nao usar - isso
    %esuel1 weightLS,bound,upsilon,kappa, smetodo,tordem,
%     S_old=solversaturation(S_old,flowrate,d_t,esuel1,wells,flowresult,...
%         f_elem,S_cont,weightLS,bound,smetodo,tordem,...
%         upsilon,kappa,nw,no,auxflag);
    
        [S_old]= firstorderstandard(S_old,flowrate,flowresult,f_elem,d_t,wells,S_cont,auxflag,nw,no);
        
        X = sprintf('Calculo do campo de satura��o pelo m�todo: %s ','FOU');
        disp(X)
        
    
    %% reporte de produ��o
    [VPI,oilrecovery,cumulateoil,watercut]=reportproduction(S_old,...
        wells,f_elem,cont,VPI,oilrecovery,cumulateoil,watercut,flowresult,d_t);
    
    %% calculo o passo de vpi ou passo de tempo
    %t_old=VPI
    
    vpi_old=VPI(cont);
    
    %visualiza��o
    if vpi_old >= time2
        step = 1000*time2;
        postprocessorT(p,S_old,step);
        time2 = time2 + 0.01;
    end
%     % armezena a produ��o em arquivo .dat
%     fprintf(fid3,'%13.3e %13.3e %13.3e %13.3e \n',VPI(cont), ...
%         oilrecovery(cont),cumulateoil(cont),watercut(cont));
%% mudei aqui  
%  vpi_old=1e+40;
    %postprocessor(p,S_old,cont+1)
    cont=cont+1
    
end
