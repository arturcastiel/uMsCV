function [parameter,p_old,tol,nit,nflagno]=preNLFV(kmap,benchmark,bedge)
global elem

%% calculo dos parametros ou constantes (ksi)
[parameter]=coefficientLPSangle(kmap);
% adequa��o dos flags de contorno
nflagno= contflagno(benchmark,bedge);
%% dados inicializa��o m�todos dos volumes finitos n�o linear
p_old=1e1*ones(size(elem,1),1); % inicializando a press�o
tol=1e-14;                        % tolerancia
nit=10000;                       % # itera��es de Picard
end