%FV-MsRB - MultiScale Restriction smoothed Basis method
%2d -Impl�cita na Press�o Expl�cita na Satura��o (IMPES)
% 21/09/2016 por Artur Castiel Reis de Souza

if phasekey == 1 && strcmp(pmethod, 'tpfa')
    onephaseflowTPFA
elseif phasekey == 1 && strcmp(pmethod, 'mpfad')
    onephaseflowMPFAD
end