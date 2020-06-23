function [ pms ] = neumanm(A,B, coarseelem , edgesOnCoarseBoundary, fluxMs,pc )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global refCenterInCoaseElem npar inedge elemloc coarseelem dictionary coarseDiricht coarseElemCenter

pms = sparse(size(B,1),1);
%%Adicionar a influência em ambos os lados do fluxo



%NAO SE USA DUAS CONDICOES DE DIRICHLET 
% A CONDICAO DE DIRICHLET EH IMPOSTA AUTOMATICAMENTE
%PORTANTO DEVE-SE RETIRAR A IMPOSICAO DA CONDICAO DE DIRICHLET NOS VERTICES
%NESTES ELEMENTOS
%%Separar em matrizes menores e resolver o sistema

for index = 1:npar    
    if index == 16
        1
    end
    ce=coarseelem{index};
    peqA = A(coarseelem{index},:);
    peqA = peqA(:,coarseelem{index});
    peqB = B(coarseelem{index});

   % if sum(coarseDiricht == index) == 0 & (coarseDiricht ~= 0)
    if ~coarseDiricht(index)
        %impondo dirichlet se os elementos nao estiverem na fronteira ja
        %com valor de dirichlet. problema ja bem posto
        peqA(refCenterInCoaseElem(index),:) = 0;
        peqA(refCenterInCoaseElem(index),refCenterInCoaseElem(index)) = 1;
        
        %somar os fluxos em peqB antes de rodar para testar se esta de fato
        %conservativo, boa maneira de debugar
        
        peqB(refCenterInCoaseElem(index)) = pc(ce(refCenterInCoaseElem(index)));
    end
    
    
    peqSol =   peqA\peqB;

    pms(coarseelem{index}) = peqSol;
end

%%Adicionar as influências das vazoes calculadas

%%Resolver os sistemas menores


%% Devolver a resposta para um mesmo vetor e retornar

end

