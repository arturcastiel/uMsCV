function[N,F,V,S_old,S_cont]=presaturationMod(wells)

%depois de V tirei %weightLS%
global satlimit elem bedge bcflag
% adequa��o das faces por elemento e por n� em sentido anti-horario
% F: vetor de faces em cada elemento
% V: faces ordenados "convenientemente" ao rededor de um n�
% N: faces ordenados ao rededor de um no em sentido anti-horario

[F,V,N]=elementface;

%[esuel, esuel_coord,A,bound] = esurnelem;

% calculo dos pesos para o m�todo minimos quadrados
%[weightLS] = weights(A,esuel,esuel_coord);

% Condi�ao inicial satura�ao
S_old = zeros(size(elem,1),1);

S_old(:)=satlimit(2);
% Condi�ao de contorno da satura�ao

S_cont=1-satlimit(1);

% adequa��o dos po�os
if max(wells)~=0
        %% minha modificacao
        if  wells(iw,end-1) <= 500   
            S_old(wells(iw,1)) = 1;
        end
    
    
%     
%     for i=1:size(wells,1)
%         if wells(i,2)==1
%             S_old(wells(i,1))=S_cont;
%         end
%     end
end

    for ifacont=1:size(bedge,1)
        x=bcflag(:,1)==bedge(ifacont,5);
        r=find(x==1);
        rr=bcflag(r,2);
        if rr~=0 && bcflag(r,1)>200
            S_cont=1-satlimit(1);
        end
        
    end

end