function [ M, I ] = globalmatrixmpfad_bc( M,I)

global coord elem esurn1 esurn2  bedge inedge  centelem elemarea bcflag wells
% auxmobility1=mobility(1:size(inedge,1),1);
% auxmobility2=mobility((size(inedge,1)+1):(size(inedge,1)+size(bedge,1)),1);
% mobility(1:size(bedge,1),1)=auxmobility2;
% mobility((size(bedge,1)+1):(size(inedge,1)+size(bedge,1)),1)=auxmobility1;
%-----------------------inicio da rutina ----------------------------------%
%Constr�i a matriz global.
gamma = 1e+10;

% % adequa��o da matriz nos po�os produtores
if max(wells)~=0
    for iw = 1:size(wells,1)
        if wells(iw,2)==2 %produtor
            M(wells(iw,1),:)=0*M(wells(iw,1),:);
            M(wells(iw,1),wells(iw,1))=1;
            I(wells(iw,1))=0;
        end
    end
end

% %% Adding the influence of the wells on the matrix
if size(wells,2) > 1
    wellSolv;
    if ~isempty(wells)
        %Saturation flags
        for index = 1:size(flagsSatElem,1)
            I(flagsSatElem{index})=  (flagsSatValue(index)/flagsSatTotalArea(index))  * elemarea(flagsSatElem{index});
        end
        
        for index = 1: size(flagsInjElem,1)
            M(flagsInjElem{index},:) = 0* M(flagsInjElem{index},:);
            %flagsInjElem{index},:)
            for ii = 1 : size(flagsInjElem{index},2)
                M(flagsInjElem{index}(ii),flagsInjElem{index}(ii)) = 1;
                I(flagsInjElem{index}(ii)) = flagsInjValue(index);
            end
            %for ii = 1: size(flags
        end
        
        %         for index = 1: size(flagsProdElem,1)
        %             %M(flagsProdElem{index},:) =  M(flagsProdElem{index},:) + gamma * ones(  ;
        %
        %             %flagsInjElem{index},:)
        %             for ii = 1 : size(flagsProdElem{index},2)
        %                 M(flagsProdElem{index}(ii),flagsProdElem{index}(ii)) = M(flagsProdElem{index}(ii),flagsProdElem{index}(ii)) + gamma;
        %                 I(flagsProdElem{index}(ii)) = I(flagsProdElem{index}(ii))  + gamma * flagsProdValue(index);
        %             end
        %             %for ii = 1: size(flags
        %         end
        
        
        for index = 1: size(flagsProdElem,1)
            M(flagsProdElem{index},:) = 0* M(flagsProdElem{index},:);
            %flagsInjElem{index},:)
            for ii = 1 : size(flagsProdElem{index},2)
                M(flagsProdElem{index}(ii),flagsProdElem{index}(ii)) = 1;
                I(flagsProdElem{index}(ii)) = flagsProdValue(index);
            end
            %for ii = 1: size(flags
        end
    end
end



end

