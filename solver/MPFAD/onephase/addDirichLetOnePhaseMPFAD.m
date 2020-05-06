function [ M, I ] = addDirichLetOnePhaseMPFAD(M,I)

global wells elemarea 
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%% Adding the influence of the wells on the matrix
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

