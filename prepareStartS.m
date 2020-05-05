% finding case
if strcmp(osMode,'windows')
    superFile = strcat(pwd,['\start\' nameFile]);
elseif strcmp(osMode,'linux')
    superFile = strcat(pwd,['/start/' nameFile]);
end
A = openCell(superFile);
% converting case to os
A = convertInput(A,0);
%copying case to folder
if strcmp(osMode,'windows')
    sDir = split(superFile,'\');
    sFolder = [pwd '\' sDir{end-1} '\' sDir{end}(1:end-4) '-windows.dat'];
elseif strcmp(osMode,'linux')
    sDir = split(superFile,'\');
    sFolder = [pwd '/' sDir{end-1} '/' sDir{end}(1:end-4) '-linux.dat'];
end
writeCell(A,sFolder);
if strcmp(osMode,'windows')
    superFile = strcat(pwd,['\start\' [nameFile(1:end-4) '-windows.dat']]);
    command = ['copy ' superFile ' ' pwd '\start.dat'];
elseif strcmp(osMode,'linux')
    superFile = strcat(pwd,['/start/'  [nameFile(1:end-4) '-linux.dat']]);
    command = ['copy ' superFile ' ' pwd '/start.dat'];
end
system(command);
