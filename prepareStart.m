if strcmp(osMode,'windows')
    superFolder = strcat(pwd,'\start\');
    command = ['copy ' superFolder nameFile ' ' pwd '\start.dat'];
    system(command);
elseif strcmp(osMode,'linux')
    superFolder = strcat(pwd,'/start/');
    command = ['cp ' superFolder nameFile ' ' pwd '/start.dat'];
    system(command);
end