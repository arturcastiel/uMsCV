%% Preparing output dat
global superFolder osMode

if strcmp(osMode,'windows')
    sF = '\';
elseif strcmp(osMode,'linux')
    sF = '';
end
caseType = mshfile{1}(1:end-4);
folder = [sF 'results'];
name1 = [sF 'Original.dat'];
name2 = [sF 'Multiscale.dat'];
name3 = [sF 'Header.dat'];
name4 = [sF 'CurveOriginal.dat'];
name5 = [sF 'CurveMs.dat'];
if itOn == 0
    corText = '-offCor';
else
    corText = '-onncor';
end
caseTypeT = strcat(sF ,caseType,'-',num2str(npar),'cv',corText);
superFolder = strcat(pwd,folder,caseTypeT);


%pointers to files
fileID1 = strcat(superFolder,name1);
fileID2 = strcat(superFolder,name2);
fileID3 = strcat(superFolder,name3);
fileID4 = strcat(superFolder,name4);
fileID5 = strcat(superFolder,name5);

fid3 = fopen(fileID4,'w');
fid4 = fopen(fileID5,'w');
mkdir(superFolder);