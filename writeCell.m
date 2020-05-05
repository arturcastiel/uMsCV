function writeCell(A,fileName)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
fid = fopen(fileName, 'w');
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid,'%s', A{i});
        break
    else
        fprintf(fid,'%s\n', A{i});
    end
end
fclose(fid);
end

