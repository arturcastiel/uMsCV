function [elem] = getelemf(verifymshfile, nnode)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
global osMode
%Open the *.msh file
if strcmp(osMode, 'windows')
    fS = '\';
elseif strcmp(osMode, 'linux')
    fS = '/';
end
verifymshfile = [pwd fS 'grids' fS verifymshfile];
A = openCell(verifymshfile);

flag = 0;
for index = nnode+9 : (size(A,2) - 2)
    tline = str2num(A{index});
    if (flag == 0) & ((tline(2) ~= 15) &  (tline(2) ~= 1))
        flag = index;
    end
end
elem = zeros( ((size(A,2) - 2) - flag) + 1, 4);
% ii = 1;
% for index = flag : (size(A,2) - 2)
%     auxvec = str2num(A{index});
%     elem(ii,:) = auxvec(end-3:end);
%     ii = ii + 1;
% end
% if all(elem(:,1) == 1)
%    elem(:,1:3) =   elem(:,2:4);
%    elem(:,4) = 0;    
% end
A = A(flag:end-2);
%type_flag = false(size(A,2),1);
for ii = 1:size(A,2)
   line = str2num(A{ii}); 
   if line(2) == 2
       elem(ii,1:3) = line(end-2:end);
   elseif line(2) == 3
       elem(ii,1:4) = line(end-3:end);
   end    
end

end