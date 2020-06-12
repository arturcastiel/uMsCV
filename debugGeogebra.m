

QQ = pTransF;

[x,y] = find(QQ ~=0);
mat = [x,-y];


% n = size(QQ,1);
% 
% mat = zeros(sum(sum(ref)), 2);
% index = 1;
% for x = 1:n
%    for y = 1:n
%        if ref(x,y)
%           mat(index,:)= [x,y];
%           index = index+1;
%        end
%    end    
% end