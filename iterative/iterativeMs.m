function [p_last] = iterativeMs(TransF, F, ac,bc, OP, OR, po)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global superFolder
index = 1;
maxoutit = 7;
maxinit = 500;
maxoutit = 10;
maxinit = 2000;

tol = 1e-3;

tol = 1e-3;
in_tol = 0.8e-3;



tol =  1e-2;
in_tol =  2 * 1e-3;


%%
tol =  1e-1;
in_tol =   1e-4;
%%
tol =  1e-1;
in_tol =  1 * 1e-2;

 tol =  1e-2;
 in_tol =   1e-4;

%in_tol = 0.3*1e-2;
%in_tol = 0.3*1e-4;
%symrcm(TransF)

%ameba
 tol =  1e-3;
 in_tol =   1e-6;
 
 % SPE
  tol =  1e-2;
 in_tol =   1e-4;

 
  % SPE
  tol =  1e-3;
 in_tol =   1e-6;
 
 
 
% in_tol =   1e-3;
%  tol =  1e-2;
%  in_tol =   1e-4;
p_old = po;

%ac1 = OP'*TransF*OP;
%ac1i = ac1^-1;
aci = ac^-1;

[L,U] = ilu(TransF,struct('type','nofill','droptol',1e-5));
%[L,U] = ilu(TransF);
%OR = OP';

%tol = 0.5*1e-6;

r_old = F - (TransF*p_old);

dlmwrite(strcat(superFolder,'\res.dat'),[]);
dlmwrite(strcat(superFolder,'\control.dat'),[]);

while index <= maxoutit  & (norm(r_old)> tol)  %(max(abs(r_old))> tol)
    disp(['Iterativo:' num2str(index)])
    
    %dp1 = OP* aci * OR'*r_old;
    
   % dp1 = OP * (aci)* OR*r_old;
    %dp1 = OP * (aci)* OP'*r_old;
    dp1 = OP * (aci)* OR*r_old;

   %dp1 = OP*((OR*r_old)'\ac );  
    
    
    rp1 = r_old - TransF*dp1 ;
    
       
  %[p_old,fl1,rr1,it1,rv1]=bicgstab(M_old,RHS_old,1e-10,1000,L,U);
   %[dp2,fl1,rr1,it1,rv1]=gmres(TransF,rp1,10,1e-4,20,L,U, rp1);

    [dp2,flag,relres,iter,resvec] = bicgstab(TransF,rp1,in_tol,maxinit,L,U, rp1);
    dlmwrite(strcat(superFolder,'\res.dat'),resvec,'-append');
    dlmwrite(strcat(superFolder,'\control.dat'),[iter,size(resvec,1)],'-append');
    disp('suavizacao')
    disp(iter)
    
    %dp2 = TransF\rp1;
    %spparms('spumoni',2)
    p_old = p_old + dp1 + dp2;
    %p_old = p_old +  dp1;
    r_old = F - (TransF*p_old);
    index = index + 1;
end

1
%
r_old = F - (TransF*p_old);
%dp1 = OP* ac * OP'*r_old;
%dp1 = OP * (aci)* OR*r_old;
dp1 = OP * (aci)* (OP')*r_old;

p_old = p_old + dp1;


%     
% 
% 
% r_old = F - (TransF*p_old);
% %dp1 = OP* ac * OP'*r_old;
% dp1 = OP * (ac^-1)* OR*r_old;
% p_old = p_old + dp1;
%      
     
     %     
     
% %     
% %     
% %     rp1 = r_old - TransF*dp1 ;
% %    
%     
%        
%   %[p_old,fl1,rr1,it1,rv1]=bicgstab(M_old,RHS_old,1e-10,1000,L,U);
%     %[dp2,fl1,rr1,it1,rv1]=gmres(TransF,rp1,10,1e-4,1000,L,U);
%     
%     dp2 = TransF\rp1;
%     p_old = p_old + dp1 + dp2;
%     index = index + 1;
%     
    


p_last = p_old;
end

