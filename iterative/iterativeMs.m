function [p_last] = iterativeMs(TransF, F, ac,bc, OP, OR, po)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
index = 1;
%symrcm(TransF)



p_old = po;

ac1 = OP'*TransF*OP;
ac1i = ac1^-1;
aci = ac^-1;

[L,U] = ilu(TransF,struct('type','nofill','droptol',1e-4));

%OR = OP';
tol = 0.5*1e-3;
r_old = F - (TransF*p_old);
while index <= 5 & (max(abs(r_old))> tol)

    
    %dp1 = OP* aci * OR'*r_old;
    
    dp1 = OP * (ac^-1)* OR*r_old;
    %dp1 = OP*((OR*r_old)'\ac );  
    
    
    rp1 = r_old - TransF*dp1 ;
   
    
       
  %[p_old,fl1,rr1,it1,rv1]=bicgstab(M_old,RHS_old,1e-10,1000,L,U);
   %[dp2,fl1,rr1,it1,rv1]=gmres(TransF,rp1,10,1e-4,8000,L,U, rp1);
    dp2 = bicgstabl(TransF,rp1,tol,200,L,U, rp1);
   
    %dp2 = TransF\rp1;
    %spparms('spumoni',2)
    p_old = p_old + dp1 + dp2;
    %p_old = p_old +  dp2;
    r_old = F - (TransF*p_old);
    index = index + 1;
end
% 
r_old = F - (TransF*p_old);
%dp1 = OP* ac * OP'*r_old;
dp1 = OP * (ac^-1)* OR*r_old;
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

