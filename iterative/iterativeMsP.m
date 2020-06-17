function [p_last] = iterativeMsP(TransF, F, ac,bc, OP, OR, po)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
index = 1;
lperm = createPerMat(symrcm(TransF));


lTransF = lperm * TransF * lperm';

lF = lperm * F;

lOP = lperm * OP;
lOR = lperm*OR';


p_old = lperm*po;

ac1 = lOR'*lTransF*lOP;
ac1i = ac1^-1;

[L,U] = ilu(lTransF,struct('type','nofill','droptol',1e-4));

%OR = OP';
while index <= 2
    r_old = lF - (lTransF*p_old);
    
    %dp1 = OP* aci * OR'*r_old;
    
    dp1 = lOP * (ac1i)* lOR*r_old;
    %dp1 = OP*((OR*r_old)'\ac );  
    
    
    rp1 = r_old - lTransF*dp1 ;
   
    
       
  %[p_old,fl1,rr1,it1,rv1]=bicgstab(M_old,RHS_old,1e-10,1000,L,U);
   [dp2,fl1,rr1,it1,rv1]=gmres(lTransF,rp1,10,1e-4,8000,L,U, rp1);
    
   
    %dp2 = TransF\rp1;
    %spparms('spumoni',2)
    p_old = p_old + dp1 + dp2;
    %p_old = p_old +  dp2;

    index = index + 1;
end
% 
r_old = lF - (lTransF*p_old);
%dp1 = OP* ac * OP'*r_old;
dp1 = lOP * ac1i * lOR*r_old;
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
    


p_last = lperm'*p_old;
end

