function nflag= contflagno(bedge)
global  bcflag coord

benchmark = 'homogeneo';
% determina se o flag do n� interior esta na fronteira de Neumann ou Dirichlet

nflag=50000*ones(size(coord,1),2); % flags 7.15934 � uma constante

for ifacont=1:size(bedge,1)
    a=coord(bedge(ifacont,1),:);
    x=a(1,1);
    y=a(1,2);
    %%
    switch benchmark
        case {'homogeneo', 'heterogeneo','nikitin','durlofsky','lamine',...
                'shuec','altamenteheterogeneo','nikitin1','pinchout'}
            x=bcflag(:,1)==bedge(ifacont,4);
            r=find(x==1);
            nflag(bedge(ifacont,1),2)=bcflag(r,2);
            nflag(bedge(ifacont,1),1)=bcflag(r,1);
        case 'crumpton'
            
            %%
            alfa=1000;
            if x<0 || x==0
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=(2*sin(y)+cos(y))*alfa*x + sin(y)+3*alfa;
            else
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=exp(x)*sin(y)+3*alfa;
            end
        case 'crumptonhyman'
            %%
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=exp(x*y);
        case 'gaowu2'
            %%
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=0.5*((sin((1-x)*(1-y))/(sin(1)))+(1-x)^3*(1-y)^2);
            
        case 'gaowu1'
            %%
            if (x<0.5 || x==0.5)
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=(1+(x-0.5)*(0.1+8*pi*(y-0.5)))*exp(-20*pi*((y-0.5)^2));
                
            else
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=exp(x-0.5)*exp(-20*pi*((y-0.5)^2));
                
            end
            
        case  'gaowu3'
            %%
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=exp(-20*pi*((x-0.5)^2 + (y-0.5)^2));
        case 'gaowu4'
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=sin(pi*x)*sin(pi*y);
        case 'lepotier'
            %%
            
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=sin(pi*x)*sin(pi*y);
            
        case 'lipnikov1'
            %%
            if x<0.5 || x==0.5
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=1-2*y^2+4*x*y+6*x+2*y;
                
            else
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=-2*y^2+1.6*x*y-0.6*x+3.2*y+4.3;
                
            end
        case 'edwards'
            %%
            a1=1;
        f=((4*a1)/(((1/50)-2)*(1/10)+1));
        b2=((1/10)-1)*f;
        c2=f;
        d2=-c2*(1/10);
        c1=(1/50)*(1/10)*f;
        d1=d2;
            if x<0.5
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=c1*x^2+d1*y^2;
                
            else
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=1+b2*x+c2*x^2+d2*y^2;
                
            end
         case 'shenyuan16'
            
             nflag(bedge(ifacont,1),1)=101;
             nflag(bedge(ifacont,1),2)=sin(pi*x)*sin(pi*y);    
        case 'lipnikov2'
            %%
            if (y==0.4444)||(y==0.5555)
                nflag(bedge(ifacont,1),1)=102;
                nflag(bedge(ifacont,1),2)=2;
            elseif (y==0)||(y==1)
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=0;
            end
            if (x==0.4444)||(x==0.5555)
                nflag(bedge(ifacont,1),1)=102;
                nflag(bedge(ifacont,1),2)=2;
            elseif (x==0)||(x==1)
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=0;
            end
        case 'guangwei'
            %%
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=0;
        case 'guangwei1'
            %%
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=0;
            
        case 'gaowu5'
            %%
            if (((0<x || x==0 )&& (x <0.2 || x==0.2)) && y==0) || (((0<y || y==0 )&& (y <0.2 || y==0.2)) && x==0)
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=1;
            elseif (((0.8<x || x==0.8 )&& (x <1 || x==1)) && y==1) || (((0.8<y || y==0.8 )&& (y <1 || y==1)) && x==1)
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=0;
            elseif (((0.3<x || x==0.3 )&& (x <1 || x==1)) && y==0) || (((0.3<y || y==0.3 )&& (y <1 || y==1)) && x==0)
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=0.5;
            elseif (((0<x || x==0 )&& (x <0.7 || x==0.7)) && y==1) || (((0<y || y==0 )&& (y <0.7 || y==0.7)) && x==1)
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=0.5;
            else
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=0.5;
            end
            
        case 'gaowu6'
            %%
            nflag(bedge(ifacont,1),1)=201;
            nflag(bedge(ifacont,1),2)=0;
            
        case 'gaowu7'
            %%
            delta=0.2;
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=-x-delta*y;
        case 'gaowu8'
            %%
            delta=0.2;
            phi1=y-delta*(x-0.5)-0.475;
            phi2=phi1-0.05;
            % dominio 1
            if phi1<0 || phi1==0
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=-phi1;
                %dominio
            elseif phi1>0 && phi2<0
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=-100*phi1;
            elseif phi2>0 || phi2==0
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=-phi2-5;
            end
            
            
        case 'gaowu9'
            %%
            delta=0.2;
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=2-x-delta*y;
            
        case {'benchmar5_7','benchmar5_6'}
            %%
            
            nflag(bedge(ifacont,1),1)=101;
            nflag(bedge(ifacont,1),2)=0;
            
        case 'edqueiroz'
            %%
            if bedge(ifacont,4)==101
                nflag(bedge(ifacont,1),1)=101;
                nflag(bedge(ifacont,1),2)=0;
            else
                nflag(bedge(ifacont,1),1)=102;
                nflag(bedge(ifacont,1),2)=2;
            end
            
    end
end