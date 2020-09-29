function [F,V,N]=elementface2(nflag)
global bedge elem coord esurn1 esurn2 nsurn1 nsurn2 inedge Nregion Vregion npar
%

si=size(inedge,1);
F=zeros(size(elem,1),6);

ny=0;
for y=1:size(esurn2,1)-1
    n=esurn2(y+1)-esurn2(y);
    if n>ny
        ny=n;
    end
end

N=zeros(size(coord,1),ny);
V=zeros(2,ny,size(coord,1));

for i=1:size(elem,1)
   t=1;
   for j=1:size(inedge,1)
       if (inedge(j,3)==i)||(inedge(j,4)==i)
           F(i,t)=j;
           t=t+1;
       end
   end
   for j=1:size(bedge,1)
      if bedge(j,3)==i
          F(i,t)=j+si;
          t=t+1;
      end
   end
end 

for y=1:size(coord,1)
    if (isempty(find(bedge(:,1)==y))==0)||(isempty(find(bedge(:,2)==y))==0)
        isbedge=1; % se for de contorno.
    else
        isbedge=0; % se for interno.
    end
    ve=esurn1(esurn2(y)+1:esurn2(y+1));
    if isbedge==0
        for g=1:size(ve,1)-1
            for k=1:size(F,2)
                a=F(ve(g),k);
                f1=find(F(ve(g+1),:)==a);
                f=F(ve(g+1),f1);
                if isempty(f1)==0
                    break
                end
            end
            v(g)=f;
        end
        for k=1:size(F,2)
            a=F(ve(size(ve,1)),k);
            f1=find(F(ve(1),:)==a);
            f=F(ve(1),f1);
            if isempty(f1)==0
                break
            end
        end
        v(size(ve,1))=f;
    else
        for g=1:size(ve,1)-1
            for k=1:size(F,2)
                a=F(ve(g),k);
                f1=find(F(ve(g+1),:)==a);
                f=F(ve(g+1),f1);
                if isempty(f1)==0
                    break
                end
            end
            v(g+1)=f;
        end
        f1=find(F(ve(1),:)>size(inedge,1));
        if size(f1,2)>1
            eb=ve(1);
            pe=find(elem(eb,:)==y);
            cf1=elem(eb,pe:4);
            cf2=elem(eb,1:pe-1);
            cf=horzcat(cf1,cf2);
            bf1=bedge(F(ve(1),f1(1))-size(inedge,1),1:2);
            bf2=bedge(F(ve(1),f1(2))-size(inedge,1),1:2);
            if bf1(1)==y
                bff1=bf1(2);
            else
                bff1=bf1(1);
            end
            if bf2(1)==y
                bff2=bf2(2);
            else
                bff2=bf2(1);
            end
            ff1=intersect(cf,bff1);
            ff2=intersect(cf,bff2);
            if isempty(ff1)==0
                f2(1)=find(cf==ff1);
            else
                f2(1)=0;
            end
            if isempty(ff2)==0
                f2(2)=find(cf==ff2);
            else
                f2(2)=0;
            end
            if f2(1)>f2(2)
                f1=f1(2);
            else
                f1=f1(1);
            end 
        end
        v(1)=F(ve(1),f1);
        f1=find(F(ve(size(ve,1)),:)>size(inedge,1));
        if size(f1,2)>1
            eb=ve(size(ve,1));
            pe=find(elem(eb,:)==y);
            cf1=elem(eb,pe:4);
            cf2=elem(eb,1:pe-1);
            cf=horzcat(cf1,cf2);
            bf1=bedge(F(ve(size(ve,1)),f1(1))-size(inedge,1),1:2);
            bf2=bedge(F(ve(size(ve,1)),f1(2))-size(inedge,1),1:2);
            if bf1(1)==y
                bff1=bf1(2);
            else
                bff1=bf1(1);
            end
            if bf2(1)==y
                bff2=bf2(2);
            else
                bff2=bf2(1);
            end
            ff1=intersect(cf,bff1);
            ff2=intersect(cf,bff2);
            if isempty(ff1)==0
                f2(1)=find(cf==ff1);
            else
                f2(1)=0;
            end
            if isempty(ff2)==0
                f2(2)=find(cf==ff2);
            else
                f2(2)=0;
            end
            if f2(1)<f2(2)
                f1=f1(2);
            else
                f1=f1(1);
            end 
        end
        v(size(ve,1)+1)=F(ve(size(ve,1)),f1);
    end
    n=size(v,2);
    N(y,1:n)=v;
    if isbedge==1 % Se for de contorno
        V(1,1:n-1,y)=v(1:n-1);
        V(2,1:n-1,y)=v(2:n);
    else
        V(1,1:n,y)=v;
        V(1,n,y)=v(1);
        V(2,1:n-1,y)=v(2:n);
        V(2,n,y)=v(n);
    end   
    clear v;
end



Nregion = zeros(size(N,1),size(N,2),npar);
%tmp = zeros(size(N));

for ii = 1: npar
    tmp = zeros(size(N));
    tmp = N;
   
    if ii == 51
        1+1;
    end
    tmp(~filtEdges2(N,ii)) = 0; 
 
   
    for jj = 1:size(tmp,1)
        tmp(jj,:) = zeroswift(tmp(jj,:));   % ordenação sentido anti-horário     
    end
    %% 
    Nregion(:,:,ii) = tmp;
end


Vregion = zeros(size(V,1),size(V,2),size(V,3),npar);



for region = 1:npar
   
    for ii = 1:size(coord,1)
        auxmat = V(:,:,ii);
        auxvec = filtEdges2(V(:,:,ii),region);
        auxvec = auxvec(1,:) .*auxvec(2,:);
        auxvec = [auxvec ; auxvec];
        
        Vregion(:,:,ii,region) = zeroSwiftMat( auxmat.*auxvec);
        
    end
    
    
end



end



function out = filtEdges2(matEd, region)
    global inedge bedge elemloc
    %out = zeros(size(matEd));
    auxIn1 = zeros(size(matEd));
    auxIn2 = zeros(size(matEd));

    auxBed = zeros(size(matEd));
    
    refBedge= matEd > size(inedge,1);
    refInedge = (matEd <= size(inedge,1)) & matEd ~= 0;
    
    
    
    auxIn1(refInedge) = elemloc( inedge( matEd(refInedge),3));
    auxIn2(refInedge) = elemloc( inedge( matEd(refInedge),4));
    auxBed(refBedge )= elemloc(bedge(matEd(refBedge)- size(inedge,1),3));
    
    
    out = (auxIn1 == region) | (auxIn2 == region) | (auxBed== region) ;
end

function out = zeroswift(vec)
    b = find(vec==0);
    a = find(vec~=0);
    if (isempty(a)==0)&&(isempty(b)==0)
        c = a>b(1);
        if any(c)~=0
            ref = min(find(c==1));
            prim = vec(a(ref));
            while vec(1)~=prim
               vec = circshift(vec,1); 
            end
        end
        pos = find(vec~=0);
        vec(1:size(pos,2)) = vec(pos);
        vec(size(pos,2)+1:end) = 0;
    end
    out = vec;
    %out = circshift(vec,-1);
end

function out = antihorarioOrder(vec, node)
    global bedge inedge
    ldist = @(x,y) (x.^2 + y.^2)^(0.5);
    ref = size(inedge,1);
    bedge_ref = vec > ref;
    inedge_ref = (~bedge_ref) & (vec ~= 0); 
    
    all_nodes_b = reshape(bedge(vec(bedge_ref) - ref,1:2),1,[]);
    all_nodes_in = reshape(inedge(vec(inedge_ref),1:2),1,[]);
    node = mode([all_nodes_b, all_nodes_in]);
    
    if all(vec ~= 0)
       1+1 
    end
    
        
    end


function mat = zeroSwiftMat(mat)
    ref = mat(1,:) ~= 0;
    auxmat = zeros(2, size(find(ref),2));
    auxmat(1,:) = mat(1,ref);
    auxmat(2,:) = mat(2,ref);    
    mat = zeros(size(mat));
    mat(1,1:size(auxmat,2)) = auxmat(1,:);
    mat(2,1:size(auxmat,2)) = auxmat(2,:);
    
end

function vec = zeroSwift(vecInp,node,region)
    sz = size(vecInp,2);
    vec = zeros(1,sz);    
    p = vecInp(vecInp ~=0);  
%     tmp = zeros(size(p,1),4);
%     if no == 2
%         1+1;
%     end
%     
    
%     refP = refConv(p);
    
if ~isempty(p)
 
    p = ordEdges(node,p',region);
    
   % tmp(refP(:,1) == 1,1) = bedge(
    
    
    vec(1: 1+ size(p,1)-1) = p;

    
end
end

