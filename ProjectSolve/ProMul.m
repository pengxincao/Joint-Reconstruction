function x= ProMul(x0,W,varargin)
% 多约束的投影求解
% varargin O1 R1 O2 R2 ...
% ||W.*(x0-O)||<R

n=length(varargin)/2;
O=varargin(1:2:end);
R=varargin(2:2:end);


TF=true;  %假设均在约束内

for i=1:n 
    F=W.*(x0-O{i});
    DD=sum(F.*F,'all');
    if DD>R{i}^2 
        TF=false;
        break;  %有一个不符合，即重新投影
    end
end

if TF  
    x=x0;
else
    P=cell([1,n]);
    
    for i=1:n 
        P{i}=@(y)Py(y,R{i},W.*O{i});
    end
    
    x=SloveMul(W,x0,P);
end


end