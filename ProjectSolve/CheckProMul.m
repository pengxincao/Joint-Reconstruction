function [TF,D] = CheckProMul(x0,W,varargin)
%UNTITLED6 此处显示有关此函数的摘要
%   此处显示详细说明


n=length(varargin)/2;
O=varargin(1:2:end);
R=varargin(2:2:end);


TF=true(1,n);  %假设均在约束内
D=zeros(1,n);
for i=1:n 
    F=W.*(x0-O{i});
    DD=sum(F.*F,'all');
    if DD>R{i}^2 
        TF(i)=false;
        D(i)=sqrt(DD)-R{i};
    end
end

end

