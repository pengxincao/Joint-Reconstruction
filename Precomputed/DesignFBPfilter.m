function H = DesignFBPfilter(d,n)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

L=2*n+1;
w=zeros(L,1);
w(2:n+1)=(1:n)/(d*L);
w(n+2:end)=w(n+1:-1:2);

A=1/(2*d);
P=(A.^2)./(A.^2-w.^2);
R=P.*exp(-P);
S=max(R);

H=w.*R/S;
end

