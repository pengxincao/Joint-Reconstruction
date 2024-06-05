function T = LpSolve(Input,lambda,p)
%UNTITLED9 
%   
S=sign(Input);
Y=abs(Input);
T=Y;
t1=(2*lambda*(1-p))^(1/(2-p))+lambda*p*(2*lambda*(1-p))^((p-1)/(2-p));
v=Y<t1;
T(v)=0;

for i=1:30  
    T(~v)=Y(~v)-lambda*p*(T(~v).^(p-1));
end
T=S.*T;





end

