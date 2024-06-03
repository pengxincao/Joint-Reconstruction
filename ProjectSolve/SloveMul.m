function [x,dx] =SloveMul(A,s,P)
[Lx,Ly]=size(s);

T=A.*A+ones([Lx,Ly]);
T=1./T;

%y = Py(A.*s,R,y0);
y=PyMul(A.*s,P);  
u=zeros([Lx,Ly]);
x_old=zeros([Lx,Ly]);
dx=zeros([1,100]);

for i=1:100

x=T.*(s+A.*(y-u));
yv=A.*x+u;

%y=Py(yv,R,y0);
y=PyMul(yv,P);  

r=A.*x-y;
u=u+r;

dx(i)=norm(x-x_old,'fro');
x_old=x;
end


end