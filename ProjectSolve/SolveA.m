function [x,dx] = SolveA(A,s,R,y0)
%UNTITLED2 
%   
%   min||x-s||
%   ||Ax-y0||<R
[Lx,Ly]=size(s);

T=A.*A+ones([Lx,Ly]);
T=1./T;

y = Py(A.*s,R,y0);
u=zeros([Lx,Ly]);
x_old=zeros([Lx,Ly]);
dx=zeros([1,100]);
for i=1:100

x=T.*(s+A.*(y-u));
yv=A.*x+u;
y=Py(yv,R,y0);
r=A.*x-y;
u=u+r;

dx(i)=norm(x-x_old,'fro');
x_old=x;
end



end

