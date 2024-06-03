function z=PyMul(x0,P)  

z=x0;
[Lx,Ly]=size(x0);
n=length(P);
u=cell([1,n]);
x=cell([1,n]);
for i=1:n 
    u{i}=zeros([Lx,Ly]);
end


z_old=z;
r=1;

while r>0.001
    xbar=x0;
    for i=1:n  
        x{i}=P{i}(z-u{i}); 
        xbar=xbar+x{i}+u{i};
    end
    z=xbar/(n+1);
    
   for i=1:n  
        u{i}=u{i}+x{i}-z;
   end
   r=norm(z-z_old,'fro')/norm(z_old,'fro');
   z_old=z;
   
end

end