function image = ADMMSolve(x0,filter,P,Maxiter,name)
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
TF=isempty(name);
if(~TF)
path=['data/',name];
status = mkdir(path);
end

x=x0;
y=x0;
u=x-y;

for i=1:Maxiter  
    v1=y-u;
    x=filter(v1);
    x(x0==0)=0;
    v2=x+u;
    y=P(v2);
    u=u+x-y;

 if(~TF)   
    disp([name,' i= ',num2str(i)]);  
    save([path,'/record',num2str(i),'.mat'],'x','y','u');
    imageW=(x-min(x(:)))/(max(x(:))-min(x(:)));
    imwrite(imageW,[path,'/image',num2str(i),'.bmp']);
 end
end

image=x;
end

