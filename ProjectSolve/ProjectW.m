function x= ProjectW(x0,O,W,R)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
%   (W.*（x-o）).^2<R^2
%   将x0向此集合投影
b=x0-O;

F=W.*b;
DD=sum(F.*F,'all');

if DD<=R*R
    x=x0;
else
    
 [x,~] = SolveA(W,x0,R,W.*O);
    
    
end



end

