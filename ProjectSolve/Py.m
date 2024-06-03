function yo = Py(y,R,y0)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
v=y-y0;
d=norm(v,'fro');
if d<=R
    yo=y;
else
    yo=y0+v*R/d;
end


end

