function yo = Py(y,R,y0)
%UNTITLED3 
%   
v=y-y0;
d=norm(v,'fro');
if d<=R
    yo=y;
else
    yo=y0+v*R/d;
end


end

