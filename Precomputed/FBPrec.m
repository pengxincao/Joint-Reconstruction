function image = FBPrec(y,A)
%UNTITLED8 此处显示有关此函数的摘要
%   此处显示详细说明
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%滤波器设计

order = max(64,2^nextpow2(2*size(y,1)));
n=order/2;
H = DesignFBPfilter(A.d_dec,n);

[yp] = FBPfilterH(y,H);

image=FBPBackProjectC(yp,A.thetapi,A.d_pixel,A.N,A.d_dec,A.Lpos,A.subpixel);
image=pi*image/(length(A.thetapi));
%image=pi*x;
%image=line2map(image);
%image=nomlize(image,[],[]);
end


