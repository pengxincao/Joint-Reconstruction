function [yp] = FBPfilterH(yp,H)
%UNTITLED15 
%   
NumDec=size(yp,1);
L=size(H,1);


yp(L,1)=0;
yp = fft(yp);
yp= bsxfun(@times, yp, H);
yp = ifft(yp,'symmetric');
yp(NumDec+1:end,:)=[];
end