function F = SART_F(A,SARTdata,image)
%UNTITLED5 
%   
F=SARTForward(image,A.thetapi,A.d_pixel,A.d_dec,A.Lpos,SARTdata{1},SARTdata{2},SARTdata{3},SARTdata{4},SARTdata{5});
end