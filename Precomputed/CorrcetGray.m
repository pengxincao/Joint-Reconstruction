function image = CorrcetGray(yp,A,SARTdata,Cx)

[Lx,Ly]=size(yp);

Ix1=Cx;
Ix2=Cx.^2;
Ix3=Cx.^3;
Ix4=Cx.^4;
Ix5=sqrt(Cx);
temp=Ix5;
temp(temp==0)=1;
%Ix6=1./temp;

F1=SARTForward(Ix1,A.thetapi,A.d_pixel,A.d_dec,A.Lpos,SARTdata{1},SARTdata{2},SARTdata{3},SARTdata{4},SARTdata{5});
F2=SARTForward(Ix2,A.thetapi,A.d_pixel,A.d_dec,A.Lpos,SARTdata{1},SARTdata{2},SARTdata{3},SARTdata{4},SARTdata{5});
F3=SARTForward(Ix3,A.thetapi,A.d_pixel,A.d_dec,A.Lpos,SARTdata{1},SARTdata{2},SARTdata{3},SARTdata{4},SARTdata{5});
F4=SARTForward(Ix4,A.thetapi,A.d_pixel,A.d_dec,A.Lpos,SARTdata{1},SARTdata{2},SARTdata{3},SARTdata{4},SARTdata{5});
F5=SARTForward(Ix5,A.thetapi,A.d_pixel,A.d_dec,A.Lpos,SARTdata{1},SARTdata{2},SARTdata{3},SARTdata{4},SARTdata{5});
%F6=SARTForward(Ix6,A.thetapi,A.d_pixel,A.d_dec,A.Lpos,SARTdata{1},SARTdata{2},SARTdata{3},SARTdata{4},SARTdata{5});


Y=reshape(yp,[Lx*Ly,1]);
X1=reshape(F1,[Lx*Ly,1]);
X2=reshape(F2,[Lx*Ly,1]);
X3=reshape(F3,[Lx*Ly,1]);
X4=reshape(F4,[Lx*Ly,1]);
X5=reshape(F5,[Lx*Ly,1]);
%X6=reshape(F6,[Lx*Ly,1]);
%X=[X1,X2,X3,X4,X5,X6];
X=[X1,X2,X3,X4,X5];
b = regress(Y,X);

%image=b(1)*Ix1+b(2)*Ix2+b(3)*Ix3+b(4)*Ix4+b(5)*Ix5+b(6)*Ix6;
image=b(1)*Ix1+b(2)*Ix2+b(3)*Ix3+b(4)*Ix4+b(5)*Ix5;
mask=Cx==0;
image(mask)=0;

end