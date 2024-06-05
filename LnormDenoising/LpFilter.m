function x = LpFilter(x0,w,normP,lambda, kappa)
%UNTITLED5 
%   
Mask=x0==0;
betamax = 1e6;
beta = 2*max(lambda);
x=x0;
sizeI2D=size(x0);
if ~iscell(w)
    w={w};
end

n=length(w);
Fw=cell(1,n);
Denormin2=zeros(sizeI2D);

for i=1:n
wconv=rot90(w{i},2);
Fw{i}= psf2otf(wconv,sizeI2D);
Denormin2=Denormin2+abs(Fw{i}).^2;
end

Normin1=fft2(x0);
P=cell(1,n);

while beta<betamax
Denormin1=1+beta*Denormin2;
% Step 1
PL=zeros(sizeI2D);
for i=1:n  
    P{i}=imfilter(x,w{i},'circular');
    PL=PL+P{i}.^2;
end
PL=sqrt(PL);

switch normP
    case 0
        DL=L0Solve(PL,sqrt(lambda/beta));
    case 1
        DL=sign(PL).*max(abs(PL)-lambda/(2*beta),0);
        %disp('p==1');
    otherwise
        DL = LpSolve(PL,lambda/(2*beta),normP);
        %disp('p else');
end
   D=L2vector(DL,PL,P);
%Step 2
Normin2=zeros(sizeI2D);
for i=1:n 
    Normin2= Normin2+conj(Fw{i}).*fft2(D{i});
end

S=(Normin1+beta*Normin2)./(Denormin1);
x=ifft2(S,'symmetric' );
% x(Mask)=0;
% x(x>1)=1;
% x(x<0)=0;
% 更新beta
beta=beta*kappa;
end



end

