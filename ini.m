%%% Low contrast phantom reconstruction experiment %%%
addpath Precomputed\;
addpath BasicTools\;
addpath LnormDenoising\;

load('SynFdata.mat');
load('Matrix.mat', 'A');
load('Matrix.mat', 'SARTdata');
load('t.mat', 'Cx')

I = FBPrec(SynFdata,A);
imshow(I,[0.00002,0.0001])
Iini=ImageNorm(I,0.00002,0.0001);  

Syn=LpFilter(Iini,{[-1,2,-1],[-1;2;-1]},0,0.003,1.5);

S=zeros(3304,1);
S(1003:1003+99)=((1003:1003+99)-1002)/100;
S(2203:2302)=1-((2203:2302)-2202)/100;
S(1102:2202)=1;

CxC = CorrcetGray(SynFdata,A,SARTdata,Cx);  % Correct the grayscale

F0 = SART_F(A,SARTdata,CxC);  % Forward projection
F=F0.*(1-S)+SynFdata.*S;
Iini0 = FBPrec(F,A);
Iini=ImageNorm(Iini0,0.00002,0.0001);


N01=750;
EF=exp(F)/N01+exp(2*F)*8/(N01^2);
L=3304;Lx=1300;
maskT=L/2-Lx/2+1:L/2+Lx/2;
F2=F(maskT,:);
N02=2000;
EF2=exp(F2)/N02+exp(2*F2)*3/(N02^2);
EF(maskT,:)=EF2;

H = DesignFBPfilter(11,4096);
FBPHt=ifft(H,'symmetric');
FBPH=[FBPHt(end-4095:end);FBPHt(1:end-4096)];  % Calculate hFBP
VarH=FBPH.*FBPH;
EFH = imfilter(EF,VarH,'symmetric');

VarX=FBPBackProjectC(EFH,A.thetapi,A.d_pixel,A.N,A.d_dec,A.Lpos,A.subpixel);
sigmaX=sqrt(VarX);
sigmaX=sigmaX*pi/(2500*0.00008);
W=1./sigmaX;
W(SARTdata{1}==0)=0;
W=W*20;  % Calculate the weight matrix W