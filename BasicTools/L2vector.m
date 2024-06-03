function D=L2vector(DL,PL,P) 
n=length(P);
D=P;
PL(PL==0)=1;
K=DL./PL;

for i=1:n 
    D{i}=K.*P{i}; 
end


end