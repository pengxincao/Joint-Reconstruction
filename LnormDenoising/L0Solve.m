function D=L0Solve(P,k)
%UNTITLED7 
%   

D=P;
D(abs(P)<k)=0;

end

