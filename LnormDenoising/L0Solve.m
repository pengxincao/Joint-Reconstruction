function D=L0Solve(P,k)
%UNTITLED7 此处显示有关此函数的摘要
%   此处显示详细说明

D=P;
D(abs(P)<k)=0;

end

