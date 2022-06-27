
function [Delta,Sum_Delta] = FuncDeltaLM(Jacobian,Error,Lambda)

[nRowNum,nColumNum] = size(Jacobian);
C = sparse(1:nColumNum,1:nColumNum,1)*Lambda;
A = Jacobian'*Jacobian+C;
B = -Jacobian'*Error;
clear Jacobian Error Lambda C;
Delta = A\B;
Sum_Delta = Delta'*Delta;
