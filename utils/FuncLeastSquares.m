
function [PVector,Reason,Info,Sum_Error, Errors] = FuncLeastSquares(xVector,PVector,Feature,K,FixVa)
Errors = {}; 
%Initial Error
nRowNum = length(xVector.PID);
[Error, Sum_Error,reprojectionErrors,Errors]= FuncDiff(xVector,PVector,Feature,K,Errors);
Sum_Error = Sum_Error/nRowNum;
fprintf('Initial Error is %.8f\n', Sum_Error);
Sum_Delta = 22;
MaxIter = 30;
MinError = 1e-8;
MinDelta = 1e-10;

[ID1,ID2,nJacobian] = FuncGetJacobianID(xVector,PVector,Feature);

Iter = 0;
while Sum_Error>MinError && Sum_Delta>MinDelta && Iter<=MaxIter;
    [Jacobian] = FuncJacobian(xVector,PVector,Feature,K,ID1,ID2,nJacobian);
    [DeltaP,DeltaF,Sum_Delta] = FuncDelta(Jacobian,Error,PVector,Feature(:,2),FixVa);
    [PVector] = FuncUpdate(PVector,DeltaP,DeltaF);
[Error, Sum_Error,reprojectionErrors,Errors]= FuncDiff(xVector,PVector,Feature,K,Errors);
    Sum_Error = Sum_Error/nRowNum;
% Sum_Error = Sum_Error;
    Iter = Iter+1;
    fprintf('Iterations %d Error %.8f\n', Iter,Sum_Error);
    
end;

if Sum_Error<MinError;
    Reason = 1;
elseif Sum_Delta<MinDelta;
    Reason = 2;
elseif Iter>MaxIter;
    Reason = 3;
else
    Reason = 4;
end;

Info = sparse([]);
if Iter>0;
    Info = Jacobian'*Jacobian;
end;