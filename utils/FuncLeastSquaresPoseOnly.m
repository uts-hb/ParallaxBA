
function [PVector,Reason,Info] = FuncLeastSquaresPoseOnly(xVector,PVector,Feature,K)

%Initial Error
[nRowNum,nColumnNum] = size(xVector.PID);
[Error, Sum_Error]= FuncDiff(xVector,PVector,Feature,K);
Sum_Error = Sum_Error/nRowNum;
fprintf('Initial Error is %.8f\n', Sum_Error);
Sum_Delta = 22;
MaxIter = 21;
MinError = 0.000001;
MinDelta = 0.000000001;
% MinDelta = 0;

Iter = 0;
while Sum_Error>MinError && Sum_Delta>MinDelta && Iter<=MaxIter;
%     [Jacobian] = FuncJacobian(xVector,PVector,Feature,K);
    [Jacobian] = FuncJacobianPoseOnly(xVector,PVector,Feature,K);
    [Delta,Sum_Delta] = FuncDelta(Jacobian,Error);
%     [Delta,Sum_Delta] = FuncDeltaCG(Jacobian,Error,Iter);
    [PVector] = FuncUpdatePoseOnly(PVector,Delta);
    [Error, Sum_Error]= FuncDiff(xVector,PVector,Feature,K);
    Sum_Error = Sum_Error/nRowNum;
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