
function [PVector,Reason,Info] = FuncLeastSquaresLM(xVector,PVector,Feature,K)

%Initial Error
[nRowNum,nColumnNum] = size(xVector.PID);
[Error, Sum_Error]= FuncDiff(xVector,PVector,Feature,K);
Sum_Error = Sum_Error/nRowNum;
fprintf('Initial Error is %.8f\n', Sum_Error);
Sum_Delta = 22;
ErrorPre = Sum_Error;
MaxIter = 222;
MinError = 0.001;
MinDelta = 0.0000000000000000000000000000001;
Factor = 10;
% MinDelta = 0;

Iter = 0;
while Sum_Error>MinError && Sum_Delta>MinDelta && Iter<=MaxIter;
    [Jacobian] = FuncJacobian(xVector,PVector,Feature,K);
    if Iter == 0;
%         Lambda = abs(mean(diag(Jacobian'*Jacobian)))*0.001;
%         Lambda = mean(abs(diag(Jacobian'*Jacobian)))*0.001;
        Lambda = mean(diag(Jacobian'*Jacobian))*0.001;
    end;
%     [Delta,Sum_Delta] = FuncDeltaLM(Jacobian,Error,Lambda);
    [Delta,Sum_Delta] = FuncDeltaLM(Jacobian,Error,Lambda);
    [PVector] = FuncUpdate(PVector,Delta);
    [Error, Sum_Error]= FuncDiff(xVector,PVector,Feature,K);
    Sum_Error = Sum_Error/nRowNum;
    Iter = Iter+1;
    fprintf('Iterations %d Error %.8f\n', Iter,Sum_Error);
    if Sum_Error>ErrorPre;
        Delta = -Delta;
        [PVector] = FuncUpdate(PVector,Delta);
        Lambda = Factor*Lambda;
        Sum_Delta = 22;
    else
        Lambda = Lambda/Factor;
        ErrorPre = Sum_Error;
    end;
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