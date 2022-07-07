
function [PVector,Reason,Info] = FuncLeastSquaresLM2(xVector,PVector,Feature,K)

[nRowNum,nColumnNum] = size(xVector.PID);
MaxIter = 5000;
Factor = 2;
T1 = 0.01;
T2 = 0.00000000000001;
T3 = 0.001;
T4 = 0.000001;
Stop = 0;
Reason = 0;
Iter = 0;

[Error, Sum_Error]= FuncDiff(xVector,PVector,Feature,K);
Sum_Error = Sum_Error/nRowNum;
fprintf('Initial Error is %.8f\n', Sum_Error);
ErrorPre = Sum_Error;
[Jacobian] = FuncJacobian(xVector,PVector,Feature,K);
Lambda = T3*max(diag(Jacobian'*Jacobian));
G = sqrt((Jacobian'*Error)'*(Jacobian'*Error));
G2 = sqrt(G'*G);
if G2<=T1;
    Stop = 1;
    Reason = 1;
end;

while Stop~=1 && Iter<=MaxIter;
    [Delta,Sum_Delta,G,B] = FuncDeltaLM2(Jacobian,Error,Lambda);
    P2 = FuncGetP2(PVector);
    if Sum_Delta<=T2*P2;
        Stop = 1;
        Reason = 2;
    else
        [PVector] = FuncUpdate(PVector,Delta);
        [Error, Sum_Error]= FuncDiff(xVector,PVector,Feature,K);
        Sum_Error = Sum_Error/nRowNum;
        Iter = Iter+1;
        fprintf('Iterations %d Error %.8f\n', Iter,Sum_Error);
        P = (ErrorPre-Sum_Error)/(Delta'*(Lambda*Delta+B));
        if P>0;
            if G2<=T1;
                Stop = 1;
                Reason = 1;
            elseif Sum_Error<=T4;
                Stop = 1;
                Reason = 3;
            else
                Lambda = Lambda*max(1/3,1-(2*P-1)^3);
%                 A = max(1/3,1-(2*P-1)^3)
                [Jacobian] = FuncJacobian(xVector,PVector,Feature,K);
                Factor = 2;
                ErrorPre = Sum_Error;
            end;
        else    
            Delta = -Delta;
            [PVector] = FuncUpdate(PVector,Delta);
            Lambda = Factor*Lambda;
            Factor = Factor*2;
        end;
    end;
end;
    
Info = sparse([]);
if Iter>0;
    Info = Jacobian'*Jacobian;
end;