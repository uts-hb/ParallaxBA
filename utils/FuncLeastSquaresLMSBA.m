
function [PVector,Reason,Info] = FuncLeastSquaresLMSBA(xVector,PVector,Feature,K,FixVa)

nRowNum = length(xVector.PID);
MaxIter = 50;
Factor = 2;
t = 1e-3;
e1 = 1e-12;
e2 = 1e-12;
e3 = 1e-12;
e4 = 0;
Stop = 0;
Reason = 0;
Iter = 0;

[ID1,ID2,nJacobian] = FuncGetJacobianID(xVector,PVector,Feature);

[Error, Sum_Error]= FuncDiffSBA(xVector,PVector,Feature,K);
Sum_Error2 = Sum_Error/nRowNum;
fprintf('Initial Error is %.8f\n', Sum_Error2);
ErrorPre = Sum_Error;
[Jacobian] = FuncJacobian(xVector,PVector,Feature,K,ID1,ID2,nJacobian);
A = Jacobian'*Jacobian;
Lambda = t*max(diag(A));
G = Jacobian'*Error;
g = max(abs(G));
if g<=e1;
    Stop = 1;
    Reason = 1;
end;

while Stop~=1 && Iter<=MaxIter;
    Iter = Iter+1;
    P = -1;
    while Stop~=1 && P<=0;
    [DeltaP,DeltaF,Sum_Delta] = FuncDeltaLMSBA(A,G,Lambda,PVector,Feature(:,2),FixVa);
    P2 = FuncGetP2(PVector,FixVa);
    if Sum_Delta<=e2*(P2+e2);
        Stop = 1;
        Reason = 2;
    else
        [PVector] = FuncUpdate(PVector,DeltaP,DeltaF);
        [Error, Sum_Error]= FuncDiffSBA(xVector,PVector,Feature,K);
        Delta = [DeltaP;DeltaF];
        P = (ErrorPre-Sum_Error)/(Delta'*(Lambda*Delta+G));
        if P>0;
            if sqrt(ErrorPre)-sqrt(Error)<e4*sqrt(ErrorPre);
                Stop = 1;
                Reason = 3;
            end;
            [Jacobian] = FuncJacobian(xVector,PVector,Feature,K,ID1,ID2,nJacobian);
            A = Jacobian'*Jacobian;
            G = Jacobian'*Error;
            g = max(abs(G));
            if Stop ==1 || g<=e1;
                Stop = 1;
                Reason = 1;
            end;
            Lambda = Lambda*max(1/3,1-(2*P-1)^3);
            Factor = 2;
        else    
            DeltaP = -DeltaP;
            DeltaF = -DeltaF;
            [PVector] = FuncUpdate(PVector,DeltaP,DeltaF);
            Lambda = Factor*Lambda;
            Factor = Factor*2;
        end;
    end;
    end;
    if sqrt(Sum_Error)<=e3;
        Stop = 1;
        Reason = 4;
    end;
    if P>0;
        Sum_Error2 = Sum_Error/nRowNum;
        fprintf('Iterations %d Error %.8f\n', Iter, Sum_Error2);
        ErrorPre = Sum_Error;
    end;
end;
    
Info = sparse([]);
if Iter>0;
    Info = Jacobian'*Jacobian;
end;