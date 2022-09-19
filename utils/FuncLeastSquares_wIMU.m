
function [PVector,Reason,Info,Sum_Error, Errors, reprojectionErrors_initial] = FuncLeastSquares_wIMU(xVector,PVector,Feature,K, camera, measurements)
Errors = {}; 
%Initial Error
nRowNum = length(xVector.PID);
[Error, Sum_Error,reprojectionErrors,Errors]= FuncDiff(xVector,PVector,Feature,K,Errors);
%% ============================================================= %%
% Error State for IMU

measurements = measurements(~cellfun('isempty',measurements));
[imu_Error, imu_Sum_Error]= FuncDiffvw_v2(PVector.Pose, camera, measurements);
Error = [Error; imu_Error]; 
Sum_Error = Error'*Error;
Sum_Error = Sum_Error/(nRowNum + (size(imu_Error,1)/6));

%% ============================================================= %%

Sum_Error = Sum_Error/nRowNum;
fprintf('Initial Error is %.8f\n', Sum_Error);
reprojectionErrors_initial = Sum_Error;
Sum_Delta = 22;
MaxIter = 20;
MinError = 1e-8;
MinDelta = 1e-10;

[ID1,ID2,nJacobian] = FuncGetJacobianID_wIMU(xVector,PVector,Feature);

Iter = 0;
while Sum_Error>MinError && Sum_Delta>MinDelta && Iter<=MaxIter;
    [Jacobian] = FuncJacobian_wIMU(xVector,PVector,Feature,K,ID1,ID2,nJacobian,camera,measurements);
%% ============================================================= %%
        w_ID1 = [1:size(xVector.u,1)];
        w_ID2 = [1:size(xVector.u,1)];
        w_val = ones(1,size(xVector.u,1));

        v_Cov = 4e-2;
        w_Cov = 4e-2;
        weight = 1000; 

        w_LID1 = [size(xVector.u,1)+1:size(xVector.u,1)+size(imu_Error,1)];
        w_LID2 = [size(xVector.u,1)+1:size(xVector.u,1)+size(imu_Error,1)];
        w_ID1 = [w_ID1,w_LID1];
        w_ID2 = [w_ID2,w_LID2];
        w_Val = [w_val, weight*inv(v_Cov)*ones(1,size(imu_Error,1))];

        W = sparse(w_ID1,w_ID2,w_Val);
%% ============================================================= %%
    [DeltaP,DeltaF,Sum_Delta] = FuncDelta_wIMU(Jacobian,Error,PVector,Feature(:,2),W);
    [PVector] = FuncUpdate(PVector,DeltaP,DeltaF);
    [Error, Sum_Error, reprojectionErrors, Errors]= FuncDiff(xVector,PVector,Feature,K,Errors);
    

    [imu_Error, imu_Sum_Error]= FuncDiffvw_v2(PVector.Pose, camera, measurements);

    Error = [Error; imu_Error];
    Sum_Error = Error'*Error;
    Sum_Error = Sum_Error/(nRowNum + (size(imu_Error,1)/6));

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
    Info = Jacobian'*W*Jacobian;
end;