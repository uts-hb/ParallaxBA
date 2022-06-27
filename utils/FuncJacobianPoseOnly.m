
function [Jacobian] = FuncJacobianPoseOnly(xVector,PVector,Feature,K)

[nRowNumu,nColumNum] = size(xVector.u);
[nRowNumP,nColumNum] = size(PVector.Pose);
[nRowNumF,nColumNum] = size(PVector.Feature);
[nRowNumO,nColumNum] = size(PVector.Omega);
Jacobian1 = sparse(nRowNumu,nRowNumP+nRowNumF+nRowNumO);

for j=1:nRowNumu/2;
    PID = xVector.PID(j,1);
    FID = xVector.FID(j,1);
    if PID == Feature(FID,3);
        Pi = PVector.Pose(6*(PID-1)+1:6*PID)';
        Xj = PVector.Feature(2*FID-1:2*FID,1);
        [dudRi,dudPT] = FuncdudPi(K,Pi,Xj);
        Jacobian1(2*j-1:2*j,6*(PID-1)+1:6*(PID-1)+3) = dudRi;
        Jacobian1(2*j-1:2*j,nRowNumP+2*FID-1:nRowNumP+2*FID) = dudPT;
        clear Pi Xj dudRi dudPT;
    elseif PID == Feature(FID,4);
        Pi = PVector.Pose(6*(Feature(FID,3)-1)+1:6*Feature(FID,3))';
        Pk = PVector.Pose(6*(PID-1)+1:6*PID)';
        Xj = cat(1,PVector.Feature(2*FID-1:2*FID,1),PVector.Omega(FID,1));
        [dudRk,dudTi,dudTk,dudXj] = FuncdudPk(K,Pi,Pk,Xj);
        Jacobian1(2*j-1:2*j,6*(PID-1)+1:6*(PID-1)+3) = dudRk;
        Jacobian1(2*j-1:2*j,6*(Feature(FID,3)-1)+4:6*Feature(FID,3)) = dudTi;
        Jacobian1(2*j-1:2*j,6*(PID-1)+4:6*PID) = dudTk;
        Jacobian1(2*j-1:2*j,nRowNumP+2*FID-1:nRowNumP+2*FID) = dudXj(:,1:2);
        Jacobian1(2*j-1:2*j,nRowNumP+nRowNumF+FID) = dudXj(:,3);
        clear Pi Pk Xj dudRk dudTi dudTk dudPT;
    else
        Pi = PVector.Pose(6*(Feature(FID,3)-1)+1:6*Feature(FID,3))';
        Pk = PVector.Pose(6*(Feature(FID,4)-1)+1:6*Feature(FID,4))';
        Pl = PVector.Pose(6*(PID-1)+1:6*PID)';
        Xj = cat(1,PVector.Feature(2*FID-1:2*FID,1),PVector.Omega(FID,1));
        [dudRl,dudTi,dudTk,dudTl,dudXj] = FuncdudPl(K,Pi,Pk,Pl,Xj);
        Jacobian1(2*j-1:2*j,6*(PID-1)+1:6*(PID-1)+3) = dudRl;
        Jacobian1(2*j-1:2*j,6*(Feature(FID,3)-1)+4:6*Feature(FID,3)) = dudTi;
        Jacobian1(2*j-1:2*j,6*(Feature(FID,4)-1)+4:6*Feature(FID,4)) = dudTk;
        Jacobian1(2*j-1:2*j,6*(PID-1)+4:6*PID) = dudTl;
        Jacobian1(2*j-1:2*j,nRowNumP+2*FID-1:nRowNumP+2*FID) = dudXj(:,1:2);
        Jacobian1(2*j-1:2*j,nRowNumP+nRowNumF+FID) = dudXj(:,3);
        clear Pi Pk Pl Xj dudRl dudTi dudTk dudTl dudXj;
    end;
end;

% Jacobian = cat(2,Jacobian1(1:end,7:9),Jacobian1(1:end,11:end));
Jacobian = cat(2,Jacobian1(1:end,7:11),Jacobian1(1:end,13:nRowNumP));%Z
% Jacobian = cat(2,Jacobian1(1:end,7:10),Jacobian1(1:end,12:end));%Y
% Jacobian = Jacobian1(1:end,7:end);