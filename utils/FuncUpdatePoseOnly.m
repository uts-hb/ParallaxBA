
function [PVector] = FuncUpdatePoseOnly(PVector,Delta)

% [nRowNumP,nColumNumP] = size(Pose);

% PoseDelta = reshape(Delta(1:6*(nRowNumP-1)),nRowNumP-1,6);
% FeatureDelta = reshape(Delta(6*(nRowNumP-1)+1:end),[],3);

% Delta = cat(1,Delta(1:5),0,Delta(6:end));
% PoseDelta = reshape(Delta(1:6*(nRowNumP-1)),6,nRowNumP-1)';
% FeatureDelta = reshape(Delta(6*(nRowNumP-1)+1:end),3,[])';

% Delta=Delta
% PoseDelta = PoseDelta
% FeatureDelta = FeatureDelta
% pause

[nRowNumP,] = size(PVector.Pose);
[nRowNumF,] = size(PVector.Feature);
[nRowNumO,] = size(PVector.Omega);

PoseDelta = cat(1,Delta(1:5),0,Delta(6:end));%Z
% PoseDelta = cat(1,Delta(1:5),0,Delta(6:nRowNumP-7));%Z
% PoseDelta = cat(1,Delta(1:4),0,Delta(5:nRowNumP-7));%Y
% PoseDelta = cat(1,Delta(1:3),0,Delta(4:nRowNumP-7));
% FeatureDelta = Delta(nRowNumP-6:nRowNumP+nRowNumF-7);
% OmegaDelta = Delta(nRowNumP+nRowNumF-6:end);

PVector.Pose(7:end,1) = PVector.Pose(7:end,1)+PoseDelta;
% PVector.Feature(:,1) = PVector.Feature(:,1)+FeatureDelta;
% PVector.Omega(:,1) = PVector.Omega(:,1)+OmegaDelta;

% if Feature(418,3)<0;
%     Feature(418,3)=0.00001;
% end;
% Feature(:,3) = abs(Feature(:,3));


