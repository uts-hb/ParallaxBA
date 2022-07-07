close all;
clc;
clear;
addpath('utils');
GTName = strcat('DataPrepareBA/Malaga/','GT_P0_PA.mat');
load(GTName);

%% =========================Setup=============================== %%
%intrinsic matrix for "MALAGA" dataset 
K = textread('calSBA_malaga.txt');

% Number of Image (1~170)
ImageNum=20;

% input_value = 'Ground_truth';
input_value = 'Estimated';

%% ==========================Parallax BA======================== %%
%ID,kind(2or3),M_anchor,A_anchor,Times,Pose_i,uv,...
% Feature = [];
Feature = zeros(10000,180);
xVector.u = []; xVector.PID = []; xVector.FID = [];
PVector.Pose = []; PVector.Feature = []; PVector.ID = []; PVector.Info = sparse([]);

%% Load Files
tic;
for i=1:ImageNum;
    file=strcat('DataPrepareBA/Malaga/Image',int2str(i));
    load(file);
    fprintf('%s\n', file);
    xVector = FuncGetxVector(xVector,Image,i);
    [PVector,Feature] = malaga_FuncGetInitial3_02(PVector,Feature,Image,i,K,input_value,GT_P0);
end;
input_PVector = PVector; 

%% Deleting the Feature only observed once 
Feature = Feature(find(Feature(:,1) ~=0),:);
input_Feature = Feature; 
delete_Id = find(Feature(:,4) ==0); 
Feature = Feature(1:delete_Id(1)-1,:);
PVector.Feature = PVector.Feature(1:3*(delete_Id(1)-1)); 
xVector.PID = xVector.PID(1:length(xVector.PID)- length(delete_Id),:); 
xVector.FID = xVector.FID(1:end-length(delete_Id),:); 
xVector.u = xVector.u(1:end-2*length(delete_Id),:); 

%% Choose Variavle to Fix
% FixVa = 4;
FixVa = 3; % Fix Z
% PVector.Pose(12,1) = PVector.Pose(6,1)+1;
% FixVa = 2; % Fix Y
% PVector.Pose(11,1) = 1;
% FixVa = 1; % Fix X
% PVector.Pose(10,1) = 1;

tic;
% Least Squares
[PVector,Reason,Info,Sum_Error, Errors_par, reprojectionErrors_PBA_initial] = FuncLeastSquares(xVector,PVector,Feature,K,FixVa);
%
%% Levenberg-Marquardt Iteration SBA
%     [PVector,Reason,Info] = FuncLeastSquaresLMSBA(xVector,PVector,Feature,K,FixVa);
output_PVector = PVector; 
BATime = toc;
fprintf('Time Use %d\n\n', BATime);


%% Print Reason

fprintf('Reason is %d\n', Reason);

PVector.ID = Feature(:,1:4);
PVector.Info = Info;

%% plotting the result of Parallax BA
Pose = reshape(PVector.Pose,6,[])';
% FeaturePos = reshape(PVector.Feature,3,[])';

figure('Name','Parallax BA');
plot3(GT_P0(:,4),GT_P0(:,5),GT_P0(:,6),'-+r');
axis equal;
hold on;
grid on;
plot3(Pose(:,4),Pose(:,5),Pose(:,6),'-*g');
axis equal;

%% Calculating Reprojection Error and RMSE
% RMSE_feat_parallax = sqrt(mean(true_feat' - feat_pos').^2);
RMSE_pose_parallax = sqrt(mean(GT_P0(1:ImageNum,:)-Pose).^2);
reprojectionErrors_PBA_final= (Errors_par{end}'*Errors_par{end})/(size(Errors_par{end},1)/2);


%% ==========================Standard BA======================== %%
PVector = input_PVector; 
Feature = input_Feature; 

%% Converting data suitable for Standard BA MATLAB 
data = FuncInputBA_malaga(PVector,Feature,K,GT_P0,input_value);

%% BA MATLAB function
Errors_st = {};
tic;
[xyzRefinedPoints,refinedPoses,reprojectionErrors_standard,Errors_st,iter] = ...
    bundleAdjustment(data.xyzPoints,data.pointTracks,data.cameraPoses,data.intrinsics,Errors_st,'FixedViewIDs',[1,5],'RelativeTolerance', 1e-10,'MaxIterations',500);
toc

%% Plotting Result of Standard BA
figure('Name','Standard BA');
plot3(GT_P0(:,4),GT_P0(:,5),GT_P0(:,6),'-+r');
% plot(GT_P0(:,4),GT_P0(:,6),'-+r');

Pose = [];
for i = 1 : size(refinedPoses,1)
    Pose(i,:) = [v_InvRotMatrixYPR22(refinedPoses.AbsolutePose(i,1).Rotation), refinedPoses.AbsolutePose(i,1).Translation ];
end

hold on;
grid on;
plot3(Pose(:,4),Pose(:,5),Pose(:,6),'-*b');
axis equal;


%% Calculating Reprojection Error and RMSE
RMSE_pose_standard = sqrt(mean(Pose-GT_P0(1:ImageNum,:)).^2);

Errors_standard = [];
for i = 1 : size(Errors_st{end},2)
    Errors_standard = [Errors_standard;Errors_st{end}(:,i)];
end
reprojectionErrors_SBA_final = (Errors_standard'*Errors_standard)/(size(Errors_standard,1)/2);

Errors_standard = [];
for i = 1 : size(Errors_st{1},2)
    Errors_standard = [Errors_standard;Errors_st{1}(:,i)];
end
reprojectionErrors_SBA_initial = (Errors_standard'*Errors_standard)/(size(Errors_standard,1)/2);








