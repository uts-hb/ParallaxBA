close all;
clc;
clear;
addpath('utils');
load IMU_prop.mat


%% =============================Setup================================== %%

%% "Starry Night" dataset with 80 features and 500 features
% Data = 'starry_80';
Data = 'starry_500';
load(strcat(Data,'.mat'));

%% 4 different initialization methods
% input_value = 'Ground_truth';              % Initialization 3
input_value = 'Estimated';               % Initialization 6
% input_value = 'Initialization_4';       
% input_value = 'Initialization_5';        

%% Number of Image (1~500)
ImageNum =20;

%% Intrinsic matrix for "Starry Night" dataset
K = textread('calSBA_starry.txt');

%% Data Association using "Starry Night" dataset
[state, state_k, groundT, measurement, indices, feat_ob] = dataAssociation(ImageNum, groundTruthStates, groundTruthMap, msckfState_imuOnly, measurements,camera);


%% ==============================Parallax BA============================ %%

%% Creating Image File
FuncCreateImg(groundT, state, measurements,state_k,indices, input_value);


%% Loading data from the Image file
tic;
GTName = strcat('DataPrepareBA/Starry/','GT_PO_PA.mat');
load(GTName);                           % Loading Ground Truth

Feature = zeros(10000,180);
xVector.u = []; xVector.PID = []; xVector.FID = [];
PVector.Pose = []; PVector.Feature = []; PVector.ID = []; PVector.Info = sparse([]);


for i=1:ImageNum;
    file=strcat('DataPrepareBA/Starry/','Image',int2str(i));
    load(file);
    fprintf('%s\n', file);
    xVector = FuncGetxVector(xVector,image,i);
    [PVector,Feature] = FuncGetInitial3_02(PVector,Feature,image,i,K);
end;


%% To calculate True parallax parameterized value from ground truth%%
if strcmp(input_value, 'Ground_truth')
    temp = [];
    for i = 1 : length(feat_ob)
        vector_1 = groundTruthMap(:,feat_ob(i))- GT_P0(Feature(feat_ob(i),3), 4:6)';
        [Phi,Theta] = FuncV2PT(vector_1);
        V1 = [Phi,Theta];
        vector_2 = groundTruthMap(:,feat_ob(i))- GT_P0(Feature(feat_ob(i),4), 4:6)';
        Omega = FuncV2O(vector_1,vector_2);
        
        PVector.Feature(3*(feat_ob(i)-1)+1:3*(feat_ob(i)-1)+3) = [V1';Omega];
    end
end

%% Calculate Feature position with main anchor and associated anchor (Used for feature position input for SBA)
feat_pos = FuncCalFeatPos(feat_ob, PVector, Feature);
true_feat = groundTruthMap(:, feat_ob);

LoadTime = toc;
fprintf('Time Use %d\n\n', LoadTime);


%% Choose Variable to Fix
if ImageNum > 9
    FixVa = 4;
else
    FixVa = 5;
end
%     FixVa = 3; % Fix Z
%     PVector.Pose(12,1) = PVector.Pose(6,1)+1;
%     FixVa = 2; % Fix Y
%     PVector.Pose(11,1) = 1;
%     FixVa = 1; % Fix X
%     PVector.Pose(10,1) = 1;
tic;


%% Least Squares
[PVector,Reason,Info,Sum_Error, Errors_par, reprojectionErrors_PBA_initial] = FuncLeastSquares(xVector,PVector,Feature,K,FixVa);


%% Levenberg-Marquardt
%         [PVector,Reason,Info] = FuncLeastSquaresLMSBA(xVector,PVector,Feature,K,FixVa);

BATime = toc;
fprintf('Time Use %d\n\n', BATime);


%% plotting the result of Parallax BA
Pose = reshape(PVector.Pose,6,[])';
Pose(:,4:6) = Pose(:,4:6) - Pose(1,4:6);
GT_P0_ = GT_P0;
GT_P0_(:,4:6) = GT_P0_(:,4:6) - GT_P0_(1,4:6); 

if ImageNum > 9
    Pose(:,4:6) = (Pose(:,4:6)/Pose(10,6))*GT_P0_(10,6);
else
    Pose(:,4:6) = (Pose(:,4:6)/Pose(10,6))*GT_P0_(10,6);
end

FeaturePos = reshape(PVector.Feature,3,[])';

figure('Name','Parallax BA');
plot3(GT_P0_(:,4),GT_P0_(:,5),GT_P0_(:,6),'-+r');
axis equal;
hold on;
grid on;
plot3(Pose(:,4),Pose(:,5),Pose(:,6),'-*g');
axis equal;


%% Calculating Reprojection Error and RMSE for PBA
refined_feat_pos = FuncCalFeatPos(feat_ob, PVector, Feature);
RMSE_feat_parallax = sqrt(mean(true_feat' - refined_feat_pos').^2);
ARMSE_feat_parallax = mean(RMSE_feat_parallax(1,:));
RMSE_pose_parallax = sqrt(mean(GT_P0-Pose).^2);
ARMSE_rotation_parallax = mean(RMSE_pose_parallax(1,1:3));
ARMSE_position_parallax =  mean(RMSE_pose_parallax(1,4:6));
reprojectionErrors_PBA_final= (Errors_par{end}'*Errors_par{end})/(size(Errors_par{end},1)/2);


%% ==============================Standard BA============================ %%

%% Converting data suitable for Standard BA MATLAB 
data = FuncInputBA(camera, state, groundT, feat_ob, measurement, input_value, feat_pos, refined_feat_pos, indices, Pose);



%% BA MATLAB function
Errors_st = {};
[xyzRefinedPoints,refinedPoses,reprojectionErrors_standard,Errors_st,iter] = ...
    bundleAdjustment(data.xyzPoints,data.pointTracks,data.cameraPoses,data.intrinsics,Errors_st,'FixedViewIDs',[1,5],'RelativeTolerance', 1e-10,'MaxIterations',500);


%% Plotting Result of Standard BA
figure('Name','Standard BA');
plot3(GT_P0(:,4),GT_P0(:,5),GT_P0(:,6),'-+r');
% plot(GT_P0(:,4),GT_P0(:,6),'-+r');
axis equal;

Pose = [];
for i = 1 : size(refinedPoses,1)
    Pose(i,:) = [v_InvRotMatrixYPR22(refinedPoses.AbsolutePose(i,1).Rotation), refinedPoses.AbsolutePose(i,1).Translation ];
end

hold on;
grid on;
plot3(Pose(:,4),Pose(:,5),Pose(:,6),'-*b');
axis equal;


%% Calculating Reprojection Error and RMSE for SBA

RMSE_feat_standard = sqrt(mean((groundTruthMap(:, feat_ob(:))' - xyzRefinedPoints).^2));
ARMSE_feat_standard = mean(RMSE_feat_standard(1,:));
RMSE_pose_standard = sqrt(mean(Pose-GT_P0).^2);
ARMSE_rotation_standard = mean(RMSE_pose_standard(1,1:3));
ARMSE_position_standard =  mean(RMSE_pose_standard(1,4:6));

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

