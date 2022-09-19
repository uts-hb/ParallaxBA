close all;
clc;
clear;
addpath('utils');
load IMU_prop.mat

%%=========================Setup===============================%%
Data = 'starry_40';
load(strcat(Data,'.mat'));

% input_value = 'Ground_truth';
input_value = 'Estimated';

% Number of Image
ImageNum =50;

%Load Data for Parallax BA
K = textread('calSBA_starry.txt');         % intrinsic matrix

% Data Association using "Starry Night" dataset
[state, state_k, groundT, measurement, indices, feat_ob] = dataAssociation(ImageNum, groundTruthStates, groundTruthMap, msckfState_imuOnly, measurements,camera);


%%==========================Parallax BA========================%%

%% Creating Image File
FuncCreateImg(groundT, state, measurements,state_k,indices, input_value);

%% Loading file from the Image
GTName = strcat('DataPrepareBA/Starry/','GT_PO_PA.mat');
load(GTName);                           % Loading Ground Truth

%ID,kind(2or3),M_anchor,A_anchor,Times,Pose_i,uv,...
% Feature = [];
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
        %         temp = [temp;V1';Omega];
    end
end

%% Calculate Feature position from main anchor and associated anchor for input of SBA(input to Parallax BA)
feat_pos = FuncCalFeatPos(feat_ob, PVector, Feature);
true_feat = groundTruthMap(:, feat_ob);

% fprintf('Time Use %d\n\n', LoadTime);


figure('Name','Parallax BA with IMA');
Pose = reshape(PVector.Pose,6,[])';
plot3(Pose(:,4),Pose(:,5),Pose(:,6),'-*y');
hold on;
%% Least Squares
tic;
[PVector,Reason,Info,Sum_Error, Errors_par, reprojectionErrors_PBA_initial] = FuncLeastSquares_wIMU(xVector,PVector,Feature,K, camera, measurements);

%% Levenberg-Marquardt Iteration SBA
%         [PVector,Reason,Info] = FuncLeastSquaresLMSBA(xVector,PVector,Feature,K,FixVa);
BATime = toc;
fprintf('Time Use %d\n\n', BATime);

%% plotting the result of Parallax BA
Pose = reshape(PVector.Pose,6,[])';
FeaturePos = reshape(PVector.Feature,3,[])';
feat_pos = FuncCalFeatPos(feat_ob, PVector, Feature);

% figure('Name','Parallax BA with IMA');
plot3(GT_P0(:,4),GT_P0(:,5),GT_P0(:,6),'-+r');
% axis equal;
hold on;
grid on;
plot3(Pose(:,4),Pose(:,5),Pose(:,6),'-*g');
axis equal;

RMSE_feat_parallax = sqrt(mean(true_feat' - feat_pos').^2);
ARMSE_feat_parallax = mean(RMSE_feat_parallax(1,:));

RMSE_pose_parallax = sqrt(mean(GT_P0-Pose).^2);
ARMSE_rotation_parallax = mean(RMSE_pose_parallax(1,1:3));
ARMSE_position_parallax =  mean(RMSE_pose_parallax(1,4:6));
reprojectionErrors_PBA_final= (Errors_par{end}'*Errors_par{end})/(size(Errors_par{end},1)/2);   
%%=============================================================%%
