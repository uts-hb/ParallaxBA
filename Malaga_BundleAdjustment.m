%% BundleAdjustment.m

close all;
clc;
clear;
tic;

%% =========================Setup=============================== %%
K = textread('calSBA_malaga.txt');

ImageNum=30;

input_value = 'Ground_truth';
% input_value = 'Estimated';

GTName = strcat('DataPrepareBA/Malaga/','GT_P0_PA.mat');
load(GTName);

%ID,kind(2or3),M_anchor,A_anchor,Times,Pose_i,uv,...
% Feature = [];
Feature = zeros(2500,180);
xVector.u = []; xVector.PID = []; xVector.FID = [];
PVector.Pose = []; PVector.Feature = []; PVector.ID = []; PVector.Info = sparse([]);


%% ==========================Parallax BA======================== %%
%% Load Files

for i=1:ImageNum;
    file=strcat('DataPrepareBA/Malaga/Image',int2str(i));
    load(file);
    fprintf('%s\n', file);
    xVector = FuncGetxVector(xVector,Image,i);
    [PVector,Feature] = malaga_FuncGetInitial3_02(PVector,Feature,Image,i,K,input_value,GT_P0);
end;
input_PVector = PVector; 
LoadTime = toc;
fprintf('Time Use %d\n\n', LoadTime);

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
%% ============================================================= %%

%% ==========================Standard BA======================== %%
num_feat = size(PVector.ID,1);

Feature = zeros(60000,180);
xVector.u = []; xVector.PID = []; xVector.FID = [];
PVector.Pose = []; PVector.Feature = []; PVector.ID = []; PVector.Info = sparse([]);

%% Converting data suitable for Standard BA MATLAB 

% Data.intrinsics
data.intrinsics = cameraIntrinsics([K(1,1), K(2,2)],[K(1,3), K(2,3)],[1024,768]);

% Data.pointTracks
points = [0,0];
viewIDs = [-1];
for i = 1 : num_feat
    data.pointTracks(1,i) = pointTrack(viewIDs,points);
end

for i=1:ImageNum;
    file=strcat('DataPrepareBA/Malaga/Image',int2str(i));
    load(file);
    fprintf('%s\n', file);
    xVector = FuncGetxVector(xVector,Image,i);
    [PVector,Feature] = malaga_FuncGetInitial3_02(PVector,Feature,Image,i,K,input_value,GT_P0);
    
    for j = 1 : size(Image,1)-1
        
        if data.pointTracks(1,Image(j+1,1)).ViewIds(1,1) == 0
            data.pointTracks(1,Image(j+1,1)).ViewIds(end) = i;
            data.pointTracks(1,Image(j+1,1)).Points = Image(j+1,2:3);
        else
            data.pointTracks(1,Image(j+1,1)).ViewIds(end+1) = i;
            data.pointTracks(1,Image(j+1,1)).Points = [data.pointTracks(1,Image(j+1,1)).Points; Image(j+1,2:3)];
        end
    end
    
end;

Pose = reshape(PVector.Pose,6,[])';

feat_pos = [];

for i = 1:num_feat
    Xj = FuncXj(PVector.Feature(3*(i-1)+1),PVector.Feature(3*(i-1)+2));
    R = RMatrixYPR22(PVector.Pose(6*(Feature(i,3))-5),PVector.Pose(6*(Feature(i,3))-4),PVector.Pose(6*(Feature(i,3))-3));
    
    tm = PVector.Pose(6*((Feature(i,3)-1))+4:6*((Feature(i,3)-1))+6)';
    if Feature(i,4) ~= 0
        ta = PVector.Pose(6*((Feature(i,4)-1))+4:6*((Feature(i,4)-1))+6)';
    else
        ta = tm;
    end
    
    phi = acos(dot(Xj, ((ta-tm)/norm(ta-tm))'));
    
    if Feature(i,4) ~= 0
        depth = (sin(PVector.Feature(3*(i-1)+3)+phi)*norm(ta-tm))/sin(PVector.Feature(3*(i-1)+3));
    else
        depth = 0;
    end
    
    feat_pos(:,end+1)= depth*Xj + tm';
end

temp_id = [];
for i = 1 : length(feat_pos)
    if feat_pos(:,i) - tm' == 0
        temp_id(end+1)=i;
    end
end

if ~isempty(temp_id)
    feat_pos(:,temp_id) = [];
    data.pointTracks =data.pointTracks(1:temp_id(1,1)-1);
end
data.xyzPoints = feat_pos' ;

% Data.cameraposes
ViewId = [];
AbsolutePose = [];
for i = 1:size(Pose,1)
    
    if i == 1
        ViewId = [ViewId;uint32(i)];
        AbsolutePose = [AbsolutePose; rigid3d(RMatrixYPR22(GT_P0(i,1),GT_P0(i,2),GT_P0(i,3)),GT_P0(i,4:6))];
        
    else
        ViewId = [ViewId;uint32(i)];
        if strcmp(input_value, 'Ground_truth')
            AbsolutePose = [AbsolutePose; rigid3d(RMatrixYPR22(GT_P0(i,1),GT_P0(i,2),GT_P0(i,3)),GT_P0(i,4:6))];
        end
        if strcmp(input_value, 'Estimated')
            AbsolutePose = [AbsolutePose;rigid3d(RMatrixYPR22(Pose(i,1),Pose(i,2),Pose(i,3)),Pose(i,4:6))];
            %             AbsolutePose = [AbsolutePose;rigid3d(RMatrixYPR22(Pose(i,1),Pose(i,2),Pose(i,3)), Pose(i,4:6))];
        end
    end
end
data.cameraPoses = table(ViewId, AbsolutePose);


%% BA MATLAB function
Errors_st = {};
tic;
[xyzRefinedPoints,refinedPoses,reprojectionErrors_standard,Errors_st,iter] = ...
    bundleAdjustment(data.xyzPoints,data.pointTracks,data.cameraPoses,data.intrinsics,Errors_st,'FixedViewIDs',[1],'RelativeTolerance', 1e-10,'MaxIterations',500);
toc

%% Plotting Result of Standard BA
figure('Name','Standard BA');
plot3(GT_P0(:,4),GT_P0(:,5),GT_P0(:,6),'-+r');
% plot(GT_P0(:,4),GT_P0(:,6),'-+r');

Pose = [];
for i = 1 : size(refinedPoses,1)
    Pose(i,:) = [me_InvRotMatrixYPR22(refinedPoses.AbsolutePose(i,1).Rotation), refinedPoses.AbsolutePose(i,1).Translation ];
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

%% ============================================================= %%







