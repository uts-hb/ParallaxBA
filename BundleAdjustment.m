close all;
clc;
clear;
addpath('utils');

%% =============================Setup============================== %%
% method = 'parallaxBA';
method = 'standardBA';
% method = 'comparison'; % Run both BA at same time 
% % 
Data = 'starry_80';
% Data = 'starry_500';

% Number of Image 
ImageNum =200;

%Load Data for Parallax BA 
K = textread('calSBA_starry.txt');         % intrinsic matrix

% Load Data for Standard BA 
load(strcat(Data,'.mat'));
load IMU_prop.mat
%% =============================Data Association============================== %%
state_k = [500:499+ImageNum];         % State number of camera states
k_Start = state_k(1,1);
k_Second = state_k(1,2);
k_End = state_k(1,end);

% Input estimated camera poses into the state.camState vector
state.camState{1,1}.p_C_G = groundTruthStates{1, k_Start}.camState.p_C_G;
state.camState{1,1}.q_C_G =groundTruthStates{1, k_Start}.camState.q_CG;
state.camState{1,1}.eul = me_InvRotMatrixYPR22(quatToRotMat(state.camState{1,1}.q_C_G));
state.camState{1,1}.state_k = k_Start;
for i = 1 : k_End-k_Start
    C_IG = quatToRotMat(msckfState_imuOnly{1, i+k_Start}.imuState.q_IG);
    state.camState{1,i+1}.p_C_G = msckfState_imuOnly{1, i+k_Start}.imuState.p_I_G + C_IG' * camera.p_C_I;
    state.camState{1,i+1}.q_C_G = quatLeftComp(camera.q_CI) * msckfState_imuOnly{1, i+k_Start}.imuState.q_IG;
    state.camState{1,i+1}.eul = me_InvRotMatrixYPR22(quatToRotMat(state.camState{1,i+1}.q_C_G));
    state.camState{1,i+1}.state_k = i + k_Start;
end

% Ground Truth for Camera States and feature position
feature = [1:size(groundTruthMap, 2)];        % Observed feature position
for i = 1 : length(state_k)
    groundT.camState{1,i}.p_C_G = groundTruthStates{1, state_k(1,i)}.camState.p_C_G;
    groundT.camState{1,i}.q_C_G = groundTruthStates{1, state_k(1,i)}.camState.q_CG;
    groundT.camState{1,i}.eul = me_InvRotMatrixYPR22(quatToRotMat(groundT.camState{1,i}.q_C_G));
end
for i = 1: length(feature)
    groundT.feature{1,i} = groundTruthMap(:,feature(1,i)) ;
end

% Measurement from each image in camera coordinate system
% Indices cell indicates the camera pose ID that is observing each feature
measurement = cell(1, numel(feature));
indices =  cell(1, numel(feature)); 
for i =  1 : length(feature)
    for j = 1 : length(state_k)
        measurement{1,i}(:,end+1) = measurements{1,state_k(1,j)}.y(:,feature(1,i));
        %         measurement{1,i}(:,~all(isnan(measurement{1,i})));
    end
    indices{1,i} = find(~isnan(measurement{1,i}(1,:)) == 1);
    if length(indices{1,i}) < 2
        indices{1,i} = [];
    end
end

% Feature ID which is observed during the trajectory
feat_ob = [];
for i = 1: length(indices)
    if ~isempty(indices{1,i})
        feat_ob(:,end+1) = i;
    end
end
save feat_ob.mat feat_ob
%% =============================Parallax BA============================== %%
if strcmp(method, 'parallaxBA') | strcmp(method, 'comparison')
    
    %% Creating Image File 
    GT_P0 = [];
    for i = 1 : length(groundT.camState)
        GT_P0(end+1, :) = [me_InvRotMatrixYPR22(quatToRotMat(groundT.camState{1,i}.q_C_G)), groundT.camState{1,i}.p_C_G'];   
    end
    image = [];
    for i = 1:length(state_k)
        image = [];
        if i == 1
            image(1,:) = [me_InvRotMatrixYPR22(quatToRotMat(groundT.camState{1,i}.q_C_G)), groundT.camState{1,i}.p_C_G'];
            c = 1;
        else
%             image(1,:) = [me_InvRotMatrixYPR22(quatToRotMat(groundT.camState{1,i}.q_C_G)), groundT.camState{1,i}.p_C_G'];
            image(1,:) = [me_InvRotMatrixYPR22(quatToRotMat(state.camState{1,i}.q_C_G)), state.camState{1,i}.p_C_G'];
            c = 1;
        end
        for j = 1 : length(measurements{1,state_k(1,i)}.y)
            if ~isempty(indices{1,j})
                if ~isnan(measurements{1,state_k(1,i)}.y(1,j))
                    image(c+1,:) = [j, measurements{1,state_k(1,i)}.y(:,j)', zeros(1,3)];
                    c = c+1;
                end
            end
        end
        imageName = strcat('DataPrepareBA/image',int2str(i),'.mat');
        save(imageName, 'image');
    end

    save('DataPrepareBA/GT_PO_PA.mat', 'GT_P0')

    %% Loading file from the Image
    tic;
    GTName = strcat('DataPrepareBA/','GT_PO_PA.mat');
    load(GTName);                           % Loading Ground Truth
    
    %ID,kind(2or3),M_anchor,A_anchor,Times,Pose_i,uv,...
    % Feature = [];
    Feature = zeros(10000,180);
    xVector.u = []; xVector.PID = []; xVector.FID = [];
    PVector.Pose = []; PVector.Feature = []; PVector.ID = []; PVector.Info = sparse([]);
    
    for i=1:ImageNum;
        file=strcat('DataPrepareBA/','image',int2str(i));
        load(file);
        fprintf('%s\n', file);
        xVector = FuncGetxVector(xVector,image,i);
        [PVector,Feature] = FuncGetInitial3_02(PVector,Feature,image,i,K);
    end;
    
    input_PVector = PVector; 
    
    
    
    %%To calculate True parallax parameterized value from ground truth%%
    %%================================================================%%
%     temp = []; 
%     for i = 1 : length(feat_ob)
%         vector_1 = groundTruthMap(:,feat_ob(i))- GT_P0(Feature(feat_ob(i),3), 4:6)';    
%         [Phi,Theta] = FuncV2PT(vector_1);
%         V1 = [Phi,Theta]; 
%         vector_2 = groundTruthMap(:,feat_ob(i))- GT_P0(Feature(feat_ob(i),4), 4:6)';
%         Omega = FuncV2O(vector_1,vector_2);
%         
%         PVector.Feature(3*(feat_ob(i)-1)+1:3*(feat_ob(i)-1)+3) = [V1';Omega];  
% %         temp = [temp;V1';Omega]; 
%     end 

    %%================================================================%%


    
    
    
    
    %% Calculate Feature position from main anchor and associated anchor (input to Parallax BA)
    %%================================================================%%
    feat_pos = []; 
    for i = 1:length(feat_ob)
        Xj = FuncXj(PVector.Feature(3*(feat_ob(1,i)-1)+1),PVector.Feature(3*(feat_ob(1,i)-1)+2));  
        R = RMatrixYPR22(PVector.Pose(6*(Feature(feat_ob(1,i),3))-5),PVector.Pose(6*(Feature(feat_ob(1,i),3))-4),PVector.Pose(6*(Feature(feat_ob(1,i),3))-3));
        
        tm = PVector.Pose(6*((Feature(feat_ob(1,i),3)-1))+4:6*((Feature(feat_ob(1,i),3)-1))+6)'; 
        ta = PVector.Pose(6*((Feature(feat_ob(1,i),4)-1))+4:6*((Feature(feat_ob(1,i),4)-1))+6)'; 
        
        phi = acos(dot(Xj, ((ta-tm)/norm(ta-tm))'));
        depth = (sin(PVector.Feature(3*(feat_ob(1,i)-1)+3)+phi)*norm(ta-tm))/sin(PVector.Feature(3*(feat_ob(1,i)-1)+3));
        
        feat_pos(:,end+1)= depth*Xj + tm';
    end 
    
    true_feat = groundTruthMap(:, feat_ob);
     %%================================================================%%
    
  
    LoadTime = toc;
    fprintf('Time Use %d\n\n', LoadTime);
    
    %% Choose Variavle to Fix
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
    % PVector.Pose(10,1) = 1;
    tic;
    
    % Least Squares
    [PVector,Reason,Info,Sum_Error, Errors_par] = FuncLeastSquares(xVector,PVector,Feature,K,FixVa);

    %% Levenberg-Marquardt Iteration SBA
%         [PVector,Reason,Info] = FuncLeastSquaresLMSBA(xVector,PVector,Feature,K,FixVa);
    BATime = toc;
    fprintf('Time Use %d\n\n', BATime);
    
    Pose = reshape(PVector.Pose,6,[])';
    FeaturePos = reshape(PVector.Feature,3,[])';
    
    %% Calculate Feature position from main anchor and associated anchor (output from Parallax BA)
    %%================================================================%%
%     feat_pos = []; 
%     for i = 1:length(feat_ob)
%         Xj = FuncXj(PVector.Feature(3*(feat_ob(1,i)-1)+1),PVector.Feature(3*(feat_ob(1,i)-1)+2));  
%         R = RMatrixYPR22(PVector.Pose(6*(Feature(feat_ob(1,i),3))-5),PVector.Pose(6*(Feature(feat_ob(1,i),3))-4),PVector.Pose(6*(Feature(feat_ob(1,i),3))-3));
%         
%         tm = PVector.Pose(6*((Feature(feat_ob(1,i),3)-1))+4:6*((Feature(feat_ob(1,i),3)-1))+6)'; 
%         ta = PVector.Pose(6*((Feature(feat_ob(1,i),4)-1))+4:6*((Feature(feat_ob(1,i),4)-1))+6)'; 
%         
%         phi = acos(dot(Xj, ((ta-tm)/norm(ta-tm))'));
%         depth = (sin(PVector.Feature(3*(feat_ob(1,i)-1)+3)+phi)*norm(ta-tm))/sin(PVector.Feature(3*(feat_ob(1,i)-1)+3));
%         
%         feat_pos(:,end+1)= depth*Xj + tm';
%     end 
%     
%     true_feat = groundTruthMap(:, feat_ob);
     %%================================================================%%

    %% plotting the result of Parallax BA 
    figure('Name','Parallax BA');
    plot3(GT_P0(:,4),GT_P0(:,5),GT_P0(:,6),'-+r');
    axis equal;
    hold on;
    grid on;
    plot3(Pose(:,4),Pose(:,5),Pose(:,6),'-*g');
    axis equal;
%             legend('Ground Truth','BA');

%     plot3(feat_pos(1,:),feat_pos(2,:),feat_pos(3,:), '*b');
%     plot3(true_feat(1,:),true_feat(2,:),true_feat(3,:), '+r')
    
    RMSE_feat_parallax = sqrt(mean(true_feat' - feat_pos').^2); 
    RMSE_pose_parallax = sqrt(mean(GT_P0-Pose).^2);
    reprojectionErrors_parallax = (Errors_par{end}'*Errors_par{end})/(size(Errors_par{end},1)/2); 
end

% =============================Standard BA============================== %%
if strcmp(method, 'standardBA') | strcmp(method, 'comparison')

    data.intrinsics = cameraIntrinsics([camera.f_u, camera.f_v],[camera.c_u, camera.c_v],[640,480]);
    
   ViewId = []; 
   AbsolutePose = []; 
    for i = 1:length(state.camState)
        if i == 1
            ViewId = [ViewId;uint32(i)];
%             data.cameraPoses.AbsolutePose(i,1) = rigid3d((quatToRotMat(groundT.camState{1,i}.q_C_G)),groundT.camState{1,i}.p_C_G');
            AbsolutePose = [AbsolutePose; rigid3d((quatToRotMat(groundT.camState{1,i}.q_C_G)),groundT.camState{1,i}.p_C_G')];

        else
            ViewId = [ViewId;uint32(i)];
%             AbsolutePose = [AbsolutePose;rigid3d((quatToRotMat(groundT.camState{1,i}.q_C_G)),groundT.camState{1,i}.p_C_G')];
            AbsolutePose = [AbsolutePose;rigid3d((quatToRotMat(state.camState{1,i}.q_C_G)),state.camState{1,i}.p_C_G')];
%             AbsolutePose = [AbsolutePose;rigid3d(RMatrixYPR22(Pose(i,1),Pose(i,2),Pose(i,3)), Pose(i,4:6))];
        end
    end
%     AbsolutePose = AbsolutePose'
    data.cameraPoses = table(ViewId, AbsolutePose);
    
    % data.pointTrack = {};
    for i = 1 : length(feat_ob)
        points = [];
        viewIDs = [];
        viewIDs = indices{1,feat_ob(1,i)};
        for j = 1 : length(viewIDs)
            points(end+1,:) = measurement{1,feat_ob(1,i)}(:,viewIDs(1,j))';
        end
        data.pointTracks(1,i) = pointTrack(viewIDs, points); 
    end
    
    data.xyzPoints = [];
    
    for i = 1:length(feat_ob)
            data.xyzPoints = [data.xyzPoints ; groundT.feature{1, feat_ob(1,i)}'] ;
%         data.xyzPoints = [data.xyzPoints ; feat_pos(:,i)'] ;
        % Too pool result of feature estimation due to the lack
    end
    
    %% BA MATLAB function
    Errors_st = {}; 
     [xyzRefinedPoints,refinedPoses,reprojectionErrors_standard,Errors_st,iter] = ...
        bundleAdjustment(data.xyzPoints,data.pointTracks,data.cameraPoses,data.intrinsics,Errors_st,'FixedViewIDs',[1,10],'RelativeTolerance', 1e-10,'MaxIterations',500);

    %% Plotting Result 
    GT_P0 = [];
    for i = 1 : length(groundT.camState)
        GT_P0(i,:) = [groundT.camState{1, i}.eul, groundT.camState{1, i}.p_C_G'];
    end
 
    figure('Name','Standard BA');
    plot3(GT_P0(:,4),GT_P0(:,5),GT_P0(:,6),'-+r');
    % plot(GT_P0(:,4),GT_P0(:,6),'-+r');
    axis equal;
    
    Pose = [];
    for i = 1 : size(refinedPoses,1)
        Pose(i,:) = [me_InvRotMatrixYPR22(refinedPoses.AbsolutePose(i,1).Rotation), refinedPoses.AbsolutePose(i,1).Translation ];
    end
    hold on;
    grid on;
    plot3(Pose(:,4),Pose(:,5),Pose(:,6),'-*b');
    axis equal;
%     legend('Ground Truth','BA');

    RMSE_feat_standard = sqrt(mean((groundTruthMap(:, feat_ob(:))' - xyzRefinedPoints).^2));
    RMSE_pose_standard = sqrt(mean(Pose-GT_P0).^2);
    Errors_standard = [];
    for i = 1 : size(Errors_st{end},2)
        Errors_standard = [Errors_standard;Errors_st{end}(:,i)];
    end
    
    reprojectionErrors_final = (Errors_standard'*Errors_standard)/(size(Errors_standard,1)/2);
    %     plot3(xyzRefinedPoints(:,1),xyzRefinedPoints(:,2),xyzRefinedPoints(:,3), '*b');
    %     plot3(true_feat(1,:),true_feat(2,:),true_feat(3,:), '+r')
    Errors_standard = [];
    for i = 1 : size(Errors_st{1},2)
        Errors_standard = [Errors_standard;Errors_st{1}(:,i)];
    end
    reprojectionErrors_initial = (Errors_standard'*Errors_standard)/(size(Errors_standard,1)/2);
    
    xyzRefinedPoints =xyzRefinedPoints';
end

