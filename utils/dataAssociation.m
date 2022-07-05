function [state, state_k,groundT,measurement,indices,feat_ob] = dataAssociation(ImageNum, groundTruthStates, groundTruthMap, msckfState_imuOnly, measurements,camera)

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