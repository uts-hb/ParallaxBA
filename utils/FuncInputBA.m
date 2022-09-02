function data = FuncInputBA(camera, state, groundT, feat_ob, measurement, input_value, feat_pos, refined_feat_pos, indices, Pose, temp_Pose)

data.intrinsics = cameraIntrinsics([camera.f_u, camera.f_v],[camera.c_u, camera.c_v],[640,480]);

ViewId = [];
AbsolutePose = [];

for i = 1:length(state.camState)
    if i == 1
        ViewId = [ViewId;uint32(i)];
        AbsolutePose = [AbsolutePose; rigid3d((quatToRotMat(groundT.camState{1,i}.q_C_G)),groundT.camState{1,i}.p_C_G')];
        
    else
        ViewId = [ViewId;uint32(i)];
        if strcmp(input_value, 'Ground_truth')  || strcmp(input_value, 'Initialization_4')
            AbsolutePose = [AbsolutePose;rigid3d((quatToRotMat(groundT.camState{1,i}.q_C_G)),groundT.camState{1,i}.p_C_G')];
        elseif strcmp(input_value, 'Estimated')
            AbsolutePose = [AbsolutePose;rigid3d((quatToRotMat(state.camState{1,i}.q_C_G)),state.camState{1,i}.p_C_G')];
        elseif strcmp(input_value, 'Initialization_5')
            Pose(:,4:6) = temp_Pose(:,4:6) + Pose(1,4:6);
            AbsolutePose = [AbsolutePose;rigid3d(RMatrixYPR22(Pose(i,1),Pose(i,2),Pose(i,3)), Pose(i,4:6))];
            
        end
    end
end

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
    if strcmp(input_value, 'Ground_truth')
        data.xyzPoints = [data.xyzPoints ; groundT.feature{1, feat_ob(1,i)}'] ;
    elseif strcmp(input_value, 'Estimated') || strcmp(input_value, 'Initialization_4')
        data.xyzPoints = [data.xyzPoints ; feat_pos(:,i)'] ;
    elseif strcmp(input_value, 'Initialization_5')
        data.xyzPoints = [data.xyzPoints ; refined_feat_pos(:,i)'] ;
    end
end
