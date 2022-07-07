 function data = FuncInputBA_malaga(PVector,Feature,K,GT_P0,input_value)


% Data.intrinsics
data.intrinsics = cameraIntrinsics([K(1,1), K(2,2)],[K(1,3), K(2,3)],[1280,720]); 

% Data.pointTracks
for i = 1 : size(Feature,1)
    points = [];
    viewIDs = [];    
    for j = 1 : Feature(i,5)
        viewIDs(end+1)  = Feature(i, 3+3*j); 
        points(end+1,:) = Feature(i, 4+3*j: 5+3*j);
    end 
    data.pointTracks(1,i) = pointTrack(viewIDs, points);
end 

% Data.xyzpoints
Pose = reshape(PVector.Pose,6,[])';
feat_pos = [];
tm_temp =[]; 
num_feat = size(Feature, 1); 

for i = 1:num_feat
    Xj = FuncXj(PVector.Feature(3*(i-1)+1),PVector.Feature(3*(i-1)+2));
    R = RMatrixYPR22(PVector.Pose(6*(Feature(i,3))-5),PVector.Pose(6*(Feature(i,3))-4),PVector.Pose(6*(Feature(i,3))-3));
    
    tm = PVector.Pose(6*((Feature(i,3)-1))+4:6*((Feature(i,3)-1))+6)';
    tm_temp = [tm_temp; tm];
    
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
    if feat_pos(:,i) - tm_temp(i,:)' == 0
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

        ViewId = [ViewId;uint32(i)];
        if strcmp(input_value, 'Ground_truth')
            AbsolutePose = [AbsolutePose; rigid3d(RMatrixYPR22(GT_P0(i,1),GT_P0(i,2),GT_P0(i,3)),GT_P0(i,4:6))];
        elseif strcmp(input_value, 'Estimated')
            AbsolutePose = [AbsolutePose; rigid3d(RMatrixYPR22(Pose(i,1),Pose(i,2),Pose(i,3)),Pose(i,4:6))];
        end

end
data.cameraPoses = table(ViewId, AbsolutePose);
