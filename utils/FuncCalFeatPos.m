function feat_pos = FuncCalFeatPos(feat_ob, PVector, Feature)


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

