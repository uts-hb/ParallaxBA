function [PVector,Feature] = FuncGetInitial3(PVector,Feature,Image,PoseID,K)

if PoseID == 1;
	PVector.Pose(6*(PoseID-1)+1:6*PoseID,1) = Image(1,1:6)';
else
    dR = RMatrixYPR22(Image(1,1),Image(1,2),Image(1,3));
    R1 = RMatrixYPR22(PVector.Pose(6*(PoseID-2)+1),PVector.Pose(6*(PoseID-2)+2),PVector.Pose(6*(PoseID-2)+3));
    R = dR*R1;
    [PVector.Pose(6*(PoseID-1)+1,1),PVector.Pose(6*(PoseID-1)+2,1),PVector.Pose(6*(PoseID-1)+3,1)] = InvRotMatrixYPR22(R);
    PVector.Pose(6*(PoseID-1)+4:6*PoseID,1) = PVector.Pose(6*(PoseID-2)+4:6*(PoseID-1),1)-R'*Image(1,4:6)';
end;

[nRowNum,nColumNum] = size(Image);
[nRowNumF,nColumNumF] = size(Feature);
 %ID,kind(2or3),M_anchor,A_anchor,Times,Pose_i,uv,...
for i=1:nRowNum-1;
    FID = Image(i+1,1);
    if FID > nRowNumF;
        Feature(FID,1) = FID;
        Feature(FID,2) = 2;
        Feature(FID,3) = PoseID;
        Feature(FID,5) = 1;
        Feature(FID,6) = PoseID;
        Feature(FID,7:8) = Image(i+1,2:3);
        V1 = Funcuv2V(Image(i+1,2:3),PVector.Pose(6*(PoseID-1)+1:6*(PoseID-1)+3,1),K);
        [PVector.Feature(3*FID-2,1),PVector.Feature(3*FID-1,1)] = FuncV2PT(V1);
%         clear V1;
        continue;
        
    elseif Feature(FID,5) == 0;
        Feature(FID,1) = FID;
        Feature(FID,2) = 2;
        Feature(FID,3) = PoseID;
        Feature(FID,5) = 1;
        Feature(FID,6) = PoseID;
        Feature(FID,7:8) = Image(i+1,2:3);
        V1 = Funcuv2V(Image(i+1,2:3),PVector.Pose(6*(PoseID-1)+1:6*(PoseID-1)+3),K);
        [PVector.Feature(3*FID-2,1),PVector.Feature(3*FID-1,1)] = FuncV2PT(V1);
%         clear V1;
        continue;
    
    elseif Feature(FID,5) == 1;
        Feature(FID,2) = 3;
        Feature(FID,4) = PoseID;
        Feature(FID,5) = 2;
        Feature(FID,9) = PoseID;
        Feature(FID,10:11) = Image(i+1,2:3);
        V1 = Funcuv2V(Feature(FID,7:8),PVector.Pose(6*(Feature(FID,3)-1)+1:6*(Feature(FID,3)-1)+3),K);
        V2 = Funcuv2V(Image(i+1,2:3),PVector.Pose(6*(PoseID-1)+1:6*(PoseID-1)+3),K);
        PVector.Feature(3*FID,1) = FuncV2O(V1,V2);
%         clear V1 V2;
        continue;
    
    elseif Feature(FID,5) >= 2;
        Omegai = [];
%         if PVector.Feature(3*FID,1)<0; 
        if PVector.Feature(3*FID,1)<0.02;           
            V2 = Funcuv2V(Image(i+1,2:3),PVector.Pose(6*(PoseID-1)+1:6*(PoseID-1)+3),K);
            for j=1:Feature(FID,5);            
                V1 = Funcuv2V(Feature(FID,3*(j-1)+7:3*(j-1)+8),PVector.Pose(6*(Feature(FID,3*(j-1)+6)-1)+1:6*(Feature(FID,3*(j-1)+6)-1)+3),K);
                Omegai(j) = FuncV2O(V1,V2);
            end;
            [Max,ID] = max(Omegai);
            if Max>PVector.Feature(3*FID,1)*1.2 || Max>0.02;
                PVector.Feature(3*FID,1) = Max;
                Feature(FID,3) = Feature(FID,3*(ID-1)+6);
                Feature(FID,4) = PoseID;
                V1 = Funcuv2V(Feature(FID,3*(ID-1)+7:3*(ID-1)+8),PVector.Pose(6*(Feature(FID,3)-1)+1:6*(Feature(FID,3)-1)+3),K);
                [PVector.Feature(3*FID-2,1),PVector.Feature(3*FID-1,1)] = FuncV2PT(V1);
            end;
        end;
    
        Feature(FID,3*Feature(FID,5)+6) = PoseID;
        Feature(FID,3*Feature(FID,5)+7:3*Feature(FID,5)+8) = Image(i+1,2:3);
        Feature(FID,5) = Feature(FID,5)+1;
%         clear V1 V2 Max ID Omegai;
        continue;
    
    end;
end;