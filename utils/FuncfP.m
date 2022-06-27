
function [uvcomp] = FuncfP(xVector,PVector,Feature,K)

nRowNum = length(xVector.u);
uvcomp = zeros(nRowNum,1);

for j=1:nRowNum/2;
    PID = xVector.PID(j,1);
    FID = xVector.FID(j,1);
    
    if PID == Feature(FID,3);
        Pi = PVector.Pose(6*(PID-1)+1:6*PID)';
        Fj = PVector.Feature(3*FID-2:3*FID-1,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        RZ = [cos(Pi(1)),sin(Pi(1)),0;
              -sin(Pi(1)),cos(Pi(1)),0;
              0,0,1];
        RY = [cos(Pi(2)),0,-sin(Pi(2));
              0,1,0;
              sin(Pi(2)),0,cos(Pi(2))];
        RX = [1,0,0;
              0,cos(Pi(3)),sin(Pi(3));
              0,-sin(Pi(3)),cos(Pi(3))];
        Xj = [sin(Fj(1))*cos(Fj(2));
              sin(Fj(2));
              cos(Fj(1))*cos(Fj(2))];        
        R = RX*RY*RZ;
        x=K*R*Xj;
        u = [x(1)/x(3);
             x(2)/x(3)];        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        uvcomp(2*j-1:2*j) = u;

    elseif PID == Feature(FID,4);
        Pi = PVector.Pose(6*(Feature(FID,3)-1)+1:6*Feature(FID,3))';
        Pk = PVector.Pose(6*(PID-1)+1:6*PID)';
        Fj = PVector.Feature(3*FID-2:3*FID,1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        RZ = [cos(Pk(1)),sin(Pk(1)),0;
              -sin(Pk(1)),cos(Pk(1)),0;
              0,0,1];
        RY = [cos(Pk(2)),0,-sin(Pk(2));
              0,1,0;
              sin(Pk(2)),0,cos(Pk(2))];
        RX = [1,0,0;
              0,cos(Pk(3)),sin(Pk(3));
              0,-sin(Pk(3)),cos(Pk(3))];
        Xj = [sin(Fj(1))*cos(Fj(2));
              sin(Fj(2));
              cos(Fj(1))*cos(Fj(2))];        
        R = RX*RY*RZ;
        
        Tik = Pk(4:6)'-Pi(4:6)';
        Dik = norm(Tik);
        Omega2 = acos(Xj'*Tik/Dik);
        Xk = sin(Fj(3)+Omega2)*Dik*Xj-sin(Fj(3))*Tik;
        x=K*R*Xk;
        u = [x(1)/x(3);
             x(2)/x(3)];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        uvcomp(2*j-1:2*j) = u;
        
    else
        Pi = PVector.Pose(6*(Feature(FID,3)-1)+1:6*Feature(FID,3))';
        Pk = PVector.Pose(6*(Feature(FID,4)-1)+1:6*Feature(FID,4))';
        Pl = PVector.Pose(6*(PID-1)+1:6*PID)';
        Fj = PVector.Feature(3*FID-2:3*FID,1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        RZ = [cos(Pl(1)),sin(Pl(1)),0;
              -sin(Pl(1)),cos(Pl(1)),0;
              0,0,1];
        RY = [cos(Pl(2)),0,-sin(Pl(2));
              0,1,0;
              sin(Pl(2)),0,cos(Pl(2))];
        RX = [1,0,0;
              0,cos(Pl(3)),sin(Pl(3));
              0,-sin(Pl(3)),cos(Pl(3))];
        Xj = [sin(Fj(1))*cos(Fj(2));
              sin(Fj(2));
              cos(Fj(1))*cos(Fj(2))];        
        R = RX*RY*RZ;
        
        Tik = Pk(4:6)'-Pi(4:6)';        
        Dik = norm(Tik);
        Omega2 = acos(Xj'*Tik/Dik);
        Xl = sin(Fj(3)+Omega2)*Dik*Xj-sin(Fj(3))*(Pl(4:6)'-Pi(4:6)');
        x=K*R*Xl;
        u = [x(1)/x(3);
             x(2)/x(3)];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        uvcomp(2*j-1:2*j) = u;
        
    end;
end;
            
