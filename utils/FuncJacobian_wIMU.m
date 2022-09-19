
function [Jacobian] = FuncJacobian_wIMU(xVector,PVector,Feature,K,ID1,ID2,nJacobian,camera,measurements)

nRowNumu = length(xVector.u);
Val = zeros(1,nJacobian);
aa = 1;

for j=1:nRowNumu/2;
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
        dRZdA = [-sin(Pi(1)),cos(Pi(1)),0;
                 -cos(Pi(1)),-sin(Pi(1)),0;
                 0,0,0];
        dRYdB = [-sin(Pi(2)),0,-cos(Pi(2));
                 0,0,0;
                 cos(Pi(2)),0,-sin(Pi(2))];
        dRXdG = [0,0,0;
                 0,-sin(Pi(3)),cos(Pi(3));
                 0,-cos(Pi(3)),-sin(Pi(3))];
        dudx = [1/x(3),0,-x(1)/x(3)^2;
                0,1/x(3),-x(2)/x(3)^2;];
        dxdA = K*RX*RY*dRZdA*Xj;    
        dudA = dudx*dxdA;%%     
        dxdB = K*RX*dRYdB*RZ*Xj;
        dudB = dudx*dxdB;%%
        dxdG = K*dRXdG*RY*RZ*Xj;
        dudG = dudx*dxdG;%%
        dxdXj = K*R;
        dXjdPT = [cos(Fj(1))*cos(Fj(2)),-sin(Fj(1))*sin(Fj(2));
                  0,cos(Fj(2));
                  -sin(Fj(1))*cos(Fj(2)),-cos(Fj(1))*sin(Fj(2))];
        dudPT = dudx*dxdXj*dXjdPT;%%
        dudRi = cat(2,dudA,dudB,dudG);     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Val(aa:aa+9) = [dudRi(1,:),dudRi(2,:),dudPT(1,:),dudPT(2,:)];
        aa = aa+10;

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
        DotP = Xj'*Tik;
        Omega2 = acos(DotP/Dik);
        Xk = sin(Fj(3)+Omega2)*Dik*Xj-sin(Fj(3))*Tik;
        x=K*R*Xk;
        dRZdA = [-sin(Pk(1)),cos(Pk(1)),0;
                 -cos(Pk(1)),-sin(Pk(1)),0;
                 0,0,0];
        dRYdB = [-sin(Pk(2)),0,-cos(Pk(2));
                 0,0,0;
                 cos(Pk(2)),0,-sin(Pk(2))];
        dRXdG = [0,0,0;
                 0,-sin(Pk(3)),cos(Pk(3));
                 0,-cos(Pk(3)),-sin(Pk(3))];
        dudx = [1/x(3),0,-x(1)/x(3)^2;
                0,1/x(3),-x(2)/x(3)^2;];
        dxdA = K*RX*RY*dRZdA*Xk;    
        dudA = dudx*dxdA;%%     
        dxdB = K*RX*dRYdB*RZ*Xk;
        dudB = dudx*dxdB;%%
        dxdG = K*dRXdG*RY*RZ*Xk;
        dudG = dudx*dxdG;%%
        dxdXk = K*R;
        dSinOdO = cos(Fj(3)+Omega2);
        dXkdO = dSinOdO*Dik*Xj-cos(Fj(3))*Tik;
        dudO = dudx*dxdXk*dXkdO;%%
        dSinOdO2 = dSinOdO;
        dO2dCosO2 = -1/sqrt(1-(DotP/Dik)^2);
        dDotPdTik = Xj';
        dDikdD2 = 1/(2*Dik);
        dD2dTik = 2*Tik';
        dTikdTi = -eye(3);
        dTikdTk = eye(3);
        dDotPdTi = dDotPdTik*dTikdTi;
        dDotPdTk = dDotPdTik*dTikdTk;
        dDikdTi = dDikdD2*dD2dTik*dTikdTi;
        dDikdTk = dDikdD2*dD2dTik*dTikdTk;
        dCosO2dTi = (dDotPdTi*Dik-dDikdTi*DotP)/Dik^2;
        dCosO2dTk = (dDotPdTk*Dik-dDikdTk*DotP)/Dik^2;
        dSinOdTi = dSinOdO2*dO2dCosO2*dCosO2dTi;
        dSinOdTk = dSinOdO2*dO2dCosO2*dCosO2dTk;
        dXkdTi = Xj*(dSinOdTi*Dik+dDikdTi*sin(Fj(3)+Omega2))-sin(Fj(3))*dTikdTi;
        dXkdTk = Xj*(dSinOdTk*Dik+dDikdTk*sin(Fj(3)+Omega2))-sin(Fj(3))*dTikdTk;
        dudTi = dudx*dxdXk*dXkdTi;%%
        dudTk = dudx*dxdXk*dXkdTk;%%
        dCosO2dXj = Tik'/Dik;
        dXjdPT = [cos(Fj(1))*cos(Fj(2)),-sin(Fj(1))*sin(Fj(2));
                  0,cos(Fj(2));
                  -sin(Fj(1))*cos(Fj(2)),-cos(Fj(1))*sin(Fj(2))];
        
        dSinOdPT = dSinOdO2*dO2dCosO2*dCosO2dXj*dXjdPT;
        dXkdPT = Dik*(Xj*dSinOdPT+sin(Fj(3)+Omega2)*dXjdPT);
        dudPT = dudx*dxdXk*dXkdPT;%%
        dudRk = cat(2,dudA,dudB,dudG);
        dudXj = cat(2,dudPT,dudO);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Val(aa:aa+23) = [dudRk(1,:),dudRk(2,:),dudTi(1,:),dudTi(2,:),dudTk(1,:),dudTk(2,:),dudXj(1,:),dudXj(2,:)];
        aa = aa+24;
 
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
        Til = Pl(4:6)'-Pi(4:6)';
        Dik = norm(Tik);
        DotP = Xj'*Tik;
        Omega2 = acos(DotP/Dik);
        Xl = sin(Fj(3)+Omega2)*Dik*Xj-sin(Fj(3))*Til;
        x=K*R*Xl;
        dRZdA = [-sin(Pl(1)),cos(Pl(1)),0;
                 -cos(Pl(1)),-sin(Pl(1)),0;
                 0,0,0];
        dRYdB = [-sin(Pl(2)),0,-cos(Pl(2));
                 0,0,0;
                 cos(Pl(2)),0,-sin(Pl(2))];
        dRXdG = [0,0,0;
                 0,-sin(Pl(3)),cos(Pl(3));
                 0,-cos(Pl(3)),-sin(Pl(3))];
        dudx = [1/x(3),0,-x(1)/x(3)^2;
                0,1/x(3),-x(2)/x(3)^2;];
        dxdA = K*RX*RY*dRZdA*Xl;    
        dudA = dudx*dxdA;%%     
        dxdB = K*RX*dRYdB*RZ*Xl;
        dudB = dudx*dxdB;%%
        dxdG = K*dRXdG*RY*RZ*Xl;
        dudG = dudx*dxdG;%%
        dxdXl = K*R;
        dSinOdO = cos(Fj(3)+Omega2);
        dXldO = dSinOdO*Dik*Xj-cos(Fj(3))*Til;
        dudO = dudx*dxdXl*dXldO;%%
        dSinOdO2 = dSinOdO;
        dO2dCosO2 = -1/sqrt(1-(DotP/Dik)^2);
        dDotPdTik = Xj';
        dDikdD2 = 1/(2*Dik);
        dD2dTik = 2*Tik';
        dTikdTi = -eye(3);
        dTikdTk = eye(3);
        dTildTi = -eye(3);
        dTildTl = eye(3);
        dDotPdTi = dDotPdTik*dTikdTi;
        dDotPdTk = dDotPdTik*dTikdTk;
        dDikdTi = dDikdD2*dD2dTik*dTikdTi;
        dDikdTk = dDikdD2*dD2dTik*dTikdTk;
        dCosO2dTi = (dDotPdTi*Dik-dDikdTi*DotP)/Dik^2;
        dCosO2dTk = (dDotPdTk*Dik-dDikdTk*DotP)/Dik^2;
        dSinOdTi = dSinOdO2*dO2dCosO2*dCosO2dTi;
        dSinOdTk = dSinOdO2*dO2dCosO2*dCosO2dTk;
        dXldTi = Xj*(dSinOdTi*Dik+dDikdTi*sin(Fj(3)+Omega2))-sin(Fj(3))*dTildTi;
        dXldTk = Xj*(dSinOdTk*Dik+dDikdTk*sin(Fj(3)+Omega2));
        dXldTl = -sin(Fj(3))*dTildTl;
        dudTi = dudx*dxdXl*dXldTi;%%
        dudTk = dudx*dxdXl*dXldTk;%%
        dudTl = dudx*dxdXl*dXldTl;%%
        dCosO2dXj = Tik'/Dik;
        dXjdPT = [cos(Fj(1))*cos(Fj(2)),-sin(Fj(1))*sin(Fj(2));
                  0,cos(Fj(2));
                  -sin(Fj(1))*cos(Fj(2)),-cos(Fj(1))*sin(Fj(2))];     
        dSinOdPT = dSinOdO2*dO2dCosO2*dCosO2dXj*dXjdPT;
        dXldPT = Dik*(Xj*dSinOdPT+sin(Fj(3)+Omega2)*dXjdPT);
        dudPT = dudx*dxdXl*dXldPT;%%
        dudRl = cat(2,dudA,dudB,dudG);
        dudXj = cat(2,dudPT,dudO);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        Val(aa:aa+29) = [dudRl(1,:),dudRl(2,:),dudTi(1,:),dudTi(2,:),dudTk(1,:),dudTk(2,:),dudTl(1,:),dudTl(2,:),dudXj(1,:),dudXj(2,:)];
        aa = aa+30;
        
    end;
end;

R_CI = quatToRotMat(camera.q_CI);
t_CI = camera.p_C_I;


for i = 1:size(PVector.Pose,1)/6-1

    %% Linear velocity Function
    A = R_CI' * v_RMatrixYPR22(PVector.Pose(6*(i-1)+1:6*(i-1)+3,1)');
    B = PVector.Pose(6*i+4:6*i+6, 1) - PVector.Pose(6*(i-1)+4:6*(i-1)+6, 1) -  (R_CI' * v_RMatrixYPR22(PVector.Pose(6*i+1:6*i+3,1)'))' * t_CI;
    fv = (A*B + t_CI) / measurements{1,i}.dT;

    %% Angular velocity Function
    C = R_CI' * v_RMatrixYPR22(PVector.Pose(6*i+1:6*i+3,1)');
    fw = -(C * A' - eye(3))/measurements{1,i}.dT;
    fw = crossMatToVec(fw);

    %% Rotation matrix for step k
    Pi = PVector.Pose(6*(i-1)+1:6*(i-1)+3,1)';

    % Rotation Matrix
    RZk = [cos(Pi(1)),sin(Pi(1)),0;
        -sin(Pi(1)),cos(Pi(1)),0;
        0,0,1];
    RYk = [cos(Pi(2)),0,-sin(Pi(2));
        0,1,0;
        sin(Pi(2)),0,cos(Pi(2))];
    RXk = [1,0,0;
        0,cos(Pi(3)),sin(Pi(3));
        0,-sin(Pi(3)),cos(Pi(3))];
    Rk = RXk*RYk*RZk;

    % Derivative of Rotation matrix
    dRZdAk = [-sin(Pi(1)),cos(Pi(1)),0;
        -cos(Pi(1)),-sin(Pi(1)),0;
        0,0,0];
    dRYdBk = [-sin(Pi(2)),0,-cos(Pi(2));
        0,0,0;
        cos(Pi(2)),0,-sin(Pi(2))];
    dRXdGk = [0,0,0;
        0,-sin(Pi(3)),cos(Pi(3));
        0,-cos(Pi(3)),-sin(Pi(3))];

    %% Rotation matrix for step k+1
    Pi = PVector.Pose(6*i+1:6*i+3,1)';

    % Rotation Matrix
    RZk1 = [cos(Pi(1)),sin(Pi(1)),0;
        -sin(Pi(1)),cos(Pi(1)),0;
        0,0,1];
    RYk1 = [cos(Pi(2)),0,-sin(Pi(2));
        0,1,0;
        sin(Pi(2)),0,cos(Pi(2))];
    RXk1 = [1,0,0;
        0,cos(Pi(3)),sin(Pi(3));
        0,-sin(Pi(3)),cos(Pi(3))];
    Rk1 = RXk1*RYk1*RZk1;

    % Derivative of Rotation matrix
    dRZdAk1 = [-sin(Pi(1)),cos(Pi(1)),0;
        -cos(Pi(1)),-sin(Pi(1)),0;
        0,0,0];
    dRYdBk1 = [-sin(Pi(2)),0,-cos(Pi(2));
        0,0,0;
        cos(Pi(2)),0,-sin(Pi(2))];
    dRXdGk1 = [0,0,0;
        0,-sin(Pi(3)),cos(Pi(3));
        0,-cos(Pi(3)),-sin(Pi(3))];


    %% Jacobian for linear velocity
    %dfv/dr_c_k
    dfvdAk = (R_CI' * (RXk*RYk*dRZdAk)* B) / measurements{1,i}.dT;
    dfvdBk = (R_CI' * (RXk*dRYdBk*RZk)* B) / measurements{1,i}.dT;
    dfvdGk = (R_CI' * (dRXdGk*RYk*RZk)* B) / measurements{1,i}.dT;

    %dfv/dt_c_k
    dfvdtk = -A / measurements{1,i}.dT;

    %dfv/dr_c_k+1
    dfvdAk1 = -(A * (R_CI' * (RXk1*RYk1*dRZdAk1))' * t_CI) / measurements{1,i}.dT;
    dfvdBk1 = -(A * (R_CI' * (RXk1*dRYdBk1*RZk1))' * t_CI) / measurements{1,i}.dT;
    dfvdGk1 = -(A * (R_CI' * (dRXdGk1*RYk1*RZk1))' * t_CI) / measurements{1,i}.dT;

    %dfv/dt_c_k+1
    dfvdtk1 = A / measurements{1,i}.dT;

    %% Jacobian for angular velocity

    %dfw/dr_c_k
    dfwdAk = -(C * (R_CI' * (RXk*RYk*dRZdAk))') / measurements{1,i}.dT;
    dfwdAk = crossMatToVec(dfwdAk);
    dfwdBk = -(C * (R_CI' * (RXk*dRYdBk*RZk))') / measurements{1,i}.dT;
    dfwdBk = crossMatToVec(dfwdBk);
    dfwdGk = -(C * (R_CI' * (dRXdGk*RYk*RZk))') / measurements{1,i}.dT;
    dfwdGk = crossMatToVec(dfwdGk);

    %dfw/dr_c_k+1
    dfwdAk1 = -(R_CI' * (RXk1*RYk1*dRZdAk1) * A') / measurements{1,i}.dT;
    dfwdAk1 = crossMatToVec(dfwdAk1);
    dfwdBk1 = -(R_CI' * (RXk1*dRYdBk1*RZk1) * A') / measurements{1,i}.dT;
    dfwdBk1 = crossMatToVec(dfwdBk1);
    dfwdGk1 = -(R_CI' * (dRXdGk1*RYk1*RZk1) * A') / measurements{1,i}.dT;
    dfwdGk1 = crossMatToVec(dfwdGk1);

    %% Jacobian vector

    Jv = [dfvdAk, dfvdBk, dfvdGk, dfvdtk, dfvdAk1, dfvdBk1, dfvdGk1, dfvdtk1];
    Val(aa:aa+35) = [Jv(1,:), Jv(2,:), Jv(3,:)]; 
    aa = aa + 36;
    Jw = [dfwdAk, dfwdBk, dfwdGk, dfwdAk1, dfwdBk1, dfwdGk1];
    Val(aa:aa+17) = [Jw(1,:), Jw(2,:), Jw(3,:)]; 
    aa = aa + 18; 

end

Jacobian = sparse(ID1,ID2,Val);
