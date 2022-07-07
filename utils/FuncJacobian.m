
function [Jacobian] = FuncJacobian(xVector,PVector,Feature,K,ID1,ID2,nJacobian)

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

Jacobian = sparse(ID1,ID2,Val);
