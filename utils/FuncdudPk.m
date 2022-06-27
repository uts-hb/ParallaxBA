%% Feature j is projective at Pose i+k

function [dudRk,dudTi,dudTk,dudXj] = FuncdudPk(K,Pi,Pk,Xj)

Alphak = Pk(1);   
Betak = Pk(2);
Gammak = Pk(3);
Ti = Pi(4:6)';
Tk = Pk(4:6)';

Phi = Xj(1);
Theta = Xj(2);
Omega = Xj(3);
% Omega = abs(Xj(3));

RZ = FuncRZ(Alphak);
RY = FuncRY(Betak);
RX = FuncRX(Gammak);
R = FuncR(RX,RY,RZ);
Xj = FuncXj(Phi,Theta);
Tik = Tk-Ti;
D2 = FuncD2(Tik);
Dik = FuncDik(D2);
DotP = FuncDotP(Xj,Tik);
CosO2 = FuncCosO2(DotP,Dik);
Omega2 = FuncOmega2(CosO2);
SinO = FuncSinO(Omega,Omega2);
Xk = FuncXk(SinO,Dik,Xj,Omega,Tik);
x = Funcx(K,R,Xk);

dRZdA = FuncdRZdA(Alphak);
dRYdB = FuncdRYdB(Betak);
dRXdG = FuncdRXdG(Gammak);

dudx = Funcdudx(x);
dxdA = FuncdxdA(K,RX,RY,dRZdA,Xk);
dudA = dudx*dxdA;%%
dxdB = FuncdxdB(K,RX,dRYdB,RZ,Xk);
dudB = dudx*dxdB;%%
dxdG = FuncdxdG(K,dRXdG,RY,RZ,Xk);
dudG = dudx*dxdG;%%

dxdXk = FuncdxdXj(K,R);
dSinOdO = FuncdSinOdO(Omega,Omega2);
dXkdO = FuncdXkdO(dSinOdO,Dik,Xj,Omega,Tik);
dudO = dudx*dxdXk*dXkdO;%%

dSinOdO2 = dSinOdO;
dO2dCosO2 = FuncdO2dCosO2(CosO2);
dDotPdTik = FuncdDotPdTik(Xj);
dDikdD2 = FuncdDikdD2(D2);
dD2dTik = FuncdD2dTik(Tik);
dTikdTi = -eye(3);
dTikdTk = eye(3);
dDotPdTi = dDotPdTik*dTikdTi;
dDotPdTk = dDotPdTik*dTikdTk;
dDikdTi = dDikdD2*dD2dTik*dTikdTi;
dDikdTk = dDikdD2*dD2dTik*dTikdTk;
dCosO2dTi = FuncdCosO2dTi(dDotPdTi,DotP,dDikdTi,Dik);
dCosO2dTk = FuncdCosO2dTk(dDotPdTk,DotP,dDikdTk,Dik);
dSinOdTi = dSinOdO2*dO2dCosO2*dCosO2dTi;
dSinOdTk = dSinOdO2*dO2dCosO2*dCosO2dTk;
dXkdTi = FuncdXkdTi(Xj,dSinOdTi,Dik,dDikdTi,SinO,Omega,dTikdTi);
dXkdTk = FuncdXkdTk(Xj,dSinOdTk,Dik,dDikdTk,SinO,Omega,dTikdTk);

dudTi = dudx*dxdXk*dXkdTi;%%
dudTk = dudx*dxdXk*dXkdTk;%%

dCosO2dXj = FuncdCosO2dXj(Tik,Dik);
dXjdPT = FuncdXjdPT(Phi,Theta);
dSinOdPT = dSinOdO2*dO2dCosO2*dCosO2dXj*dXjdPT;
dXkdPT = FuncdXkdPT(Dik,Xj,dSinOdPT,SinO,dXjdPT);

dudPT = dudx*dxdXk*dXkdPT;%%


dudRk = cat(2,dudA,dudB,dudG);
dudXj = cat(2,dudPT,dudO);
