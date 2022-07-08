%% Feature j is projective at Pose i+k+l

function [dudRl,dudTi,dudTk,dudTl,dudXj] = FuncdudPl(K,Pi,Pk,Pl,Xj)

Alphal = Pl(1);   
Betal = Pl(2);
Gammal = Pl(3);
Ti = Pi(4:6)';
Tk = Pk(4:6)';
Tl = Pl(4:6)';

Phi = Xj(1);
Theta = Xj(2);
Omega = Xj(3);

RZ = FuncRZ(Alphal);
RY = FuncRY(Betal);
RX = FuncRX(Gammal);
R = FuncR(RX,RY,RZ);
Xj = FuncXj(Phi,Theta);
Tik = Tk-Ti;
Til = Tl-Ti;
D2 = FuncD2(Tik);
Dik = FuncDik(D2);
DotP = FuncDotP(Xj,Tik);
CosO2 = FuncCosO2(DotP,Dik);
Omega2 = FuncOmega2(CosO2);
SinO = FuncSinO(Omega,Omega2);
Xl = FuncXk(SinO,Dik,Xj,Omega,Til);
x = Funcx(K,R,Xl);

dRZdA = FuncdRZdA(Alphal);
dRYdB = FuncdRYdB(Betal);
dRXdG = FuncdRXdG(Gammal);

dudx = Funcdudx(x);
dxdA = FuncdxdA(K,RX,RY,dRZdA,Xl);
dudA = dudx*dxdA;%%
dxdB = FuncdxdB(K,RX,dRYdB,RZ,Xl);
dudB = dudx*dxdB;%%
dxdG = FuncdxdG(K,dRXdG,RY,RZ,Xl);
dudG = dudx*dxdG;%%

dxdXl = FuncdxdXj(K,R);
dSinOdO = FuncdSinOdO(Omega,Omega2);
dXldO = FuncdXkdO(dSinOdO,Dik,Xj,Omega,Til);
dudO = dudx*dxdXl*dXldO;%%

dSinOdO2 = dSinOdO;
dO2dCosO2 = FuncdO2dCosO2(CosO2);
dDotPdTik = FuncdDotPdTik(Xj);
dDikdD2 = FuncdDikdD2(D2);
dD2dTik = FuncdD2dTik(Tik);
dTikdTi = -eye(3);
dTikdTk = eye(3);
dTildTi = -eye(3);
dTildTl = eye(3);
dDotPdTi = dDotPdTik*dTikdTi;
dDotPdTk = dDotPdTik*dTikdTk;
dDikdTi = dDikdD2*dD2dTik*dTikdTi;
dDikdTk = dDikdD2*dD2dTik*dTikdTk;
dCosO2dTi = FuncdCosO2dTi(dDotPdTi,DotP,dDikdTi,Dik);
dCosO2dTk = FuncdCosO2dTk(dDotPdTk,DotP,dDikdTk,Dik);
dSinOdTi = dSinOdO2*dO2dCosO2*dCosO2dTi;
dSinOdTk = dSinOdO2*dO2dCosO2*dCosO2dTk;
dXldTi = FuncdXkdTi(Xj,dSinOdTi,Dik,dDikdTi,SinO,Omega,dTildTi);
dXldTk = FuncdXldTk(Xj,dSinOdTk,Dik,dDikdTk,SinO);
dXldTl = FuncdXldTl(Omega,dTildTl);

dudTi = dudx*dxdXl*dXldTi;%%
dudTk = dudx*dxdXl*dXldTk;%%
dudTl = dudx*dxdXl*dXldTl;%%

dCosO2dXj = FuncdCosO2dXj(Tik,Dik);
dXjdPT = FuncdXjdPT(Phi,Theta);
dSinOdPT = dSinOdO2*dO2dCosO2*dCosO2dXj*dXjdPT;
dXldPT = FuncdXkdPT(Dik,Xj,dSinOdPT,SinO,dXjdPT);

dudPT = dudx*dxdXl*dXldPT;%%


dudRl = cat(2,dudA,dudB,dudG);
dudXj = cat(2,dudPT,dudO);
