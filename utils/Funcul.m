
function [u] = Funcul(K,Pi,Pk,Pl,Xj)

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

% if abs(Omega)<1e-22 && abs(SinO)<1e-22;
%     Xl = Xj;
% else
    Xl = FuncXk(SinO,Dik,Xj,Omega,Til);
% end;

x = Funcx(K,R,Xl);
u = Funcu(x);

% if Omega2<0||Omega2>pi/2;
%     Omega2
% end;
% 
% if Omega<0;
%     Omega
% end;