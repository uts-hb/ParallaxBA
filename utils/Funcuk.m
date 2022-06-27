
function [u] = Funcuk(K,Pi,Pk,Xj)

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

% if abs(Omega)<1e-22 && abs(SinO)<1e-22;
%     Xk = Xj;
% else
    Xk = FuncXk(SinO,Dik,Xj,Omega,Tik);
% end;

x = Funcx(K,R,Xk);
u = Funcu(x);

% if Omega2<0||Omega2>pi/2;
%     Omega2
% end;
% 
% if Omega<0;
%     Omega
% end;