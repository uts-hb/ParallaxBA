
function [u] = Funcui(K,Pi,Xj)

Alphai = Pi(1);   
Betai = Pi(2);
Gammai = Pi(3);

Phi = Xj(1);
Theta = Xj(2);

RZ = FuncRZ(Alphai);
RY = FuncRY(Betai);
RX = FuncRX(Gammai);
Xj = FuncXj(Phi,Theta);
R = FuncR(RX,RY,RZ);
x = Funcx(K,R,Xj);
u = Funcu(x);