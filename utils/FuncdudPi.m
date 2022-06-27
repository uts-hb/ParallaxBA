%% Feature j is projective at Pose i

function [dudRi,dudPT] = FuncdudPi(K,Pi,Xj)

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

dRZdA = FuncdRZdA(Alphai);
dRYdB = FuncdRYdB(Betai);
dRXdG = FuncdRXdG(Gammai);

dudx = Funcdudx(x);
dxdA = FuncdxdA(K,RX,RY,dRZdA,Xj);
dudA = dudx*dxdA;%%
dxdB = FuncdxdB(K,RX,dRYdB,RZ,Xj);
dudB = dudx*dxdB;%%
dxdG = FuncdxdG(K,dRXdG,RY,RZ,Xj);
dudG = dudx*dxdG;%%

dxdXj = FuncdxdXj(K,R);
dXjdPT = FuncdXjdPT(Phi,Theta);
dudPT = dudx*dxdXj*dXjdPT;%%

dudRi = cat(2,dudA,dudB,dudG);
