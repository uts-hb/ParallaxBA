function V = Funcuv2V(uv,R,K)

RMatrix = RMatrixYPR22(R(1),R(2),R(3));
x = [uv(1);uv(2);1];
V = K*RMatrix\x;
