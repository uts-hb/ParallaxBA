function dxdA = FuncdxdA(K,RX,RY,dRZdA,Xj)

dxdA = K*RX*RY*dRZdA*Xj;