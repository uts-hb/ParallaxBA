function dxdB = FuncdxdB(K,RX,dRYdB,RZ,Xj)

dxdB = K*RX*dRYdB*RZ*Xj;