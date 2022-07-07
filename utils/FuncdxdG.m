function dxdG = FuncdxdG(K,dRXdG,RY,RZ,Xj)

dxdG = K*dRXdG*RY*RZ*Xj;