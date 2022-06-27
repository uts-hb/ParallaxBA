function [Phi,Theta] = FuncV2PT(V)

Phi = atan2(V(1),V(3));
Theta = atan2(V(2),sqrt(V(1)^2+V(3)^2));
