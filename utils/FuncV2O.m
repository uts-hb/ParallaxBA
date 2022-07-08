function Omega = FuncV2O(V1,V2)

Omega = acos(dot(V1,V2)/(sqrt(V1(1)^2+V1(2)^2+V1(3)^2)*sqrt(V2(1)^2+V2(2)^2+V2(3)^2)));