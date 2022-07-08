
function RMatrix = RMatrixYPR22(Alpha,Beta,Gamma)

Rx = [1 0 0;0 cos(Gamma) sin(Gamma);0 -sin(Gamma) cos(Gamma)];
Ry = [cos(Beta) 0 -sin(Beta);0 1 0;sin(Beta) 0 cos(Beta)];
Rz = [cos(Alpha) sin(Alpha) 0;-sin(Alpha) cos(Alpha) 0;0 0 1];
RMatrix = Rx*Ry*Rz;

end
