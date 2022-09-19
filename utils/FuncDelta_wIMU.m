
function [DeltaP,DeltaF,Sum_Delta] = FuncDelta_wIMU(Jacobian,Error,PVector,Feature,W)

Info = Jacobian'*W*Jacobian;
E = -Jacobian'*W*Error;


nRowNumP = length(PVector.Pose);
nRowNumI = nRowNumP+length(PVector.Feature);

U = Info(1:nRowNumP,1:nRowNumP);
V = Info(nRowNumP+1:nRowNumI,nRowNumP+1:nRowNumI);
W = Info(1:nRowNumP,nRowNumP+1:nRowNumI);
EP = E(1:nRowNumP);
EF = E(nRowNumP+1:nRowNumI);

[ID1,ID2,Val] = find(V);

i = 1;
nF = 1;

while i<=length(ID1);
    if Feature(nF)==0;
        nF = nF+1;
    elseif Feature(nF)==3;
        aa = reshape(Val(i:i+8),3,3);
        aa = inv(aa);
        Val(i:i+8) = reshape(aa,9,1);
        nF = nF+1;
        i = i+9;
    elseif Feature(nF)==2;
        aa = reshape(Val(i:i+3),2,2);
        aa = inv(aa);
        Val(i:i+3) = reshape(aa,4,1);
        nF = nF+1;
        i = i+4;
    end;
end;

V = sparse(ID1,ID2,Val);

clear ID1 ID2 Val;

S = U-W*V*W';
ES = EP-W*V*EF;


S = S(7:nRowNumP,7:nRowNumP);
ES = [ES(7:nRowNumP)];

DeltaP = S\ES;



DeltaP = [zeros(6,1);DeltaP(1:nRowNumP-6)];    % First Z = 0.8770  4th Z = 0.7970

DeltaF = V*(EF-W'*DeltaP);

Sum_Delta = DeltaP'*DeltaP+DeltaF'*DeltaF;


