function [xVector,PVector,Feature,ID] = FuncDeletePts(xVector,PVector,Feature,ID)

[nRowNumO,nColumNum] = size(PVector.Omega);
[nRowNumu,nColumNum] = size(xVector.u);

m = 0;
a = 0;
FeaturePT = zeros(0);
Omega = zeros(0);
Feature2 = zeros(0);
ID2 = zeros(0);
DeletePts = zeros(0);
u = zeros(0);
PID = zeros(0);
FID = zeros(0);

for i = 1:nRowNumO;
    if PVector.Omega(i) > 0;
        m = m+1;
        Omega(m,1) = PVector.Omega(i,1);
        FeaturePT(2*m-1:2*m,1) = PVector.Feature(2*i-1:2*i,1);
        Feature2(m,:) = Feature(i,1:end);
        ID2(m,1) = ID(i,1);
    else
        a = a+1;
        DeletePts(a) = i;
    end;
end;

n = 0;
for j = 1:nRowNumu/2;
    for k = 1:a;
        flag = 0;
    	if xVector.FID(j,1) == DeletePts(k);
            flag = 1;
            break;
        end;
    end;
    if flag == 0;
        n = n+1;
        u(2*n-1:2*n,1) = xVector.u(2*j-1:2*j,1);
        PID(n,1) = xVector.PID(j,1);
        FID(n,1) = xVector.FID(j,1);
    end;
end;
            
clear PVector.Omega PVector.Feature ID Feature xVector;
PVector.Omega = Omega;
PVector.Feature = FeaturePT;
ID = ID2;
Feature = Feature2;
xVector.u = u;
xVector.PID = PID;
[nRowNumF,nColumNum] = size(FID);
for l = 1:nRowNumF;
    TT = 1;
    while(TT <= a);
        if FID(l,1) < DeletePts(TT);
            break;
        else
            TT = TT+1;
        end;
    end;
    xVector.FID(l,1) = FID(l,1)-TT+1;
end;
% xVector.FID = FID;
% xVector.FID = [1:nRowNumu/2-a,1]';
clear Omega FeaturePT ID2 Feature2 u PID FID;
        
