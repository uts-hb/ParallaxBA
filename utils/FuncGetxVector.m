function xVector = FuncGetxVector(xVector,Image,PoseID)

[nRowNum,nColumNum] = size(Image);

PID(1:nRowNum-1,1) = PoseID;
FID(1:nRowNum-1,1) = Image(2:nRowNum,1);

for i=1:nRowNum-1;
    u(2*i-1:2*i,1) = Image(i+1,2:3)';
end;

xVector.u = cat(1,xVector.u,u);
xVector.PID = cat(1,xVector.PID,PID);
xVector.FID = cat(1,xVector.FID,FID);

