
function [uv] = FuncUV(Feature)

[nRowNumF,] = size(Feature);
uv = zeros(1);
k = 0;

for j=1:nRowNumF;
    for num=1:Feature(j,4);
        k=k+1;
        uv(k,1) = Feature(j,3*(num-1)+6);
        k=k+1;
        uv(k,1) = Feature(j,3*(num-1)+7);
    end;
end;


