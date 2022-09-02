
function [DeltaP,DeltaF,Sum_Delta] = FuncDelta(Jacobian,Error,PVector,Feature,FixVa)

Info = Jacobian'*Jacobian;
E = -Jacobian'*Error;


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

%%
if FixVa==3;%% Fix Z
    S = [S(7:11,7:11),S(7:11,13:nRowNumP);
         S(13:nRowNumP,7:11),S(13:nRowNumP,13:nRowNumP)];
    ES = [ES(7:11);ES(13:nRowNumP)];    
elseif FixVa==2;%% Fix Y
    S = [S(7:10,7:10),S(7:10,12:nRowNumP);
         S(12:nRowNumP,7:10),S(12:nRowNumP,12:nRowNumP)];
    ES = [ES(7:10);ES(12:nRowNumP)];
elseif FixVa==1;%% Fix X
    S = [S(7:9,7:9),S(7:9,11:nRowNumP);
         S(11:nRowNumP,7:9),S(11:nRowNumP,11:nRowNumP)];
    ES = [ES(7:9);ES(11:nRowNumP)];  
elseif FixVa==4 
%     S = [S(7:149,7:149),S(7:149,151:nRowNumP);
%         S(151:nRowNumP,7:149),S(151:nRowNumP,151:nRowNumP)];
%     ES = [ES(7:149);ES(151:nRowNumP)];
    S = [S(7:59,7:59),S(7:59,61:nRowNumP);
        S(61:nRowNumP,7:59),S(61:nRowNumP,61:nRowNumP)];
    ES = [ES(7:59);ES(61:nRowNumP)];
%         S = [S(7:29,7:29),S(7:29,31:nRowNumP);
%         S(31:nRowNumP,7:29),S(31:nRowNumP,31:nRowNumP)];
%     ES = [ES(7:29);ES(31:nRowNumP)];
elseif FixVa==5 
    S = [S(7:29,7:29),S(7:29,31:nRowNumP);
        S(31:nRowNumP,7:29),S(31:nRowNumP,31:nRowNumP)];
    ES = [ES(7:29);ES(31:nRowNumP)];
end;
%%
DeltaP = S\ES;

%%
if FixVa==3;%% Fix Z
    DeltaP = [zeros(6,1);DeltaP(1:5);0;DeltaP(6:nRowNumP-7)];   
elseif FixVa==2;%% Fix Y
    DeltaP = [zeros(6,1);DeltaP(1:4);0;DeltaP(5:nRowNumP-7)];
elseif FixVa==1;%% Fix X
    DeltaP = [zeros(6,1);DeltaP(1:3);0;DeltaP(4:nRowNumP-7)];
elseif FixVa==4; 
%     DeltaP = [zeros(6,1);DeltaP(1:142);0;DeltaP(143:nRowNumP-7)];   %
%     First Z = 0.8770 Twenty-five Z = 0.4559
    DeltaP = [zeros(6,1);DeltaP(1:53);0;DeltaP(54:nRowNumP-7)];    % First Z = 0.8770  4th Z = 0.7970
%         DeltaP = [zeros(6,1);DeltaP(1:22);0;DeltaP(23:nRowNumP-7)];    % First Z = 0.8770  4th Z = 0.7970
elseif FixVa==5; 
        DeltaP = [zeros(6,1);DeltaP(1:23);0;DeltaP(24:nRowNumP-7)];    % First Z = 0.8770  4th Z = 0.7970
end;
%%
DeltaF = V*(EF-W'*DeltaP);

Sum_Delta = DeltaP'*DeltaP+DeltaF'*DeltaF;


