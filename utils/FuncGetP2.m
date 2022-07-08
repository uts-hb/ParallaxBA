function P2 = FuncGetP2(PVector,FixVa)

%%
if FixVa==3;%% Fix Z
    PP1 = PVector.Pose(7:11)'*PVector.Pose(7:11);
    PP2 = PVector.Pose(13:end)'*PVector.Pose(13:end);
elseif FixVa==2;%% Fix Y
    PP1 = PVector.Pose(7:10)'*PVector.Pose(7:10);
    PP2 = PVector.Pose(12:end)'*PVector.Pose(12:end);
elseif FixVa==1;%% Fix X
    PP1 = PVector.Pose(7:9)'*PVector.Pose(7:9);
    PP2 = PVector.Pose(11:end)'*PVector.Pose(11:end);
end;
%%
PP3 = PVector.Feature'*PVector.Feature;
P2 = sqrt(PP1+PP2+PP3);