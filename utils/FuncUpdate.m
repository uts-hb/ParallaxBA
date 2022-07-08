
function [PVector] = FuncUpdate(PVector,DeltaP,DeltaF)

PVector.Pose = PVector.Pose+DeltaP;
PVector.Feature = PVector.Feature+DeltaF;


