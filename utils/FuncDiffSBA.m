
function [Error, Sum_Error]= FuncDiffSBA(xVector,PVector,Feature,K)

[uvcomp] = FuncfP(xVector,PVector,Feature,K);
Error = xVector.u-uvcomp;
Sum_Error = Error'*Error;

