function dCosO2dTk = FuncdCosO2dTk(dDotPdTk,DotP,dDikdTk,Dik)

dCosO2dTk = (dDotPdTk*Dik-dDikdTk*DotP)/Dik^2;