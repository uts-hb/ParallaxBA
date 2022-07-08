function dCosO2dTi = FuncdCosO2dTi(dDotPdTi,DotP,dDikdTi,Dik)

dCosO2dTi = (dDotPdTi*Dik-dDikdTi*DotP)/Dik^2;