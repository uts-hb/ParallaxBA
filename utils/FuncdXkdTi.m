function dXkdTi = FuncdXkdTi(Xj,dSinOdTi,Dik,dDikdTi,SinO,Omega,dTikdTi)

dXkdTi = Xj*(dSinOdTi*Dik+dDikdTi*SinO)-sin(Omega)*dTikdTi;