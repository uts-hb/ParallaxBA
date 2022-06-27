function dXkdTk = FuncdXkdTk(Xj,dSinOdTk,Dik,dDikdTk,SinO,Omega,dTikdTk)

dXkdTk = Xj*(dSinOdTk*Dik+dDikdTk*SinO)-sin(Omega)*dTikdTk;