function dXkdTk = FuncdXldTk(Xj,dSinOdTk,Dik,dDikdTk,SinO)

dXkdTk = Xj*(dSinOdTk*Dik+dDikdTk*SinO);