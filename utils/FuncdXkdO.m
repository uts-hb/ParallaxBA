function dXkdO = FuncdXkdO(dSinOdO,Dik,Xi,Omega,Tik)

dXkdO = dSinOdO*Dik*Xi-cos(Omega)*Tik;