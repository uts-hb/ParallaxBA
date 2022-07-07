function dXkdPT = FuncdXkdPT(Dik,Xj,dSinOdPT,SinO,dXjdPT)

dXkdPT = Dik*(Xj*dSinOdPT+SinO*dXjdPT);