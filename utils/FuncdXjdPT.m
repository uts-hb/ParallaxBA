function dXjdPT = FuncdXjdPT(Phi,Theta)

dXjdPT = [cos(Phi)*cos(Theta),-sin(Phi)*sin(Theta);
          0,cos(Theta);
          -sin(Phi)*cos(Theta),-cos(Phi)*sin(Theta)];