
function [Alpha,Beta,Gamma] = InvRotMatrixYPR22(K)

 Beta=atan2(-(K(1,3)),sqrt(K(1,1)^2+K(1,2)^2));
 if (cos(Beta)==0)
     Alpha=0;
     Beta=pi/2;
     Gamma=atan2(K(1,2),K(2,2));
 else
     Alpha=atan2(K(1,2)/cos(Beta),K(1,1)/cos(Beta));
     Gamma=atan2(K(2,3)/cos(Beta),K(3,3)/cos(Beta));
 end

 R = [Alpha,Beta,Gamma];

end
