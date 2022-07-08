
function [Delta,Sum_Delta] = FuncDeltaCG(Jacobian,Error,Iter)

% [U,D,V] = svd(Jacobin);
% 
% [nRowNumD,nColumNumD] = size(D);
% for i=1:nColumNumD;
%     if D(i,i) == 0;
%         D(i,i) = 0;
%     else
%         D(i,i) = 1/D(i,i);
%     end;
% end;
% % D2 = D(1:nColumNumD,1:nColumNumD);
% % D(1:nColumNumD,1:nColumNumD) = inv(D2);
% 
% Delta = -V*D'*U'*Error;
% Sum_Delta = sum(abs(Delta));

% Delta = (Jacobian'*Jacobian)\(-Jacobian'*Error);
% Delta = bicg(Jacobian'*Jacobian,-Jacobian'*Error);
% Delta = bicgstab(Jacobian'*Jacobian,-Jacobian'*Error,10e-6,888+800*Iter);
% Delta = bicgstab(Jacobian'*Jacobian,-Jacobian'*Error,1e-6,222);
Delta = bicgstab(Jacobian'*Jacobian,-Jacobian'*Error);
Sum_Delta = Delta'*Delta;
