
function [Error, Sum_Error,reprojectionErrors,Errors]= FuncDiff(xVector,PVector,Feature,K,Errors)

[uvcomp] = FuncfP(xVector,PVector,Feature,K);
Error = uvcomp-xVector.u;
% Errors = {}; 

% er= reshape(Error,2,[]);
% Errors{end+1} = er; 
Errors{end+1} = Error;

for i = 1 : length(Error)/2
    errors(:,i) = Error(2*(i-1)+1:2*i,1)';
end

curMeanError = errors(1,:).^2+errors(2,:).^2;
f_ID = []; 
for i = 1 : length(Feature)
    if Feature(i,1) ~= 0
        f_ID(end+1) = Feature(i,1);
    end
end

reprojectionErrors = zeros(size(f_ID,2),1); 

for i = 1 : length(f_ID) 
    A = find(xVector.FID==f_ID(i));
    e = sqrt(curMeanError(A)); 
    reprojectionErrors(i) = sum(e) / numel(e);
end 



    

%%
% Threshold = 5;
% for i=1:length(Error);
%     if Error(i,1) >= Threshold;
%         Error(i,1) = Threshold;
%     elseif Error(i,1) <= -Threshold;
%         Error(i,1) = -Threshold;
%     end;
% end;

%%
Sum_Error = Error'*Error;

