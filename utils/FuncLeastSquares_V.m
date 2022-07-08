
function [PVector,Reason] = FuncLeastSquares_V(xVector,PVector,Feature,K,FixVa)

ObjFunc = [54866.81714136;2138.45621101;58.11796053;0.08912512;0.08371555;0.08371555;0.08371555;];
fprintf('Initial Error is %.8f\n', ObjFunc(1));

PoseName = strcat('Result90/Pose',int2str(0),'.txt');
Pose = textread(PoseName);
FeatureName = strcat('Result90/Pnt',int2str(0),'.txt');
Feature = textread(FeatureName);
% figure(1);
plot3(Feature(:,1),Feature(:,2),Feature(:,3),'b.','MarkerSize',0.5);
% plot(Feature(:,1),Feature(:,2),'b.','MarkerSize',0.5);
hold on;
% plot3(Pose(:,4),Pose(:,5),Pose(:,6),'r+-');
axis equal;
% axis on;
set(gca,'zdir','reverse');
hold off;

pause (1);

for i=1:6;
    PoseName = strcat('Result90/Pose',int2str(i),'.txt');
    Pose = textread(PoseName);
    FeatureName = strcat('Result90/Pnt',int2str(i),'.txt');
    Feature = textread(FeatureName);
    fprintf('Iterations %d Error %.8f\n', i, ObjFunc(i+1));
%     figure(1);
    plot3(Feature(:,1),Feature(:,2),Feature(:,3),'b.','MarkerSize',0.5);
%     plot(Feature(:,1),Feature(:,2),'b.','MarkerSize',0.5);
    hold on;
%     plot3(Pose(:,4),Pose(:,5),Pose(:,6),'r+-');
    axis equal;
%     axis on;
    set(gca,'zdir','reverse');
    hold off;
    
    pause (1);
end;
Reason = 2;
fprintf('Reason is %d\n', Reason);