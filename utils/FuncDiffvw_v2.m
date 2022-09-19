
function [Error, Sum_Error]= FuncDiffvw(Pose, camera, measurements)

Error = [];
R_CI = quatToRotMat(camera.q_CI);
t_CI = camera.p_C_I;

for i = 1 : size(Pose, 1)/6 -1

    A = R_CI' * v_RMatrixYPR22(Pose(6*(i-1)+1:6*(i-1)+3,1)');
    B = Pose(6*i+4:6*i+6, 1) - Pose(6*(i-1)+4:6*(i-1)+6, 1) -  (R_CI' * v_RMatrixYPR22(Pose(6*i+1:6*i+3,1)'))' * t_CI;
    fv = (A*B + t_CI) / measurements{1,i}.dT;

    C = R_CI' * v_RMatrixYPR22(Pose(6*i+1:6*i+3,1)');
    fw = -(C * A' - eye(3))/measurements{1,i}.dT;
    fw = crossMatToVec(fw);

    w_error = fw - measurements{1, i}.omega;
    v_error = fv -  measurements{1, i}.v;


    Error = [Error; v_error;w_error];


end
% Error 
Sum_Error = Error' * Error;








