function p_f_G = BA_triangulation(camStates, observations, indices)
%CALCGNPOSEST Calculate the position estimate of the feature using Gauss
%Newton optimization
%   INPUT:
%   observations: 2xM matrix of pixel values of the current landmark
%   camStates: Cell array of M structs of camera poses
%   camera: intrinsic calibration
%   OUTPUT:
%   p_f_G: 3x1 feature vector in the global frame

%K is not needed if we assume observations are not pixels but x' = (u -
%c_u)/f_u

%K = [camera.f_u 0 camera.c_u; 0 camera.f_v camera.c_v; 0 0 1];

%Get initial estimate through intersection
%Use the first 2 camStates
secondViewIdx = length(indices);

C_12 = quatToRotMat(camStates{indices(:,1)}.q_C_G)*quatToRotMat(camStates{indices(:,secondViewIdx)}.q_C_G)';
t_21_1 = quatToRotMat(camStates{indices(:,1)}.q_C_G)*(camStates{indices(:,secondViewIdx)}.p_C_G - camStates{indices(:,1)}.p_C_G);

p_f1_1_bar = triangulate(observations(:,indices(:,1)), observations(:,indices(:,secondViewIdx)),C_12, t_21_1);

%initialEst = quatToRotMat(camStates{1}.q_CG)'*p_f1_1_bar + camStates{1}.p_C_G;


xBar = p_f1_1_bar(1);
yBar = p_f1_1_bar(2);
zBar = p_f1_1_bar(3);

alphaBar = xBar/zBar;
betaBar = yBar/zBar;
rhoBar = 1/zBar;

xEst = [alphaBar; betaBar; rhoBar];

p_f_G = (1/xEst(3))*quatToRotMat(camStates{indices(:,1)}.q_C_G)'*[xEst(1:2); 1] + camStates{indices(:,1)}.p_C_G; 


    function [p_f1_1] = triangulate(obs1, obs2, C_12, t_21_1)
        % triangulate Triangulates 3D points from two sets of feature vectors and a
        % a frame-to-frame transformation

           %Calculate unit vectors 
           v_1 = [obs1;1];
           v_2 = [obs2;1];
           v_1 = v_1/norm(v_1);
           v_2 = v_2/norm(v_2);

           A = [v_1 -C_12*v_2];
           b = t_21_1;

           scalar_consts = A\b;
           p_f1_1 = scalar_consts(1)*v_1;
    end

end

