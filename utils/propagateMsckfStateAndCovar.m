function [msckfState_prop, phik] = propagateMsckfStateAndCovar(msckfState, measurements_k, noiseParams, phik, camera, state_k)

% Jacobians
Q_imu = noiseParams.Q_imu;
F = calcF(msckfState.imuState, measurements_k);
G = calcG(msckfState.imuState);

%Propagate State
msckfState_prop.imuState = propagateImuState(msckfState.imuState, measurements_k);
C_IG = quatToRotMat(msckfState_prop.imuState.q_IG);
    % Compute camera pose from current IMU pose
q_CG = quatLeftComp(camera.q_CI) * msckfState_prop.imuState.q_IG;
p_C_G = msckfState_prop.imuState.p_I_G + C_IG' * camera.p_C_I;

msckfState_prop.camState.q_CG = q_CG; 
msckfState_prop.camState.P_C_G = p_C_G;
msckfState_prop.state_k = state_k + 1; 


% State Transition Matrix
Phi = eye(size(F,1)) + F * measurements_k.dT; % Leutenegger 2013

% IMU-IMU Covariance
%     msckfState_prop.imuCovar = msckfState.imuCovar + ...
%                                 ( F * msckfState.imuCovar ...
%                                 + msckfState.imuCovar * F' ...
%                                 + G * Q_imu * G' ) ...
%                                         * measurements_k.dT;

msckfState_prop.imuCovar = Phi * msckfState.imuCovar * Phi' ...
    + G * Q_imu * G' * measurements_k.dT; % Leutenegger 2013

% Enforce PSD-ness
msckfState_prop.imuCovar = enforcePSD(msckfState_prop.imuCovar);

% Camera-Camera Covariance
msckfState_prop.camCovar = msckfState.camCovar;

% IMU-Camera Covariance
msckfState_prop.imuCamCovar = Phi * msckfState.imuCamCovar;
msckfState_prop.camStates = msckfState.camStates;


% Keystate Covariance

if isempty(msckfState.keyStates) || ~isempty(phik) 
    
    msckfState_prop.imuKeyCovar = msckfState.imuKeyCovar;
    msckfState_prop.keyCovar =  msckfState.keyCovar;
    msckfState_prop.camKeyCovar =  msckfState.camKeyCovar;
    msckfState_prop.keyStates = msckfState.keyStates;
    
end 

if ~isempty(msckfState.keyStates) && isempty(phik)
    
    msckfState_prop.imuKeyCovar = Phi * msckfState.imuKeyCovar;
    msckfState_prop.keyCovar =  msckfState.keyCovar;
    msckfState_prop.camKeyCovar =  msckfState.camKeyCovar;
    msckfState_prop.keyStates = msckfState.keyStates;
    
end

if ~isempty(phik) 
    
     phik(1:12, :) = Phi*phik(1:12, :);

end 


% %
%     P = [msckfState_prop.imuCovar, msckfState_prop.imuKeyCovar, msckfState_prop.imuCamCovar;
%         msckfState_prop.imuKeyCovar', msckfState_prop.keyCovar, msckfState_prop.camKeyCovar';
%         msckfState_prop.imuCamCovar', msckfState_prop.camKeyCovar, msckfState_prop.camCovar];
% %
%     chol(P)
%

%
%
%     % Keystate Covariance
%
end