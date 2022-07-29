function [measurements, visibility, cameraMatrices, quaternionBases, ...
    intrinsicsStruct, returnPoseType] = convertInputDataFormat(imagePoints, ...
    cameraPoses, intrinsics, isUndistorted, optimType)
%convertInputDataFormat Convert inputs to internal data structure.
%   All returned values are double.
%
%   measurements: 2-by-M packed 2-D points
%   visibility(i,j): a sparse matrix, true if point i is visible in view j
%   cameraMatrices: 6-by-V, rotation vector + translation vector
%   quaternionBases: 4-by-V, quaternions for initial rotations
%   intrinsicsStruct: a structure of camera intrinsics

% Copyright 2019 The MathWorks, Inc.

if strcmp(optimType, 'motion')
    numViews = 1;
    
    visibility   = sparse(ones(size(imagePoints,1), 1));
    measurements = double(imagePoints)';
    
    % Convert the intrinsics object to a simple structure
    % Note, the internal reprojection function uses a different definition
    % of skew factor, i.e., s = S / fc(1)
    intrinsicsStruct = struct('focalLength',double(intrinsics.FocalLength), ...
        'principalPoint', double(intrinsics.PrincipalPoint), ...
        'radialDistortion', double(intrinsics.RadialDistortion), ...
        'tangentialDistortion', double(intrinsics.TangentialDistortion), ...
        'skew', double(intrinsics.Skew / intrinsics.FocalLength(1)));

else % structure or full
    numViews  = height(cameraPoses);
    numPoints = numel(imagePoints);
    
    % visibility(i,j): true if point i is visible in view j
    visibility = zeros(numPoints, numViews);
    viewIds = cameraPoses.ViewId;
    x = zeros(numPoints, numViews);
    y = zeros(numPoints, numViews);
    for m = 1:numPoints
        trackViewIds = imagePoints(m).ViewIds;
        for n = 1:length(trackViewIds)
            viewIndex = find(viewIds == trackViewIds(n), 1, 'first');
            if isempty(viewIndex)
                error(message('vision:absolutePoses:missingViewId', trackViewIds(n)));
            end
            visibility(m, viewIndex) = 1;
            x(m, viewIndex) =  imagePoints(m).Points(n, 1);
            y(m, viewIndex) =  imagePoints(m).Points(n, 2);
        end
    end
    
    isVisible = find(visibility);
    x = x(isVisible);
    y = y(isVisible);
    
    visibility = sparse(visibility);
    
    % Measurements stores 2-D points in 1st view first, then 2nd view, ...
    measurements = double([x, y])';

    % Convert cameraParameters to cameraIntrinsics to allow parentheses-style indexing
    if isa(intrinsics, 'cameraParameters')
        intrinsics = convertToIntrinsics(intrinsics);
    end
    
    numCameras = numel(intrinsics);
    intrinsicsStruct(numCameras) = struct('focalLength',[], ...
        'principalPoint', [], ...
        'radialDistortion', [], ...
        'tangentialDistortion', [], ...
        'skew', []);
    
    for n = 1:numCameras
        % Convert the intrinsics object to a simple structure
        intrinsicsStruct(n).focalLength = double(intrinsics(n).FocalLength);
        intrinsicsStruct(n).principalPoint = double(intrinsics(n).PrincipalPoint);
        if ~isUndistorted
            % Skip if the distortion coefficients are all zeros
            if (any(intrinsics(n).RadialDistortion) || any(intrinsics(n).TangentialDistortion))
                intrinsicsStruct(n).radialDistortion = double(intrinsics(n).RadialDistortion);
                intrinsicsStruct(n).tangentialDistortion = double(intrinsics(n).TangentialDistortion);
            end
        end
        % Note, the internal reprojection function uses a different definition
        % of skew factor, i.e., s = S / fc(1)
        intrinsicsStruct(n).skew = double(intrinsics(n).Skew / intrinsics(n).FocalLength(1));
    end
end

[cameraMatrices, quaternionBases, returnPoseType] = convertToProjectionMatrices(cameraPoses, numViews);

end
%==========================================================================
% Convert cameraParameters to cameraIntrinsics
%==========================================================================
function intrinsics = convertToIntrinsics(camParam)
% Image size is not used in bundleAdjustment, but required in
% constructing the cameraIntrinsics object. When it is not available in
% the cameraParameters-type input, set it to a default value.
imageSize = camParam.ImageSize;
if isempty(imageSize)
    imageSize = [1 1];
end

intrinsics = cameraIntrinsics(camParam.FocalLength, camParam.PrincipalPoint, imageSize);
end

%==========================================================================
% Convert cameraPoses to a compact form
%==========================================================================
function [cameraMatrices, quaternionBases, returnPoseType] = convertToProjectionMatrices(cameraPoses, numViews)
% Convert camera poses to a compact form of camera projection matrices
% Use quaternion for numerical stability
if isa(cameraPoses, 'rigid3d') % In motion mode, cameraPoses is a rigid3d object
    t = double(cameraPoses.T(4,1:3));
    R = double(cameraPoses.T(1:3,1:3));
    cameraMatrices(4:6, 1) = -t*R';
    quaternionBases = vision.internal.quaternion.rotationToQuaternion(R);
    
    returnPoseType = class(t);
else
    hasAbsolutePose = ismember('AbsolutePose', cameraPoses.Properties.VariableNames);
   
    if hasAbsolutePose
        returnPoseType = class(cameraPoses.AbsolutePose(1).Translation);
    else
        returnPoseType = class(cameraPoses.Location{1});
    end
    
    quaternionBases = zeros(4, numViews); % double
    cameraMatrices  = zeros(6, numViews); % double
    
    for j = 1:numViews
        if hasAbsolutePose
            t = cameraPoses.AbsolutePose(j).Translation;
            R = cameraPoses.AbsolutePose(j).Rotation;
        else
            t = cameraPoses.Location{j};
            R = cameraPoses.Orientation{j};
        end
        
        cameraMatrices(4:6, j) = -t*R';
        quaternionBases(:, j) = vision.internal.quaternion.rotationToQuaternion(R);
    end
end
end