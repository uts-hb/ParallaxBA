function [xyzPoints, intrinsics, maxIterations, absTol, relTol, isUndistorted, ...
    verbose, returnType, pointTracks, cameraPoses, fixedCameraIndex] = ...
    validateAndParseInputs(xyzPoints, pointTracks, cameraPoses, intrinsics, ...
    optimType, funcName, varargin)
%validateAndParseInputs Parameter validation and parsing for bundle
%   adjustment functions.

% Copyright 2019 The MathWorks, Inc.

% Set input parser
defaults = struct(...
    'MaxIterations', 50,...
    'AbsoluteTolerance', 1,...
    'RelativeTolerance', 1e-5,...
    'PointsUndistorted', false, ...
    'FixedViewIDs', [], ...
    'Verbose', false);

% xyzPoints
[xyzPoints, returnType] = validateXYZPoints(xyzPoints, funcName); 

isFull      = strncmpi(optimType, 'full', 1);
isMotion    = strncmpi(optimType, 'motion', 1);
isStructure = strncmpi(optimType, 'structure', 1);

parser = inputParser;
parser.CaseSensitive = false;
parser.FunctionName = funcName;

if isMotion
    % imagePoints
    pointTracks = validateImagePoints(pointTracks, xyzPoints, funcName);
    
    % absolutePose
    validateSinglePose(cameraPoses, funcName);
    
    % intrinsics
    validateCameraIntrinsics(intrinsics, cameraPoses, true, funcName);
    
    fixedCameraIndex = 0;
elseif isFull
    % pointTracks
    validatePointTracks(pointTracks, xyzPoints, funcName);
    
    % cameraPoses
    validatePoseTableFull(cameraPoses, funcName);
    
    % intrinsics
    validateCameraIntrinsics(intrinsics, cameraPoses, false, funcName);  

    % Parameter specific to 'full' mode
    parser.addParameter('FixedViewIDs', defaults.FixedViewIDs, ...
        @(x)validateFixedViewIDs(x, funcName));
else % structure
    % pointTracks
    validatePointTracks(pointTracks, xyzPoints, funcName);
    
    % cameraPoses
    validatePoseTableStructure(cameraPoses, funcName);

    % intrinsics
    validateCameraIntrinsics(intrinsics, cameraPoses, false, funcName);
end

% Common parameters
parser.addParameter('MaxIterations', defaults.MaxIterations, ...
    @(x)validateMaxIterations(x, funcName));
parser.addParameter('AbsoluteTolerance', defaults.AbsoluteTolerance, ...
    @(x)validateTolerance(x, funcName, 'AbsoluteTolerance'));
parser.addParameter('RelativeTolerance', defaults.RelativeTolerance, ...
    @(x)validateTolerance(x, funcName, 'RelativeTolerance'));
parser.addParameter('PointsUndistorted', defaults.PointsUndistorted, ...
    @(x)vision.internal.inputValidation.validateLogical(x, 'PointsUndistorted'));
parser.addParameter('Verbose', defaults.Verbose, ...
    @(x)vision.internal.inputValidation.validateLogical(x, 'Verbose'));

parser.parse(varargin{:});

maxIterations = parser.Results.MaxIterations;
absTol        = double(parser.Results.AbsoluteTolerance);
relTol        = double(parser.Results.RelativeTolerance);
isUndistorted = parser.Results.PointsUndistorted;
verbose       = parser.Results.Verbose;

% Convert FixedViewIDs to indices in the camera poses table
if isFull
    fixedViewIDs  = uint32(parser.Results.FixedViewIDs);
    
    % Check the fixed view ID
    if ~isempty(fixedViewIDs)
        fixedCameraIndex = convertFixedViewIDsToIndices(cameraPoses.ViewId, fixedViewIDs);
    else
        fixedCameraIndex = 0;
    end
    
elseif isStructure
    % All cameras are fixed
    fixedCameraIndex = 1:height(cameraPoses);
end

end

%--------------------------------------------------------------------------
function fixedCameraIndex = convertFixedViewIDsToIndices(viewIds, fixedViewIDs)
if ~isempty(fixedViewIDs)
    % Check if fixedViewIDs are valid
    missingViewIdx = ~ismember(fixedViewIDs, viewIds);
    if any(missingViewIdx)
        missingViewIds = fixedViewIDs(missingViewIdx);
        error(message('vision:viewSet:missingViewId', ...
            missingViewIds(1)));
    end
    [~, fixedCameraIndex] = intersect(viewIds, fixedViewIDs, 'stable');
else
    fixedCameraIndex = 0;
end
end

%==========================================================================
% Validate Requried Inputs
%==========================================================================
function [xyzPoints, returnType] = validateXYZPoints(xyzPoints, funcName)
validateattributes(xyzPoints, {'single', 'double'}, ...
    {'finite', 'nonempty', 'nonsparse', '2d', 'ncols', 3}, funcName, 'xyzPoints');
returnType = class(xyzPoints); 
xyzPoints  = double(xyzPoints');
end

%--------------------------------------------------------------------------
function validatePointTracks(points, xyzPoints, funcName)

validateattributes(points, {'pointTrack'}, {'nonempty','vector'}, funcName, 'pointTracks');

% Check the size of input
if numel(points) ~= size(xyzPoints, 2)
    error(message('vision:sfm:unmatchedXYZTrack'));
end
end

%--------------------------------------------------------------------------
function points = validateImagePoints(points, xyzPoints, funcName)
points = vision.internal.inputValidation.checkAndConvertPoints(points, funcName, 'imagePoints');

% Check the size of input
if size(points, 1) ~= size(xyzPoints, 2)
    error(message('vision:sfm:unmatchedXYZImagePoints'));
end
end

%--------------------------------------------------------------------------
function validateSinglePose(pose, funcName)
validateattributes(pose, {'rigid3d'},{'scalar','nonempty'}, funcName, 'absolutePose');
end

%--------------------------------------------------------------------------
function validatePoseTableFull(poses, funcName)

% Validate the table
validateattributes(poses, {'table'},{'nonempty'}, funcName, 'cameraPoses');

% Check columns. In 'full' mode, the camera poses table can contain
% 'AbsolutePose' or 'Location' and 'Orientation'
if ismember('AbsolutePose', poses.Properties.VariableNames)
    vision.internal.inputValidation.validatePoseTableRigid3d(poses, funcName, 'cameraPoses');
else
    vision.internal.inputValidation.checkAbsolutePoses(poses, funcName, 'cameraPoses');
end
end

%--------------------------------------------------------------------------
function validatePoseTableStructure(poses, funcName)
vision.internal.inputValidation.validatePoseTableRigid3d(poses, funcName, 'cameraPoses');
end

%--------------------------------------------------------------------------
function validateCameraIntrinsics(intrinsics, cameraPoses, isScalar, funcName)
% cameraParameters is supported but not recommended or documented and it
% can only be a scalar due to restriction in its constructor.
vision.internal.inputValidation.checkIntrinsicsAndParameters( ...
    intrinsics, isScalar, funcName);

% Check the camera array
if ~isscalar(intrinsics) && numel(intrinsics) ~= height(cameraPoses)
    error(message('vision:sfm:unmatchedParamsPoses'));
end
end

%==========================================================================
% Validate Optional Parameters
%==========================================================================
function tf = validateFixedViewIDs(value, funcName)
tf = true;
if ~isempty(value)
    validateattributes(value,{'numeric'}, {'vector','integer','nonnegative'}, ...
        funcName, 'FixedViewIDs');
end
end

%--------------------------------------------------------------------------
function validateMaxIterations(maxIter, funcName)
validateattributes(maxIter,{'single', 'double'}, {'scalar','integer', 'positive'}, ...
    funcName, 'MaxIterations')
end

%--------------------------------------------------------------------------
function validateTolerance(tol, funcName, argName)
validateattributes(tol,{'single', 'double'}, {'real','nonnegative','scalar','finite'}, ...
    funcName, argName);
end