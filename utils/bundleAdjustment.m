function [xyzRefinedPoints, refinedPoses, reprojectionErrors,Errors,iter] = ...
    bundleAdjustment(xyzPoints, pointTracks, cameraPoses, intrinsics,Errors, varargin)
%bundleAdjustment Refine camera poses and 3-D points.
%   [xyzRefinedPoints, refinedPoses] = bundleAdjustment(xyzPoints,
%       pointTracks, cameraPoses, intrinsics) returns the refined 3-D
%   points and camera poses that minimize the reprojection errors. 3-D
%   points and camera poses are placed in the same world coordinate system. 
%   The refinement procedure is a variant of Levenberg-Marquardt algorithm.
%
%   Inputs:
%   -------
%   xyzPoints        - an M-by-3 matrix of 3-D [x, y, z] locations.
%
%   pointTracks      - an N-element array of pointTrack objects, where each
%                      element contains two or more points matched across
%                      multiple images.
%
%   cameraPoses      - a table containing two columns: 'ViewId' and 
%                      'AbsolutePose', typically produced by the poses
%                      method of imageviewset. The view IDs in pointTracks 
%                      refer to the view IDs in cameraPoses.
%
%   intrinsics       - a scalar or an M-element array of cameraIntrinsics
%                      objects, where M is the number of camera poses. Use 
%                      a scalar when images are captured using the same camera
%                      and a vector when images are captured by different cameras.
%
%   Outputs:
%   --------
%   xyzRefinedPoints - an M-by-3 matrix of refined point locations.
%
%   refinedPoses     - a table containing the refined camera poses.
%
%   [..., reprojectionErrors] = bundleAdjustment(...) additionally returns
%   reprojectionErrors, an N-element vector containing the mean
%   reprojection error for each 3-D world point.
%
%   [...] = bundleAdjustment(..., Name, Value) specifies additional
%   name-value pairs described below:
%
%   'MaxIterations'          A positive integer to specify the maximum
%                            number of iterations before Levenberg-Marquardt
%                            algorithm stops.
%
%                            Default: 50
%
%   'AbsoluteTolerance'      A positive scalar to specify termination
%                            tolerance of mean squared reprojection error
%                            in pixels.
%
%                            Default: 1
%
%   'RelativeTolerance'      A positive scalar to specify termination
%                            tolerance of relative reduction in
%                            reprojection error between iterations.
%
%                            Default: 1e-5
%
%   'PointsUndistorted'      A boolean to specify whether the 2-D points in
%                            pointTracks are from images without lens
%                            distortion or not.
%
%                            Default: False
%
%   'FixedViewIDs'           A vector of nonnegative integer to specify the
%                            reference cameras whose pose are fixed during
%                            optimization. This integer refers to the view
%                            IDs in cameraPoses. When it is empty, all
%                            camera poses are optimized.
%
%                            Default: []
%
%   'Verbose'                Set true to display progress information.
%
%                            Default: False
%
%   Class Support
%   -------------
%   xyzPoints must be single or double.
%
%   Example 1: Refine camera poses and 3-D points from a single camera
%   ------------------------------------------------------------------
%   % Load data for initialization
%   data = load('sfmGlobe');
%
%   % Refine the camera poses and points
%   [xyzRefinedPoints, refinedPoses] = ...
%       bundleAdjustment(data.xyzPoints, data.pointTracks, data.cameraPoses, data.intrinsics);
%
%   % Display the refined camera poses and 3-D world points
%   pcshow(xyzRefinedPoints, 'VerticalAxis', 'y', 'VerticalAxisDir', 'down', ...
%       'MarkerSize', 45);
%   hold on
%   plotCamera(refinedPoses, 'Size', 0.2);
%   hold off
%   grid on
%
%   Example 2: Structure from motion from multiple views
%   ----------------------------------------------------
%   % This example shows you how to estimate the poses of a calibrated
%   % camera from a sequence of views, and reconstruct the 3-D structure of
%   % the scene up to an unknown scale factor.
%   % <a href="matlab:helpview(fullfile(docroot,'toolbox','vision','vision.map'),'StructureFromMotionFromMultipleViewsExample')">View example</a>
%
% See also bundleAdjustmentStructure, bundleAdjustmentMotion, pointTrack, 
%   imageviewset, triangulateMultiview, cameraParameters, cameraIntrinsics.

% Copyright 2015-2019 The MathWorks, Inc.
%
% References
% ----------
% [1] M.I.A. Lourakis and A.A. Argyros (2009). "SBA: A Software Package for
%     Generic Sparse Bundle Adjustment". ACM Transactions on Mathematical
%     Software (ACM) 36 (1): 1-30.
%
% [2] R. Hartley, A. Zisserman, "Multiple View Geometry in Computer
%     Vision," Cambridge University Press, 2003.
%
% [3] B. Triggs; P. McLauchlan; R. Hartley; A. Fitzgibbon (1999). "Bundle
%     Adjustment: A Modern Synthesis". Proceedings of the International
%     Workshop on Vision Algorithms. Springer-Verlag. pp. 298-372.

[xyzRefinedPoints, refinedPoses, reprojectionErrors,Errors,iter] = ...
    sparseBA(xyzPoints, pointTracks, ...
    cameraPoses, intrinsics,Errors, 'full', mfilename, varargin{:});