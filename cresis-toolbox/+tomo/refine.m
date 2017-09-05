%Extract 3D surface of ice-bed layers.
%
% Compile with mex -v -largeArrayDims refine.cpp
%    -v: verbose
%    -largeArrayDims: required for 64 bit
% 
% Usage: surf = refine(input, surface, bottom, extra, mask, mean, variance)
%
% No input or output arguments may be missing.
%
% input: double, 3D image (Nt by Ndoa by Nx)
% surface: double, 2D surface, Ndoa by Nx matrix, contains surface ground
%   truth
% bottom: double, Nx element matrix, contains bottom ground truth
%   DEPRECATED: WILL BE REPLACED WITH "extra" ONLY
% extra: double, 3 by Ne matrix, contains bottom ground truth
%   row 1: Nx dimension index
%   row 2: Ndoa dimension index
%   row 3: Ny dimension index
% mask: double, 2D surface mask, Ndoa by Nx matrix, contains surface mask
% mean: double, Nmean size vector, image magnitude template in Nt dimension
% variance: double, Nmean size vector, image magnitude weights
%   corresponding to mean
% 
% Authors:  Mingze Xu
