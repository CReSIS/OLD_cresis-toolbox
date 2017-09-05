%Extract 3D surface of ice-bed layers.
%
% Compile with mex -v -largeArrayDims extract.cpp
%    -v: verbose
%    -largeArrayDims: required for 64 bit
% 
% Usage: surf = extract(input, surface, bottom, extra, mask, mean, variance, smooth_weight, smooth_var, smooth_slope)
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
% smooth_weight: optional (set to -1 for default), default weight is set 
%   in Instances.h
% smooth_var: optional (set to -1 for default), default weight is set in
%   Instances.h
% smooth_slope: Must be Ndoa-1 in length. This vector indicates the surface
%  slope that results in lowest cost. Element k of this vector corresponds
%  to the slope (change in rows) from column k to column k+1. For a surface
%  that is expected to be flat, this vector would be all zeros. For a
%  surface which monotonically increases from the left to right, the values
%  would be all positive. For a surface that is expected to increase 
%  and then decrease, the values would start positive and then become negative.
% 
% Authors:  Mingze Xu, John Paden
