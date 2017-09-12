%Extract 3D surface of ice-bed layers.
%
% Compile with mex -v -largeArrayDims refine.cpp
%    -v: verbose
%    -largeArrayDims: required for 64 bit
% 
% Usage: surf = extract(image, sgt, bgt, egt, mask, mean, variance, smooth_weight, smooth_var, [smooth_slope], [edge], [max_loops], [bounds])
%
% No input or output arguments may be missing unless indicated as optional
%
% All ground truth (surface, bottom, and extra) must have values which are
% less than Nt-length(mu) or they will have no effect.
%
% input: double, 3D image (Nt by Ndoa by Nx)
% sgt: double, 2D surface, Ndoa by Nx matrix, contains surface ground
%   truth
% bgt: double, Nx element matrix, contains bottom ground truth
%   DEPRECATED: WILL BE REPLACED WITH "extra" ONLY
% egt: double, 3 by Ne matrix, contains extra ground truth
%   row 1: Nx dimension index
%   row 2: Ndoa dimension index
%   row 3: Ny dimension index
% mask: double, 2D surface mask, Ndoa by Nx matrix, contains surface mask
%   A value of zero means no ice and the surface and bottom are forced to
%   be the same.
% mean: double, Nmean size vector, image magnitude template in Nt dimension
% variance: double, Nmean size vector, image magnitude weights
%   corresponding to mean
% smooth_weight: optional (set to -1 for default), default weight is set 
%   in Instances.h
% smooth_var: optional (set to -1 for default), default weight is set in
%   Instances.h
% smooth_slope: OPTIONAL: All zeroes is default. If passed in and not empty,
%  it must be Ndoa-1 in length. This vector indicates the surface slope
%  relative to the ice surface that results in lowest cost. Element k of
%  this vector corresponds to the slope (change in rows) from column k to
%  column k+1. For a surface that is expected to be flat relative to the
%  ice surface, this vector would be all zeros. For a surface which
%  monotonically increases from the left to right, the values would be all
%  positive. For a surface that is expected to increase  and then decrease,
%  the values would start positive and then become negative.
% edge: OPTIONAL: No edge constraint is the default. If passed in and not
%  empty, it must be an Ndoa by 2 matrix which provides the edge
%  conditions for the output surf. edge(:,1) specifies the surface at
%  y = 1 and edge(:,2) specifies the surface at y = 2.
% max_loops: OPTIONAL: 50 is the default. If passed in and not empty, it
%  must be a scalar. The routine then computes this many loops before
%  exiting.
% bounds: OPTIONAL: [0 Ndoa-1 0 Nx-1] is the default. If passed in and
%  not empty, it must be a 4 element double vector. It should contain
%  the start and stop DOA indexes followed by the start and stop X
%  dimension indexes. All indexes are zero-indexed. This allows a subset of
%  the data to be surface tracked by this program. Setting any element to
%  a negative value will cause the default to be used for that particular
%  element.
%
% surf: double, 2D surface, Ndoa by Nx matrix, contains ice bottom surface
%
% Authors:  Mingze Xu, John Paden
