function [phase_center] = lever_arm_example(param, tx_weights, rx_paths, ds)
% [phase_center] = lever_arm_example(param, tx_weights, rx_paths, ds)
%
% Example lever_arm function for +sim simulations. Returns antenna phase center
% locations for linear array. Centers of the transmit and receive arrays 
% are both at the origin.
%
% param: (not used)
% tx_weights: transmit amplitude weightings for each transmit phase center
%   These are amplitude weights, not power weights. Number of transmit
%   elements is assumed to be equal to the length of this array.
% rx_paths: Nc element array (only the length is used) where Nc is the
%   number of receive path phase centers
% ds: 3x1 vector describing separation between sensors (meters), used for
%   both the transmit and receive arrays. ds(1) is dx, ds(2) is dy, and
%   ds(3) is dz.
%
% phase_center = lever arm to each phase center specified by
%   tx_weights and rx_paths. 3 by Nc array where phase_center(1,:) is the
%   x-position, phase_center(2,:) is y, and phase_center(3,:) is z.
%
% See lever_arm.m for a full description of the coordinate system used.
%
% Author: John Paden, Theresa Stumpf

% Specify array positions
Nc = length(rx_paths);
physical_constants;
LArx(1,1:Nc) = ds(1) * (-(Nc-1)/2 : (Nc-1)/2);
LArx(2,1:Nc) = ds(2) * (-(Nc-1)/2 : (Nc-1)/2);
LArx(3,1:Nc) = ds(3) * (-(Nc-1)/2 : (Nc-1)/2);

Nc_tx = length(tx_weights);
LAtx(1,1:Nc_tx) = ds(1) * (-(Nc_tx-1)/2 : (Nc_tx-1)/2);
LAtx(2,1:Nc_tx) = ds(2) * (-(Nc_tx-1)/2 : (Nc_tx-1)/2);
LAtx(3,1:Nc_tx) = ds(3) * (-(Nc_tx-1)/2 : (Nc_tx-1)/2);

% =========================================================================
%% Compute Phase Centers
% =========================================================================

% Amplitude (not power) weightings for transmit side.
A = tx_weights;
magsum = sum(A);

% Weighted average of Xb, Yb and Zb components
LAtx_pc(1,1)    = dot(LAtx(1,:),A)/magsum;
LAtx_pc(2,1)    = dot(LAtx(2,:),A)/magsum;
LAtx_pc(3,1)    = dot(LAtx(3,:),A)/magsum;

phase_center = (LArx(:,rx_paths) + repmat(LAtx_pc,[1 numel(rx_paths)]))./2;

return

