function [phase_center] = lever_arm_example(param, tx_weights, rx_paths, Nc, dy)
% [phase_center] = lever_arm_example(param, tx_weights, rx_paths, Nc, dy)
%
% Example lever_arm function for simulation. Returns antenna phase center
% locations for array.
%
% param = (not used)
% tx_weights = transmit amplitude weightings for each transmit phase center
%   These are amplitude weights, not power weights.
% rx_paths = receiver phase center(s) to return the phase center for
%   (vector of positive integers that specify the rx_path)
%
% phase_center = lever arm to each phase center specified by
%   tx_weights and rx_paths
%
% See lever_arm.m for a full description of the coordinate system used.
%
% Author: John Paden, Theresa Stumpf

% Specify array positions
N = length(rx_paths);
physical_constants;
LArx(1,1:N) = 0;
LArx(2,1:N) = Nc * dy * (-(N-1)/2:(N-1)/2);
LArx(3,1:N) = 0;

LAtx = LArx(:,1:N);

% =========================================================================
%% Compute Phase Centers
% =========================================================================

% Amplitude (not power) weightings for transmit side.
A = tx_weights;
magsum       = sum(A);

% Weighted average of Xb, Yb and Zb components
LAtx_pc(1,1)    = dot(LAtx(1,:),A)/magsum;
LAtx_pc(2,1)    = dot(LAtx(2,:),A)/magsum;
LAtx_pc(3,1)    = dot(LAtx(3,:),A)/magsum;

phase_center = (LArx(:,rx_paths) + repmat(LAtx_pc,[1 numel(rx_paths)]))./2;

return

