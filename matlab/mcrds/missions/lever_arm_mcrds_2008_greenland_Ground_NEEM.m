function [phase_center] = lever_arm_mcords_2009_antarctica_TO(tx_weights, rxchannel)
% [phase_center] = lever_arm_mcrds_2008_greenland_Ground_NEEM(tx_weights, rxchannel)
%
% Returns lever arm position for 2008 Greenland Ground-based NEEM data.

% tx_weights = transmit amplitude weightings (from the radar_config file)
%   These are amplitude weights, not power weights.
% rxchannel = receive channel to return phase_center for (scalar,
%   positive integer). If left empty or undefined, every phase center
%   is returned.  Set to 0 to return the middle of the array.
%
% phase_center = lever arm to phase center of measurement defined by
%   tx_weights and rxchannel
%
% =========================================================================
% REMARKS:
% 
% 1). Lever arm refers to a (3 x 1) vector that expresses the position of 
%     each phase center relative to the INS unit in the coordinate
%     system of the plane's body (Xb, Yb, Zb).  This is a righthanded,
%     orthoganol system that agrees with aerospace convention.  +Xb points
%     from the plane's center of gravity towards its nose.  +Yb points from
%     the plane's center of gravity along the right wing.  +Zb points from
%     the plane's center of gravity down towards the Earth's surface.
% 
% 2). The lever arm of the Nth receive channel is defined using the 
%     following syntax:
%
%     LArx_N = [Xb_N; Yb_N; Zb_N]
% 
% 3). Values were determined using the following assumptions:
%     (i).  The origin of the (Xb, Yb, Zb) coordinate system is the INS
%           unit, which is assumed to be a point on the floor of the Twin
%           Otter fuselage at 225 inches aft and lying on the plane's 
%           centerline. 
% ========================================================================

if ~exist('rxchannel','var') || isempty(rxchannel)
  rxchannel = 1:4;
end

% absoulute value of components

gps.x = 0*0.0254;
gps.y = 0*0.0254;
gps.z = 0*0.0254;

LArx(1,:)   = -[ 4.365 4.365 4.365 4.365 4.365 4.365 4.365 4.365]/2 - gps.x;  % m
LArx(2,:)   = [ -3.5 -2.5 -1.5 -0.5 0.5 1.5 2.5 3.5]*0.857 - gps.y; % m
LArx(3,:)   = [ 0 0 0 0 0 0 0 0 ] - gps.z; % m

LAtx(1,:)   =  [ 4.365 4.365 ]/2 - gps.x; % m
LAtx(2,:)   =  [ -3.658 3.658 ]/2 - gps.y; % m
LAtx(3,:)   =  [ 0 0] - gps.z; % m

% Amplitude (not power) weightings for transmit side. 
if rxchannel == 0
  rxchannel = 3;
  tx_weights = ones(1,size(LAtx,2));
end
A = tx_weights;
magsum       = sum(A);

% Weighted average of Xb, Yb and Zb components
LAtx_pc(1,1)    = dot(LAtx(1,:),A)/magsum;
LAtx_pc(2,1)    = dot(LAtx(2,:),A)/magsum;
LAtx_pc(3,1)    = dot(LAtx(3,:),A)/magsum;

phase_center = (LArx(:,rxchannel) + repmat(LAtx_pc,[1 size(rxchannel,2)]))./2;

return