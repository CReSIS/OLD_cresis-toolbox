function [phase_center] = lever_arm_mcords_simulator(tx_weights, rxchannel)
% [phase_center] = lever_arm_mcords_simulator(tx_weights, rxchannel)
%
% Returns lever arm position for CSARP simulator.
%
% tx_weights = transmit amplitude weightings (from the radar_config file)
%   These are amplitude weights, not power weights.
% rxchannel = receive channel to return phase_center for (scalar,
%   positive integer)
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

LArx(1,:)   = -[ 0 0 0 0 ] - gps.x; % m
% LArx(2,:)   = -[ 0 1.5374/2 1.5374 3/2*1.5374] - gps.y; % m
LArx(2,:)   = -[ 0 1.9986/2 1.9986 3/2*1.9986] - gps.y; % m
LArx(3,:)   = -[ 0 0.05 0.1 0.15] - gps.z; % m
% LArx(3,:)   = -[ 0 0 0 0.] - gps.z; % m

LAtx(1,:)   = -[ 0 0 0 0 ] - gps.x; % m
% LAtx(2,:)   =  [ 0 1.5374/2 1.5374 3/2*1.5374] - gps.y; % m
LAtx(2,:)   =  [ 0 1.9986/2 1.9986 3/2*1.9986] - gps.y; % m
LAtx(3,:)   = -[ 0 0.05 0.1 0.15] - gps.z; % m
% LAtx(3,:)   = -[ 0 0 0 0] - gps.z; % m

% LArx(:,1)   = [0  0            0];    % inches
% LArx(:,2)   = [0 -1.5374/2    -0.05]; % inches
% LArx(:,3)   = [0 -1.5374      -0.1];  % inches
% LArx(:,4)   = [0 -3/2*1.5374  -0.15]; % inches

% LAtx(:,1)   = [0 0           0];    % inches
% LAtx(:,2)   = [0 1.5374/2   -0.05]; % inches
% LAtx(:,3)   = [0 1.5374     -0.1];  % inches
% LAtx(:,4)   = [0 3/2*1.5374 -0.15]; % inches

% LArx(1,:)   = LArx(1,:)*0.0254 - gps.x; % m, gps corrected
% LArx(2,:)   = LArx(2,:)*0.0254 - gps.y; % m, gps corrected
% LArx(3,:)   = LArx(3,:)*0.0254 - gps.z; % m, gps corrected
 
% LAtx(1,:)   = LAtx(1,:)*0.0254 - gps.x; % m, gps corrected
% LAtx(2,:)   = LAtx(2,:)*0.0254 - gps.y; % m, gps corrected
% LAtx(3,:)   = LAtx(3,:)*0.0254 - gps.z; % m, gps corrected

% Amplitude (not power) weightings for transmit side. 
if rxchannel == 0
  rxchannel = 1;
  tx_weights = ones(1,size(LAtx,2));
end
magsum       = sum(tx_weights);

% Weighted average of Xb, Yb and Zb components
LAtx_pc(1,1)    = dot(LAtx(1,:),tx_weights)/magsum;
LAtx_pc(2,1)    = dot(LAtx(2,:),tx_weights)/magsum;
LAtx_pc(3,1)    = dot(LAtx(3,:),tx_weights)/magsum;

phase_center = (LArx(:,rxchannel) + repmat(LAtx_pc,[1 size(rxchannel,2)]))./2;

return