function gps = trajectory_with_leverarm(gps,param)
% gps = trajectory_with_leverarm(gps,param)
%
% Takes a struct, gps, and updates the fields lat, lon, elev, roll,
% pitch, heading with the lever arm offset. All other fields of the struct
% are left unchanged.
%
% gps = struct
%  .lat = N by 1 vector of doubles (latitude north in degrees)
%  .lon = N by 1 vector of doubles (longitude east in degrees)
%  .elev = N by 1 vector of doubles (elevation relative to WGS-84 in meters)
%  .roll = N by 1 vector of doubles (roll in radians)
%  .pitch = N by 1 vector of doubles (pitch in radians)
%  .heading = N by 1 vector of doubles (true heading in radians,
%     relative to north)
% param = struct
%  .lever_arm_fh = string containing lever arm function name
%  .rx_path = receiver paths to use in lever arm function
%     If multiple paths are passed in, this function takes the average
%     lever arm position of all the receiver paths.
%  .tx_weights = transmit weights to use in lever arm function
%
% Author: John Paden

if isempty(param.lever_arm_fh)
  return;
end

physical_constants;

% Get the phase center of the receiver in question in BCS
[base_phase_center] = param.lever_arm_fh(param, param.tx_weights, param.rx_path);
base_phase_center = mean(base_phase_center,2);

phase_center = zeros(3,length(gps.roll));
for rline = 1:length(gps.roll)
  % Roll, pitch, heading
  rotation_mat ...
    = [cos(gps.heading(rline)), -1*sin(gps.heading(rline)), 0; sin(gps.heading(rline)), cos(gps.heading(rline)), 0; 0, 0, 1] ...
    * [cos(gps.pitch(rline)), 0, sin(gps.pitch(rline)); 0, 1, 0; -1*sin(gps.pitch(rline)), 0 cos(gps.pitch(rline))] ...
    * [1, 0, 0; 0, cos(gps.roll(rline)), -1 * sin(gps.roll(rline)); 0, sin(gps.roll(rline)), cos(gps.roll(rline))];
  
  % Middle and receiver phase centers with roll or roll/pitch only
  % compensation in NED coordinate system
  phase_center(:,rline) = rotation_mat * base_phase_center;
  
  % Convert from NED to ECEF (lv2ecef input is ENU, so we transpose the first
  % two inputs and negate the third)
  [phase_center(1,rline),phase_center(2,rline),phase_center(3,rline)] ...
    = lv2ecef(phase_center(2,rline),phase_center(1,rline),-phase_center(3,rline), ...
    gps.lat(rline)/180*pi,gps.lon(rline)/180*pi,gps.elev(rline),WGS84.ellipsoid);
end

% Convert from ECEF to geodetic
[gps.lat,gps.lon,gps.elev] = ecef2geodetic(phase_center(1,:),phase_center(2,:),phase_center(3,:),WGS84.ellipsoid);
gps.lat = gps.lat*180/pi;
gps.lon = gps.lon*180/pi;
