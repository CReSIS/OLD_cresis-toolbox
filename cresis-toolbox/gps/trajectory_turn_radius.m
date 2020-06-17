function radius = trajectory_turn_radius(speed,bank_angle)
% radius = trajectory_turn_radius(speed,bank_angle)
%
% Minimum turn radius for an aircraft for a given speed and bank angle.
%
% speed: matrix of speed (m/s)
% bank_angle: matrix of bank angles (radians) corresponding to speeds
%
% radius: matrix of minimum turn radius corresponding to inputs
%
% Author: John Paden

g = 9.80665;

radius = speed.^2 ./ (g*tan(bank_angle));
