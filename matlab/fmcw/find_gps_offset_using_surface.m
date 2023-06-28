% script find_gps_offset_using_surface
%
% This is a script which helps find the vectors.gps.time_offset field.
% It is primarily useful over ocean where the surface is known to be flat.
%

error('USE run_fmcw_set_nyquist_zone_Tsys_GPSoffset.m IN PLACE OF THIS SCRIPT UNLESS YOU DO NOT HAVE THE REQUIRED INPUTS TO RUN THAT SCRIPT.\n');

if 1
  % Only need to load things once
%   records = load('/cresis/projects/dev/cr1/records/kuband/2009_Antarctica_DC8/records_20091018_01.mat');
%   gps = load('/cresis/projects/dev/cr1/gps/2009_Antarctica_DC8/gps_20091018.mat');
  records = load('/scratch/csarp_support/records/snow/2017_Greenland_P3/records_20170309_01.mat');
  gps = load('/scratch/csarp_support/gps/2017_Greenland_P3/gps_20170309.mat');

  % Load ATM data

end

%% Plot
figure(1); clf;

% height_offset added to the red (radar) line
height_offset = -23.3+42;
gps_offset = 0;
nz_sign = 1; % Change to -1 if Nyquist zone is flipped

new_elev = interp1(gps.gps_time, gps.elev, records.gps_time + gps_offset);

hold on;
plot(gps.gps_time-gps.gps_time(1), gps.elev);
plot(records.gps_time-gps.gps_time(1) + gps_offset, +nz_sign*records.surface*3e8/2 + height_offset, 'r');

h_surface = plot(records.gps_time-gps.gps_time(1) + gps_offset, new_elev - nz_sign*records.surface*3e8/2 + height_offset, 'g');
grid on;
hold off;

return;

% Update Nyquist Zones in frames file using elevation data (NOTE: only works for sea ice where
% surface is around zero elevation):
param_fn = ct_filename_param('snow_param_2010_Greenland_DC8.xls')
param = read_param_xls(param_fn,'20100323_05','post');
fmcw_set_nyquist_zone_from_elev(param, records.elev/(3e8/2), 2);

param_fn = ct_filename_param('snow_param_2010_Greenland_DC8.xls')
param = read_param_xls(param_fn,'20100421_02','post');
fmcw_set_nyquist_zone_from_elev(param, records.elev/(3e8/2), 2);

% Iteration
gps_offset = 0;
set(h_surface,'XData',records.gps_time + gps_offset);

height_offset = -8.5;
set(h_surface,'YData',records.surface*3e8/2 + height_offset);
