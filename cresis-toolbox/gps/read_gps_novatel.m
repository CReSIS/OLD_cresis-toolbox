function gps = read_gps_novatel(fn,param)
% gps = read_gps_novatel(fn,param)
%
% This function reads the DGPS/INS data output from the Waypoint Inertial
% Explorer from NovAtel, Canada.  This is documented in the report
%   'CReSIS GPS+INS Post-Processing Guide'
% by Huan Zhao and Chris Gifford.
%
% gps = structure of the output
% param =
%  .time_reference = 'gps' or 'utc'
%
% Examples:
%  fn = '/cresis/projects/metadata/2009_Antarctica_TO/rover_diff_antarctica_01062010_1.txt';
%  gps = read_gps_novatel_rpygga(fn, struct('time_reference','gps'));
%  plot_gps(gps);
%
% Author: Huan Zhao, John Paden
%
% See also read_gps_applanix, read_gps_atm, read_gps_csv, read_gps_litton,
%   read_gps_nmea, read_gps_novatel, read_gps_reveal, read_gps_traj,
%   read_gps_txt, plot_gps


[fid,msg] = fopen(fn,'r');
if fid < 1
  fprintf('Could not open file %s\n', fn{fn_idx});
  error(msg);
end

% Format 2011 Greenland TO:
% 400989.00    3/24/2011     69 12 10.87930    -51 26 56.16094       1648.423     8.1594     0     0.012 1.1392  0     0.01243   243.8001  0  0.0309800003
% GPS-SOW      DATE          Lat deg,min,sec   Lon deg,min,sec       elev m       ?/pitch    NA    NA    ?/roll       NA    NA        Heading   NA NA

% Format:
% 444159.01  1/01/2010  -79 28 08.07002 -112 02 58.27124 \
%   1765.103  -1.0213310000 359.5417000055   0.0056421673  \
%  -0.1392850000   3.5313000679   0.0057160472  46.3112970000 \
%   1.0263999701   0.0193856582
% 1. Seconds of week (GPS time)
%  Field 1
% 2. Date (GPS time)
%  Field 2
% 3. Latitude (deg,min,sec)
%  Field 3-5
% 4. Longitude (deg,min,sec)
%  Field 6-8
% 5. Elevation (m) with respect to WGS-84 ellipsoid
%  Field 9
% 6. Pitch (deg) [value, separation?, standard deviation]
%  Field 10-12
% 7. Roll (deg) [value, separation?, standard deviation]
%  Field 13-15
% 8. Heading (deg) [value, separation?, standard deviation]
%  Field 16-18

A = textscan(fid,'%f %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'Headerlines', 20);

fclose(fid);

% Convert SOW and date into absolute time
dates = datenum(A{2});
gps.gps_time = A{1};

% Fix bug in Novatel files... at the week wrap of the seconds of week,
% the date is bad.
diff_time = diff(gps.gps_time);
median_diff_time = median(diff_time);
bad_idxs = [0; abs(diff_time - median_diff_time) > 2.5];
if ~isempty(bad_idxs)
  warning('Found some large gps time jumps in seconds of week field, discarding\n');
end
dates = dates(~bad_idxs);
gps.gps_time = gps.gps_time(~bad_idxs);

diff_dates = diff(dates);
bad_idx = find(diff_dates < 0)+1;
if ~isempty(bad_idx)
  warning('Found a week wrap in the dates/seconds of week field\n');
  dates(bad_idx:end) = dates(bad_idx:end) + 7;
end

gps.gps_time = gps_sow_to_epoch(gps.gps_time,dates);

if strcmpi(param.time_reference,'utc')
  warning('Normally these files are GPS time reference');
  % UTC time stored in file, so need to add leap seconds back in
  gps.gps_time = gps.gps_time + utc_leap_seconds(gps.gps_time(1));
elseif strcmpi(param.time_reference,'gps')
else
  error('You must specify a time reference');
end

gps.lat = sign(A{3}).*(abs(A{3}) + A{4}/60 + A{5}/3600);
gps.lon = sign(A{6}).*(abs(A{6}) + A{7}/60 + A{8}/3600);

gps.elev = A{9};

% gps.pitch = A{13}/180*pi;
% gps.roll = A{16}/180*pi;
% gps.heading = A{10}/180*pi;
gps.pitch = A{10}/180*pi;
gps.roll = A{13}/180*pi;
gps.heading = A{16}/180*pi;


gps.gps_time = reshape(gps.gps_time,[1 length(gps.gps_time)]);
gps.lat = reshape(gps.lat,[1 length(gps.lat)]);
gps.lon = reshape(gps.lon,[1 length(gps.lon)]);
gps.elev = reshape(gps.elev,[1 length(gps.elev)]);
gps.roll = reshape(gps.roll,[1 length(gps.roll)]);
gps.pitch = reshape(gps.pitch,[1 length(gps.pitch)]);
gps.heading = reshape(gps.heading,[1 length(gps.heading)]);

gps.lat = gps.lat(~bad_idxs);
gps.lon = gps.lon(~bad_idxs);
gps.elev = gps.elev(~bad_idxs);
gps.roll = gps.roll(~bad_idxs);
gps.pitch = gps.pitch(~bad_idxs);
gps.heading = gps.heading(~bad_idxs);

return;
