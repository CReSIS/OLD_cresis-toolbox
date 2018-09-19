function time = gps_sow_to_epoch(sow,param)
% time = gps_sow_to_epoch(sow,param)
%
% Converts from GPS seconds of the week to seconds since the 
% Jan 1, 1970 epoch. The assumption is that sow DO NOT WRAP
% AROUND at 604800 seconds (the number of seconds in a week).
%
% sow = GPS seconds of the week
% param = datenum OR struct that gives absolute reference to find which
%   GPS week we are in.
%  .year
%  .month
%  .day
%
% time = seconds since Jan 1, 1970 epoch
%
% Example:
%  read_gps_applanix.m
%
% See also epoch_to_datenum.m, datenum_to_epoch.m, epoch_to_gps_sow.m,
%   epoch_to_sod.m, gps_sow_to_epoch.m, utc_leap_seconds.m

% GPS weeks started counting Jan 6, 1980
start_GPS = datenum(1980,1,06);

% Current time is:
if isstruct(param)
  cur_GPS = datenum(param.year,param.month,param.day);
else
  % param is a datenum already (class is double)
  cur_GPS = param;
end

% For debugging, the broadcasted GPS week is mod 1024
if 0
  % Determine the absolute GPS week
  gps_week = floor((cur_GPS-start_GPS)/7);

  gps_week_broadcast = mod(gps_week,1024);
end

% Find the beginning of the current GPS week
ref_GPS = cur_GPS - mod(cur_GPS-start_GPS,7);

% Add seconds of week converted to days to the ref_GPS
time = ref_GPS + sow/86400;

% Convert from Matlab datenum to computer epoch
time = datenum_to_epoch(time);

return;
