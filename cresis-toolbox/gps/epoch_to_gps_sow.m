function [sow,gps_week] = epoch_to_gps_sow(time)
% [sow,gps_week] = epoch_to_gps_sow(time)
%
% Converts from seconds since the Jan 1, 1970 epoch to GPS seconds
% of the week. The mod 1024 on GPS week is not done.
%
% time = seconds since Jan 1, 1970 epoch
%
% sow = GPS seconds of the week
% gps_week = GPS week
%
% Author: John Paden
%
% See also epoch_to_datenum.m, datenum_to_epoch.m, epoch_to_gps_sow.m,
%   epoch_to_sod.m, sod_to_epoch.m, gps_sow_to_epoch.m, utc_leap_seconds.m,
%   gps_to_utc.m, utc_to_gps.m

% GPS weeks started counting Jan 6, 1980
start_GPS = datenum(1980,1,06);

% Current time is:
cur_GPS = epoch_to_datenum(time);

gps_week = floor((cur_GPS-start_GPS)/7);

% Find the beginning of the current GPS week
ref_GPS = cur_GPS - mod(cur_GPS-start_GPS,7);

sow = (cur_GPS-ref_GPS) * 86400;

return;
