function gps_time = utc_to_gps(utc_time)
% gps_time = utc_to_gps(utc_time)
% 
% Converts from ANSI-C UTC time to ANSI-C GPS time. In other words, leap
% seconds are removed.
%
% Example:
% gps_time = utc_to_gps(utc_time);
%
% Author: John Paden
%
% See also: gps_to_utc.m, utc_to_gps.m, utc_leap_seconds.m, epoch_to_datenum.m,
% datenum_to_epoch.m

gps_time = utc_time + utc_leap_seconds(utc_time(1));

return;
