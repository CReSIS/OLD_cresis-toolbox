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
% See also epoch_to_datenum.m, datenum_to_epoch.m, epoch_to_gps_sow.m,
%   epoch_to_sod.m, sod_to_epoch.m, gps_sow_to_epoch.m, utc_leap_seconds.m,
%   gps_to_utc.m, utc_to_gps.m

if isempty(utc_time)
  gps_time = [];
else
  gps_time = utc_time + utc_leap_seconds(utc_time(1));
end

return;
