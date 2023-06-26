function utc_time = gps_to_utc(gps_time)
% utc_time = gps_to_utc(gps_time)
% 
% Converts from ANSI-C GPS time to ANSI-C UTC time. In other words, leap
% seconds are inserted.
%
% Example:
% utc_time = gps_to_utc(gps_time);
%
% Author: John Paden
%
% See also epoch_to_datenum.m, datenum_to_epoch.m, epoch_to_gps_sow.m,
%   epoch_to_sod.m, sod_to_epoch.m, gps_sow_to_epoch.m, utc_leap_seconds.m,
%   gps_to_utc.m, utc_to_gps.m


if isempty(gps_time)
  utc_time = [];
else
  utc_time = gps_time - utc_leap_seconds(gps_time(1));
end

return;
