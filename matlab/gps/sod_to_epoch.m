function time = sod_to_epoch(sod,ref_time)
% time = sod_to_epoch(sod,ref_time)
%
% Seconds of day with reference time converted to epoch (ANSI-C time,
% seconds since Jan 1, 1970).
%
% sod = vector of seconds of day
% time_ref = Two formats for time_ref
%   'YYYYMMDD': Does SOD relative to the date specified in the string
%   scalar containing epoch time: Does SOD relative to the start of day
%     that this epoch time falls within
%
% time = vector of size equal to sod filled with epoch time for
%   each element in sod
%
% Author: John Paden
%
% See also epoch_to_datenum.m, datenum_to_epoch.m, epoch_to_gps_sow.m,
%   epoch_to_sod.m, sod_to_epoch.m, gps_sow_to_epoch.m, utc_leap_seconds.m,
%   gps_to_utc.m, utc_to_gps.m

if ischar(ref_time)
  %% Use the string YYYYMMDD to determine day reference
  ref_time = datenum(str2double(ref_time(1:4)),str2double(ref_time(5:6)), ...
    str2double(ref_time(7:8)));
  
else
  %% Use the scalar ref_time to determine day reference
  % Convert epoch to datenum (datenum is days since Jan 1, 1900)
  ref_time = epoch_to_datenum(ref_time(1));
  % Since time is in days, we just need to remove the decimal part
  % to get the time at the beginning of the day for the first time
  % sample.
  ref_time = floor(ref_time);
end

% Convert the ref_time from datenum to epoch and then add the seconds of
% the day.
time = epoch_to_datenum(ref_time) + sod;
