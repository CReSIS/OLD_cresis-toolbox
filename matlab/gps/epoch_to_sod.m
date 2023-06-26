function sod = epoch_to_sod(time,ref_time)
% sod = epoch_to_sod(time,ref_time)
%
% Epoch (ANSI-C time, seconds since Jan 1, 1970) converted to seconds of
% day.
%
% time = vector of numeric time
% time_ref = optional (default reference is the start of the day of the
%   time variable so that sod is seconds of day relative to the current
%   day of each element in time). Two formats for time_ref
%   'YYYYMMDD': Does SOD relative to the date specified in the string
%   scalar containing epoch time: Does SOD relative to the start of day
%     that this epoch time falls within
%
% sod = vector of size equal to time filled with seconds of day for
%   each element in time
%
% Author: John Paden
%
% See also epoch_to_datenum.m, datenum_to_epoch.m, epoch_to_gps_sow.m,
%   epoch_to_sod.m, sod_to_epoch.m, gps_sow_to_epoch.m, utc_leap_seconds.m,
%   gps_to_utc.m, utc_to_gps.m

% Convert epoch to datenum (datenum is days since Jan 1, 1900)
time = epoch_to_datenum(time);

if ~exist('ref_time','var') || isempty(ref_time)
  %% Use the day of each element of time to determine day reference
  % Since time is in days, we just need to remove the decimal part
  % to get the time at the beginning of the day for the first time
  % sample.
  ref_time = floor(time);
  
elseif ischar(ref_time)
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

% Find the number of days relative to this reference time and then
% convert this to seconds of day. ref_time is a scalar or vector equal
% in size to time depending on the type of ref_time specified
sod = 86400 * (time - ref_time);

return;
