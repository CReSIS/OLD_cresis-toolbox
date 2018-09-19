function time = epoch_to_datenum(time)
% time = epoch_to_datenum(time)
%
% Convert from C library seconds since Jan 1, 1970 epoch to Matlab's
% datenum.
%
% Author: John Paden
%
% See also epoch_to_datenum.m, datenum_to_epoch.m, epoch_to_gps_sow.m,
%   epoch_to_sod.m, gps_sow_to_epoch.m, utc_leap_seconds.m

time = datenum(1970,1,1,0,0,time);

return;
