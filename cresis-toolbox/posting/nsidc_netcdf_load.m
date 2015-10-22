function echo = nsidc_netcdf_load(fn)
% echo = nsidc_netcdf_load(fn)
%
% fn = path to radar echogram netcdf file
%
% echo = structure containing file contents (optional)
%   If there is no output argument, then variables are stored in the
%   caller's workspace (as Matlab's load function does).
%  .Data = Nt by Nx radar echogram matrix (relative power)
%  .Latitude = 1 by Nx vector, latitude (deg)
%  .Longitude = 1 by Nx vector, longitude (deg)
%  .Elevation = 1 by Nx vector, WGS-84 elevation (m)
%  .GPS_time = 1 by Nx vector, GPS time in ANSI-C standard (sec since Jan 1, 1970)
%  .Time = Nt by 1 vector, radar echogram fast-time axis
%
% Three ways to call nsidc_netcdf_load
%
% fn = 'IRACC1B_V01_20130426_01_100.nc';
% echo = nsidc_netcdf_load(fn);
%
% fn = 'IRACC1B_V01_20130426_01_100.nc';
% nsidc_netcdf_load(fn);
%
% nsidc_netcdf_load IRACC1B_V01_20130426_01_100.nc
%
% Author: John Paden
%
% See also: type "nsidc_help.m"

ncid = netcdf.open(fn,'NOWRITE');

try
  varid = netcdf.inqVarID(ncid,'amplitude');
  echo.Data = 10.^(netcdf.getVar(ncid,varid,'double')/10);
  
  varid = netcdf.inqVarID(ncid,'lat');
  echo.Latitude = netcdf.getVar(ncid,varid,'double');
  echo.Latitude = reshape(echo.Latitude, [1 length(echo.Latitude)]);
  
  varid = netcdf.inqVarID(ncid,'lon');
  echo.Longitude = netcdf.getVar(ncid,varid,'double');
  echo.Longitude = reshape(echo.Longitude, [1 length(echo.Longitude)]);
  
  varid = netcdf.inqVarID(ncid,'altitude');
  echo.Elevation = netcdf.getVar(ncid,varid,'double');
  echo.Elevation = reshape(echo.Elevation, [1 length(echo.Elevation)]);
  
  varid = netcdf.inqVarID(ncid,'time');
  GPS_time_date = netcdf.getAtt(ncid,varid,'units');
  year = str2double(GPS_time_date(15:18));
  month = str2double(GPS_time_date(20:21));
  day = str2double(GPS_time_date(23:24));
  echo.GPS_time = netcdf.getVar(ncid,varid,'double');
  echo.GPS_time = echo.GPS_time + datenum_to_epoch(datenum(year,month,day));
  echo.GPS_time = echo.GPS_time + utc_leap_seconds(echo.GPS_time(1));
  echo.GPS_time = reshape(echo.GPS_time, [1 length(echo.GPS_time)]);
  
  varid = netcdf.inqVarID(ncid,'fasttime');
  echo.Time = netcdf.getVar(ncid,varid,'double') * 1e-6;
  echo.Time = reshape(echo.Time, [length(echo.Time) 1]);
catch ME
  netcdf.close(ncid);
  rethrow(ME);
end

netcdf.close(ncid);

if nargout == 0
  assignin('caller','Data',echo.Data);
  assignin('caller','Latitude',echo.Latitude);
  assignin('caller','Longitude',echo.Longitude);
  assignin('caller','Elevation',echo.Elevation);
  assignin('caller','GPS_time',echo.GPS_time);
  assignin('caller','Time',echo.Time);
  clear echo;
end

return;
