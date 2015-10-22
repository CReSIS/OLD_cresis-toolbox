function records = check_records(fn)
%
% Checks for common problems in mcords records files that may require
% recreating the records files.
%
% Author: John Paden

if ~exist('fn','var') || isempty(fn)
  fn = 'C:\csarp_support\records\mcords2\2011_Greenland_P3\records_20110317_01.mat';
end

load(fn);

if ~isfield(records,'gps_source')
  fprintf('ERROR: no gps_source\n');
end

% close all;

if isfield(records,'lat')
  figure(1); clf;
  plot(records.lon,records.lat,'.-');
  xlabel('Longitude (deg,E');
  ylabel('Latitude (deg,N)');
  title('Check for discontinuities');

  if any(isnan(records.lat))
    fprintf('Lat contains NaN... bad!\n');
  end

  figure(2); clf;
  plot(records.elev)
  ylabel('Elevation (m)');
  title('Check for discontinuities');

  figure(3); clf;
  plot(diff(double(records.gps_time)))
  title('diff(double(records.gps\_time))), ideally around 2-4 msec');
  ylabel('EPRI (sec)');

  if any(diff(double(records.gps_time)) <= 0)
    fprintf('gps_time goes backwards... bad!\n');
  end
end

figure(4); clf;
plot(records.time)
title('records.time, seconds of day (0-86400)');
ylabel('Time (sec)');

figure(5); clf;
plot(records.gps_time,'-')
title('Seconds since Jan 1, 1970, should be ~1.Xe9');
ylabel('Time (sec)');

if any(diff(double(records.time)) <= 0)
  fprintf('time goes backwards... bad!\n');
end

figure(6); clf;
plot(diff(double(records.epri)),'.')
title('diff(double(records.epri)), ideally all ones');

if any(diff(double(records.epri)) <= 0)
  fprintf('epri goes backwards... bad!\n');
end



return;
