function [gps_time,cpu_time] = read_cpu_time(fn)
% [gps_time,cpu_time] = read_cpu_time(fn)
%
% Read in CSV file which contains GPS time and computer/CPU time. The first
% line should have the headers. Subsequent lines should be GPS time
% followed by CPU time.
%
% Note that both the GPS_TIME and CPU_TIME need to be in the same reference
% (e.g. both should be in UTC-time or both should be in GPS-atomic-time).
%
% The unix command for creating the CPU date field is:
% date +"%x %T %N"
% which correponds to:
% "DATE TIME NANOSECONDS"
%
% This file is for situations where the GPS receiver is not working.
%
% GPS_TIME,CPU_TIME
% 12/17/2018 07:23:20, 12/17/2018 06:48:21 281923350
% 12/17/2018 08:20:00, 12/17/2018 07:45:02 234235234
% 12/17/2018 09:10:00, 12/17/2018 08:35:03 325423643
% 12/17/2018 10:20:40, 12/17/2018 09:45:44 364582334
% 12/17/2018 11:30:50, 12/17/2018 10:55:55 056445653
%
% Example:
%
% fn = '~/cpu_time_20181217.csv';
% [gps_time,cpu_time] = read_cpu_time(fn);
% figure;
% plot(gps_time-gps_time(1), gps_time-cpu_time,'.');
% xlabel('Relative GPS time (sec)');
% ylabel('Drift (sec)');
%
% fns = get_filenames(fullfile(in_base_path,'CPU_TIME'),'cpu_time','','.csv',struct('recursive',true));
% [gps_time,cpu_time] = read_cpu_time(fns);
% gps_time_origin = gps_time(round(end/2));
% pp = polyfit(gps_time-gps_time_origin, gps_time-cpu_time, 5);
% h_fig = figure;
% plot(gps_time-gps_time_origin, gps_time-cpu_time,'.');
% hold on;
% grid on;
% plot(gps_time-gps_time_origin, polyval(pp, gps_time-gps_time_origin));
% plot(gps_time-gps_time_origin, cpu_time + polyval(pp, cpu_time - gps_time_origin) - gps_time,'x');
% xlabel('Relative GPS time (sec)');
% ylabel('Drift (sec)');
% correction.gps_time_origin = gps_time_origin;
% correction.pp = pp;
% set(h_fig,'UserData',pp)
%
% Author: John Paden
%
% See also:

if iscell(fn)
  % fn: cell array of filenames to load
  gps_time = [];
  cpu_time = [];
  for fn_idx = 1:length(fn)
    % Recurse on each filename to load it
    [tmp_gps_time,tmp_cpu_time] = read_cpu_time(fn{fn_idx});
    gps_time(end+(1:length(tmp_gps_time))) = tmp_gps_time;
    cpu_time(end+(1:length(tmp_cpu_time))) = tmp_cpu_time;
  end
  
else
  % fn: string of CSV filename containing cpu and gps times
  [fid,msg] = fopen(fn,'rb');
  if fid<0
    error('Failed to open %s:\n  %s', fn, msg);
  end
  
  C = textscan(fid,'%s%s','Delimiter',',','HeaderLines',1);
  
  fclose(fid);
  
  gps_time = nan(1,length(C{2}));
  cpu_time = nan(1,length(C{2}));
  for idx = 1:length(C{2})
    try
      gps_time(idx) = datenum(C{1}{idx});
      s = regexp(C{2}{idx}, '\s+', 'split');
      % Expecting s to be of the form:
      %   '12/17/2018'    '06:48:21'    '281923350'
      cpu_time(idx) = datenum([s{1} ' ' s{2}]) + str2double(s{3})/1e9/86400;
    catch ME
      ME.getReport
    end
  end
  
  % Switch from Matlab's datenum format (days since 1900) to ANSI C standard
  % (seconds since 1970)
  gps_time = datenum_to_epoch(gps_time);
  cpu_time = datenum_to_epoch(cpu_time);
end
