function gps = read_gps_arena_cpu_time(fn, param)
% gps = read_gps_arena_cpu_time(fn, param)
%
% Reads in Arena log files that includes CPU time and radar time.
%
% File format should include relTimeCntr and profileCntr.
%
% relTimeCntr:15396324470892465
% profileCntr:893
% ppsFracCntr:46704704
% ppsCntr:1539632447
% ZynqCoreTemp:69.40:C
% Board1Temp:63.00:C
% Board2Temp:63.50:C
%
% Input Args:
%   fn: string containing input Arena GPS filename
%     e.g. 20180817_094746_ARENA-CTU-ctu-gps.txt
%   param: tells the file GPS type and the year, month, day to determine
%     absolute time (GGA NMEA files just give the time of day)
%     .year
%     .month
%     .day
%     .time_reference: 'gps' or 'utc' (should always be 'utc')
%     .nmea_tag: NEMA string to identify good lines (e.g. '$GPGGA')
%     .clk: clock for radar_time counter, defaults to 10e6
%
% Output Args:
% gps: output structure with fields
%  .time: NaN (sec)
%  .lat: NaN (deg)
%  .lon: NaN (deg)
%  .elev: NaN (m)
%  .roll: NaN (rad)
%  .pitch: NaN (rad)
%  .heading: NaN (rad)
%  .relTimeCntr: radar time, 64 bit counter free-running at 10 MHz
%  .profileCntr: pulse counter
%  .ppsCntr: PPS counter
%
% Example:
%   fn = '/cresis/snfs1/dataproducts/metadata/2018_Antarctica_Ground/20181015/20181015_144047_ARENA0-awg0.txt';
%   gps = read_gps_arena_cpu_time(fn);
%   plot(gps.radar_time,gps.profileCntr);
%   datestr(epoch_to_datenum(gps.radar_time(1)));
%
% Author: John Paden
%
% See also read_gps_applanix, read_gps_atm, read_gps_csv, read_gps_litton,
%   read_gps_nmea, read_gps_novatel, read_gps_reveal, read_gps_traj,
%   read_gps_txt, plot_gps

if ~isfield(param,'clk') || isempty(param.clk)
  param.clk = 10e6;
end

[fid,msg] = fopen(fn,'r');
if fid < 0
  error('Error opening %s: %s', fn, msg);
end

relTimeCntr = [];
profileCntr = [];
ppsCntr = [];
ppsFracCntr = [];

record_idx = 1;
relTimeCntrTmp = NaN;
profileCntrTmp = NaN;
ppsFracCntrTmp = NaN;
while ~feof(fid)
  str = fgets(fid);
  [token,remain] = strtok(str,':');
  if strcmpi(token,'relTimeCntr')
    relTimeCntrTmp = str2double(remain(2:end));
  elseif strcmpi(token,'profileCntr')
    profileCntrTmp = str2double(remain(2:end));
  elseif strcmpi(token,'ppsFracCntr')
    ppsFracCntrTmp = str2double(remain(2:end));
  elseif strcmpi(token,'ppsCntr')
    ppsCntrTmp = str2double(remain(2:end));
    relTimeCntr(record_idx) = relTimeCntrTmp;
    profileCntr(record_idx) = profileCntrTmp;
    ppsCntr(record_idx) = ppsCntrTmp;
    ppsFracCntr(record_idx) = ppsFracCntrTmp;
    record_idx = record_idx + 1;
    relTimeCntrTmp = NaN;
    profileCntrTmp = NaN;
    ppsFracCntrTmp = NaN;
  end
end
fclose(fid);

if record_idx < 2
  gps.gps_time = [];
  gps.lat = [];
  gps.lon = [];
  gps.elev = [];
  gps.roll = [];
  gps.pitch = [];
  gps.heading = [];
  gps.radar_time = [];
  gps.profileCntr = [];
  gps.comp_time = [];
  return;
end

% ===================================================================
% Store outputs in structure
% ===================================================================

gps.radar_time = relTimeCntr/param.clk;
gps.profileCntr = profileCntr;
gps.comp_time = ppsCntr + ppsFracCntrTmp/param.clk;

% Create gps_time from the radar_time variable
% -------------------------------------------------------------------------
if strcmpi(param.time_reference,'utc')
  % UTC time stored in file, so need to add leap seconds back in
  if ~isempty(gps.radar_time)
    gps.gps_time = gps.radar_time + utc_leap_seconds(gps.radar_time(1));
  else
    gps.gps_time = [];
  end
else
  gps.gps_time = gps.radar_time;
end
if isfield(param,'cpu_time_fn') && exist(param.cpu_time_fn,'file')
  % Apply a correction to the GPS time
  [gps_time,cpu_time] = read_cpu_time(param.cpu_time_fn);
  if strcmpi(param.cpu_time_ref,'utc')
    % UTC time stored in file, so need to add leap seconds back in
    if ~isempty(gps_time)
      gps_time = gps_time + utc_leap_seconds(gps_time(1));
    else
      gps_time = [];
    end
  end
  if length(cpu_time) == 1
    gps.gps_time = gps.gps_time + gps_time - cpu_time;
  else
    gps.gps_time = gps.gps_time + interp1(cpu_time, gps_time - cpu_time, gps.gps_time, 'linear', 'extrap');
  end
end

% Fill in remaining fields with nan
gps.lat = nan(size(gps.radar_time));
gps.lon = nan(size(gps.radar_time));
gps.elev = nan(size(gps.radar_time));
gps.roll = nan(size(gps.radar_time));
gps.pitch = nan(size(gps.radar_time));
gps.heading = nan(size(gps.radar_time));

if all(gps.radar_time(1)==gps.radar_time)
  warning('gps.radar_time is constant. GPS 1 PPS may not have been received. Radar synchronization to GPS will not be possible with this file.');
end

if all(gps.comp_time(1)==gps.comp_time)
  warning('gps.comp_time is constant. GPS 1 PPS may not have been received.');
end
