function [param, frames, records, exec_good] = flightline_extract(param)

% [param, frames, records, exec_good] = flightline_extract(param)
%
% Function to extract or simulate a flightline to use with full simulator
%
% Author: John Paden, Hara Madhav Talasila
%
% See also sim.input_full_sim, run_load_data (example 7)

frames = [];
records = [];
exec_good = 0;

[WGS84] = physical_constants('WGS84');

%% Load params

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.sim.day_seg, datestr(now));
fprintf('=====================================================================\n');

params = read_param_xls(ct_filename_param(param,sprintf('%s_param_%s.xls',param.sim.radar_name,param.sim.season_name)));

% Find the specified day_seg in params spreadsheet
found_day_seg = 0;
for param_idx = 1:length(params)
  if strcmpi( params(param_idx).day_seg , param.sim.day_seg)
    param = merge_structs(params(param_idx), param);
    found_day_seg = 1;
    break;
  end
end

if ~found_day_seg; fprintf('day_seg not found\n'); return; end;

clear params found_day_seg param_idx;

%% Load frames

try
  fprintf('(%s) Loading frames',  datestr(now));
  frames = frames_load( param );
catch ME
  fprintf('\t-- Failed (%s)\n',  datestr(now)); return;
end
fprintf('\t-- Loaded (%s)\n',  datestr(now));

% override default variables in frames
frames.frame_idxs = 1;
frames.nyquist_zone = NaN; %1
frames.proc_mode = 0;
frames.quality = 1;

%% Load records

try
  fprintf('(%s) Loading records', datestr(now));
  records = records_load( param );
catch ME
  fprintf('\t-- Failed (%s)\n',  datestr(now)); return;
end
fprintf('\t-- Loaded (%s)\n',  datestr(now));

% assign start and stop indices of records to use in simulation
if isfield(param.sim,'start_gps_time') && isfield(param.sim,'stop_gps_time')
  start_idx = find(records.gps_time>=param.sim.start_gps_time,1,'first');
  stop_idx  = find(records.gps_time>=param.sim.stop_gps_time,1,'first');
elseif isfield(param.sim,'frame_idx') % select frame
  start_idx   = frames.frame_idxs(param.sim.frame_idx);
  stop_idx    = frames.frame_idxs(param.sim.frame_idx+1);
elseif 1  % just some default values
  start_idx = 1;
  stop_idx  = 1001;
end

rec_len = stop_idx-start_idx+1;

try % truncates records' structure (instead of read_records_aux_files)
  records = struct_truncate(records,length(records.gps_time),start_idx,stop_idx,0);
catch % fail-safe
  records.gps_time = records.gps_time(start_idx:stop_idx);
  records.lat = records.lat(start_idx:stop_idx);
  records.lon = records.lon(start_idx:stop_idx);
  records.elev = records.elev(start_idx:stop_idx);
  records.roll = records.roll(start_idx:stop_idx);
  records.pitch = records.pitch(start_idx:stop_idx);
  records.heading = records.heading(start_idx:stop_idx);
end

if rec_len~=length(records.gps_time)
  fprintf('Check start, stop idxs and length(records)\n');
  return;
end

% To use while loading simulated data
param.load_data.recs = [start_idx stop_idx];
param.load_data.imgs = param.sim.imgs;

clear start_idx stop_idx;

%%  Trajectory and Attitude (GPS+IMU)

param.gps = [];
param.gps_source = records.gps_source; % can be used for leverarm

if isfield(param.sim, 'north_along_track_en') && param.sim.north_along_track_en
  
  %% For Northward flightline
  % get along_track from trajectory without leverarm
  [along_track,~,~,~] = geodetic_to_along_track(records.lat, records.lon, records.elev); % no 'spacing'
  dx = median(diff(along_track));
  
  % Calculate unit vector towards North: geodetic (a, b) or ECEF (A, B)
  a       = [];
  a.lat   =records.lat(1);
  a.lon   =records.lon(1);
  a.elev  =records.elev(1);
  A       = nan(3,1);
  [A(1), A(2), A(3)] = geodeticD2ecef(a.lat, a.lon, a.elev, WGS84.ellipsoid);
  % B is ~1 meter North of A
  b     = a;
  b.lat = a.lat + 1/111111; % ~111111 meter per degree
  B     = nan(3,1);
  [B(1), B(2), B(3)] = geodeticD2ecef(b.lat, b.lon, b.elev, WGS84.ellipsoid);
  north_unit_vec = (B-A);
  north_unit_vec = north_unit_vec./norm(north_unit_vec);
  
  % Crete Northward way-points
  traj_ecef = bsxfun(@plus, A, bsxfun(@times, north_unit_vec, dx * (0:rec_len-1)) );
  
  % Populate alternative records info: rec(trajectory and attitude)
  rec = [];
  rec.x = traj_ecef(1,:);
  rec.y = traj_ecef(2,:);
  rec.z = traj_ecef(3,:);
  [rec.lat, rec.lon, rec.elev] = ecef2geodeticD(rec.x, rec.y, rec.z, WGS84.ellipsoid);
  rec.roll      = zeros(1,rec_len);
  rec.pitch     = zeros(1,rec_len);
  rec.heading   = zeros(1,rec_len); % North
  rec.gps_time  = records.gps_time;
  
  if 0 % check records vs rec
    figure;
    plot(records.lon, records.lat, 'x');
    hold on;
    plot(rec.lon, rec.lat, 'o');
    xlabel('Longitude'); ylabel('Latitude');
    grid on; legend({'records','rec'});
  end
  
  param.gps = rec;
  clear along_track dx a A b B north_unit_vec traj_ecef rec
  
  % Overwrite this to records file
  records = merge_structs(records, param.gps);
  
else
  
  %% For Reference flightline
  
  param.gps.lat       = records.lat;
  param.gps.lon       = records.lon;
  param.gps.elev      = records.elev;
  param.gps.roll      = records.roll;
  param.gps.pitch     = records.pitch;
  param.gps.heading   = records.heading;
  param.gps.gps_time  = records.gps_time;
  [param.gps.x, param.gps.y, param.gps.z] = geodeticD2ecef(param.gps.lat, param.gps.lon, param.gps.elev, WGS84.ellipsoid);
  
end

fprintf('(%s) Flightline \t-- Done\n',  datestr(now));

%% Target(s) location

param.target = [];

if ~isfield(param.sim,'target') || ~isfield(param.sim.target,'type') ...
    || strcmpi(param.sim.target.type, 'point')
  param.target.type = 'point'; % 'surface' %'layer'
end

switch param.target.type
  case 'point'
    mid_idx = ceil(rec_len/2);
    param.target.lat  = param.gps.lat(mid_idx);
    param.target.lon  = param.gps.lon(mid_idx);
    if isnan(param.gps.elev(mid_idx))
      param.target.elev = 0;
    elseif 1 % temporary override
      param.target.elev = 0;
    end
    clear mid_idx
end

[param.target.x, param.target.y, param.target.z] = geodeticD2ecef(param.target.lat, param.target.lon, param.target.elev, WGS84.ellipsoid);

fprintf('(%s) Target \t\t-- Done\n',  datestr(now));

%% Signal

param.signal = [];

Ntx = 1;
for idx = 1:Ntx
  param.signal.tx(Ntx).gain           = 1; %[]; % Nt x N_elev x N_azi
  param.signal.tx(Ntx).freq           = 1; %[]; % Nt x 1
  param.signal.tx(Ntx).elev_angle     = 1; %[]; % N_elev x 1
  param.signal.tx(Ntx).azimuth_angle  = 1; %[]; % N_azi x 1
end

Nrx = 1;
for idx = 1:Nrx
  param.signal.rx(Nrx).gain           = 1; %[]; % Nt x N_elev x N_azi
  param.signal.rx(Nrx).freq           = 1; %[]; % Nt x 1
  param.signal.rx(Nrx).elev_angle     = 1; %[]; % N_elev x 1
  param.signal.rx(Nrx).azimuth_angle  = 1; %[]; % N_azi x 1
end

param.sim.Ntx = Ntx;
param.sim.Nrx = Nrx;

%% Waveforms

param.load.imgs = param.sim.imgs;
[wfs,~] = data_load_wfs(param,records);

% overrides for each waveform
% for wf = 1:length(wfs)
%   [output_dir,radar_type,radar_name] = ct_output_dir(param.sim.radar_name);
%   if strcmpi(radar_type,'deramp')
%     wfs.deramp(wfs) = 1; % Default 0
%   elseif strcmpi(radar_type,'pulsed')
%     wfs.deramp(wf) = 0; % Default 0
%   end
%   wfs(wf).Nt_raw = wfs(wf).Nt;
% end

param.radar.wfs = merge_structs(param.radar.wfs,wfs);

fprintf('(%s) Waveforms \t-- Done\n',  datestr(now));

%% Closing

param.season_name = sprintf('%ssim',param.season_name);

exec_good = 1;

return;


end