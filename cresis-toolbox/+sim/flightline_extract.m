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
  north_unit_vec = north_unit_vec./norm(north_unit_vec); % vecnorm
  
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
  
  if 1
    % check records vs rec
    figure;
    plot(records.lon, records.lat, 'x');
    hold on;
    plot(rec.lon, rec.lat, 'o');
    xlabel('Longitude'); ylabel('Latitude');
    grid on; legend({'records','rec'});
    
    % check 3d plot
    figure;
    plot3(rec.x, rec.y, rec.z); hold on;
    grid on;
    
    % check northward unit vector
    figure;
    traj_ecef_diff = diff(traj_ecef,[],2); % Point2Point vectors
    traj_ecef_diff = traj_ecef_diff ./vecnorm(traj_ecef_diff); % P2P unit vectors
    test_unit_vec   = traj_ecef_diff - north_unit_vec;
    imagesc(test_unit_vec); colorbar;
    title('error should be reasonable (1e-10)');
    yticks([1, 2, 3]);
    yticklabels({'x', 'y', 'z'});
    ylabel('p2p vs north unitvector error');
    
  end
  
  % ENU
  if 0
    P = [rec.x(1); rec.y(1); rec.z(1)];
    P = traj_ecef(:,1);
    % https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_ECEF_to_ENU
    % if subscript(r) is origin or center of earth (ECEF)
    Rot_Mat_ECEF2ENU = [0 1 0; 0 0 1; 1 0 0];
    P_enu = Rot_Mat_ECEF2ENU * P;
  end
  
  param.gps = rec;
  
  % Overwrite this to records file
  records = merge_structs(records, param.gps);
  
  clear along_track dx a A b B north_unit_vec traj_ecef rec
  
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

% For simple target position method:
%
% 1. Compute along-track vector (geodetic_to_along_track)
%
% 2. Find ECEF of trajectory
%
% 3. For each point in the trajectory, find the XYZ unit vectors in the
% flight coordinate system.
%
% 4. Note that for each target, we have defined x,y,z position which is in flight
% coordinate system.
%
% 5. For each target, find the closest in along-track trajectory position
% to that target. [~,idx] = min(abs(target.x - traj.x));
%
% 6. Use XYZ FCS of this closest trajectory position to place the target.
% Convert target's coordinates to ECEF.

% Simulator:
% Operate in ECEF coordinates and just cycle through each target and each
% position computing ranges using Pythagorean's theorem

switch param.target.type
  case 'point'
    mid_idx = floor(rec_len/2);
    
    if isfield(param.sim, 'north_along_track_en') && param.sim.north_along_track_en
      % M is midpoint of northward track and MT is an orthogonal through 
      % the target (T). T is 500 meters from M towards earth's center
      M = [param.gps.x(mid_idx); param.gps.y(mid_idx); param.gps.z(mid_idx)];
      M_unit = M./norm(M);
      T = M - M_unit * 500;
      param.target.x  = T(1);
      param.target.y  = T(2);
      param.target.z  = T(3);
      [param.target.lat, param.target.lon, param.target.elev] = ecef2geodeticD(param.target.x, param.target.y, param.target.z, WGS84.ellipsoid);
      
    else
      % target is at same lat lon as mid_idx but at zero elevation
      param.target.lat  = param.gps.lat(mid_idx);
      param.target.lon  = param.gps.lon(mid_idx);
      if isnan(param.gps.elev(mid_idx))
        param.target.elev = 0;
      elseif 1 % temporary override
        param.target.elev = 0;
      end
      [param.target.x, param.target.y, param.target.z] = geodeticD2ecef(param.target.lat, param.target.lon, param.target.elev, WGS84.ellipsoid);
      
    end
    
end

if 1
  % check 3d plot
  figure;
  plot3(param.gps.x, param.gps.y, param.gps.z, '.-'); hold on;
  grid on;
  plot3(param.target.x, param.target.y, param.target.z,'o');
  
  T = [param.target.x; param.target.y; param.target.z];
  M = [param.gps.x(mid_idx); param.gps.y(mid_idx); param.gps.z(mid_idx)];
  A = [param.gps.x(1); param.gps.y(1); param.gps.z(1)];
  B = [param.gps.x(end); param.gps.y(end); param.gps.z(end)];
  pos_vect(T, A );
  pos_vect(T, M );
  pos_vect(T, B );
  
end

clear mid_idx T M A B M_unit
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