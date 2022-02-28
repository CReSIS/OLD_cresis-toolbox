function [param, records, records2, frames, exec_good] = flightline_extract(param)

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

%% Load records, frames

try
  fprintf('(%s) Loading records', datestr(now));
  records = records_load( param );
catch ME
  fprintf('\t\t-- Failed (%s)\n',  datestr(now)); return;
end
fprintf('\t\t-- Loaded (%s)\n',  datestr(now));

try
  fprintf('(%s) Loading frames',  datestr(now));
  frames = frames_load( param );
catch ME
  fprintf('\t\t-- Failed (%s)\n',  datestr(now)); return;
end
fprintf('\t\t-- Loaded (%s)\n',  datestr(now));

%% Setup records, frames

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

% override default variables in frames
frames.frame_idxs = 1;
frames.nyquist_zone = NaN; %1
frames.proc_mode = 0;
frames.quality = 1;

fprintf('(%s) Setup records, frames \t-- Done\n',  datestr(now));

%% Waveforms

param.load.imgs = param.sim.imgs;
[wfs,~] = data_load_wfs(param,records);

% overrides for each waveform
for wf = 1:length(wfs)
  
  if 0 % Set deramp flag
    [output_dir,radar_type,radar_name] = ct_output_dir(param.sim.radar_name);
    if strcmpi(radar_type,'deramp')
      wfs.deramp(wfs) = 1; % Default 0
    elseif strcmpi(radar_type,'pulsed')
      wfs.deramp(wf) = 0; % Default 0
    end
  end
  
  % Assign param.radar.wfs(wf).tx_weights = ones to use in simulator
  wfs(wf).tx_weights = ones(size(wfs(wf).tx_weights));
  
end

param.radar.wfs = merge_structs(param.radar.wfs,wfs);

if param.sim.debug_plots_en
  figure;
  for wf = 1:length(wfs)
    tmp = param.radar.wfs(wf).time_raw;
    plot(tmp/1e-6, wf*ones(size(tmp)), '.'); hold on;
  end
  grid on;
  xlabel('raw time, us');
  ylabel('wf');
  yticks(1:length(wfs));
  ylim([0 length(wfs)+1]);
end

clear wf wfs tmp;

fprintf('(%s) Waveforms \t\t-- Done\n',  datestr(now));

%% Signal % not used for now

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

clear Ntx Nrx idx;

fprintf('(%s) Signal \t\t\t-- Done\n',  datestr(now));

%%  Trajectory and Attitude (GPS+IMU)
% Creates param.gps
% Overwrites records traj with a northward traj if enabled

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
  
  if param.sim.debug_plots_en
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
    
    clear traj_ecef_diff test_unit_vec
  end
  
  param.gps = rec;

  % Overwrite simulated rajectory to records file
  records = merge_structs(records, param.gps);
  
  clear along_track dx a A b B north_unit_vec traj_ecef rec
  
  fprintf('(%s) Flightline Northward \t-- Done\n',  datestr(now));
  
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
  
  fprintf('(%s) Flightline RefTraj \t-- Done\n',  datestr(now));
  
end

%% Create phase center trajectory

% Create reference trajectory (rx_path == 0, tx_weights = [])
% tx_weights = []  or param.radar.wfs(wf).tx_weights doesn not matter
% because they will be set to ones if rx_path == 0;
trajectory_param = struct('gps_source',records.gps_source, ...
  'season_name',param.season_name,'radar_name',param.radar_name, ...
  'rx_path', 0, ...
  'tx_weights', [], ... 
  'lever_arm_fh', param.radar.lever_arm_fh);
[records2, lever_arm_val] = trajectory_with_leverarm(records,trajectory_param);

[records2.x, records2.y, records2.z] = geodeticD2ecef(records2.lat, records2.lon, records2.elev, WGS84.ellipsoid);

if param.sim.debug_plots_en
  figure; 
  plot(records.lon, records.lat, 'x'); hold on;
  plot(records2.lon, records2.lat, 'o'); 
end

clear trajectory_param

fprintf('(%s) Reference Trajectory \t-- Done\n',  datestr(now));

%% Target(s) location

% using records structure for ease. Same geodetic, ecef as param.gps

param.target = [];

if ~isfield(param.sim,'target') || ~isfield(param.sim.target,'type') ...
    || strcmpi(param.sim.target.type, 'point')
  param.target.type = 'point'; % 'surface' %'layer'
end

switch param.target.type
  case 'point'
        
    if isfield(param.sim, 'north_along_track_en') && param.sim.north_along_track_en
      % For simple target position method:
      %
      altra = [];
      % 1. Compute along-track vector (geodetic_to_along_track)
      %
      [altra.along_track, altra.lat, altra.lon, altra.elev] = ...
        geodetic_to_along_track(records2.lat, records2.lon, records2.elev);
      
      % 2. Find ECEF of trajectory
      %
      [altra.x, altra.y, altra.z] = geodeticD2ecef(altra.lat, altra.lon, altra.elev, WGS84.ellipsoid);
      
      % 3. For each point in the trajectory, find the XYZ unit vectors in the
      % flight coordinate system.
      %
      % X, Y, Z are 3x(rec_len) matrices of (rec_len) unit vectors
      X = [];
      U = [];
      N = [];
      Z = [];
      Y = [];
      
      % X is along track
      X = [diff(altra.x); diff(altra.y); diff(altra.z)]; % 3x(rec_len-1)
      X = X./vecnorm(X);
      % extend last unit vector to the last point % 3x(rec_len)
      X = [X, X(:,end)]; 
      
      
      % U is up vector
      [U(1,:), U(2,:), U(3,:)] = enu2ecef( 0, 0, 1, altra.lat, altra.lon, altra.elev, WGS84.ellipsoid);
      U = U - [altra.x; altra.y; altra.z];
      U = U./vecnorm(U); 
      
      % N is normal to the plane with X, U and Zenith
      N = cross(X, U, 1);
      N = N./vecnorm(N);
      
      % Z is Zenith
      Z = cross(N, X);
      Z = Z./vecnorm(Z);
      
      % Y is left ==> X Y Z are right-handed coord system
      Y = cross(Z, X);
      
      if param.sim.debug_plots_en
        figure;
        plot(vecnorm(Y),'x'); hold on;
      end
      
      Y = Y./vecnorm(Y);
      
      if param.sim.debug_plots_en
        plot(vecnorm(Y),'o');
      end
      
      altra.X = X;
      altra.Y = Y;
      altra.Z = Z;
        
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
      
      % A to B is the flightline
      % C is a point on flightline closest to Target T
      % T is defined using x,y,z wrt A
      % x = alongtrack distance from A to C
      % y = horizontal offset from C (zero is directly under flightpath)
      % z = vertical offset from C (height)
      
      A = [records2.x(1); records2.y(1); records2.z(1)];
      B = [records2.x(end); records2.y(end); records2.z(end)];
      Lsar = norm(B-A);
      
      % support user input
      x = Lsar/2;  % for any distance <Lsar
      y = 0; % positive(left to flighline) negative(left to flightline)
      z = -500; % positive(above flightline) negative(below flightline)
      z = -base2dec('KU', 36); % about  5 us TWTT for 750 meter
      z = -base2dec('157', 36); % about 9.89 us TWTT for 1483 meter
      
      idx_C = find(altra.along_track >= x, 1);
      C = [altra.x(idx_C); altra.y(idx_C); altra.z(idx_C)];
      
      T = C ...
        + 0 * altra.X(:,idx_C) ... % zero x offset from C
        + y * altra.Y(:,idx_C) ... % y offset from C
        + z * altra.Z(:,idx_C) ; % z offset from C

      param.target.x  = T(1);
      param.target.y  = T(2);
      param.target.z  = T(3);
      [param.target.lat, param.target.lon, param.target.elev] = ecef2geodeticD(param.target.x, param.target.y, param.target.z, WGS84.ellipsoid);
            
      if param.sim.debug_plots_en
        % check 3d plot for FCS XYZ
        figure;
        plot3(records2.x, records2.y, records2.z, '.'); hold on;
        plot3(records2.x(1), records2.y(1), records2.z(1), 'x','Color','r');
        grid on;
        unit_vect([records2.x; records2.y; records2.z], X, 50);
        unit_vect([records2.x; records2.y; records2.z], Y, 50);
        unit_vect([records2.x; records2.y; records2.z], Z, 50);
        pos_vect(T, A );
        pos_vect(T, C );
        pos_vect(T, B );
        xlabel('x ECEF');
        ylabel('y ECEF');
        zlabel('z ECEF');
        
        figure;
        test_range = vecnorm([records2.x; records2.y; records2.z]-T);
        plot(test_range,'x');
      end
      
      clear altra X U N Z Y A B Lsar x y z idx_C C T
      
    elseif 0
      % M is midpoint of northward track and MT is towards the target (T).
      % T is 500 meters from M towards earth's center
      mid_idx = floor(rec_len/2);
      M = [records2.x(mid_idx); records2.y(mid_idx); records2.z(mid_idx)];
      M_unit = M./norm(M);
      T = M - M_unit * 500;
      param.target.x  = T(1);
      param.target.y  = T(2);
      param.target.z  = T(3);
      [param.target.lat, param.target.lon, param.target.elev] = ecef2geodeticD(param.target.x, param.target.y, param.target.z, WGS84.ellipsoid);
      
      if param.sim.debug_plots_en
        % check 3d plot
        figure;
        plot3(records2.x, records2.y, records2.z, '.-'); hold on;
        plot3(param.target.x, param.target.y, param.target.z,'o');
        grid on;
        T = [param.target.x; param.target.y; param.target.z];
        M = [records2.x(mid_idx); records2.y(mid_idx); records2.z(mid_idx)];
        A = [records2.x(1); records2.y(1); records2.z(1)];
        B = [records2.x(end); records2.y(end); records2.z(end)];
        pos_vect(T, A );
        pos_vect(T, M );
        pos_vect(T, B );
      end
      clear mid_idx T M A B M_unit
      
    elseif 0
      % target is at same lat lon as mid_idx but at zero elevation
      mid_idx = floor(rec_len/2);
      param.target.lat  = records2.lat(mid_idx);
      param.target.lon  = records2.lon(mid_idx);
      if isnan(records2.elev(mid_idx))
        param.target.elev = 0;
      elseif 1 % temporary override
        param.target.elev = 0;
      end
      [param.target.x, param.target.y, param.target.z] = geodeticD2ecef(param.target.lat, param.target.lon, param.target.elev, WGS84.ellipsoid);
      clear mid_idx
      
    end
    
end % switch param.target.type

fprintf('(%s) Target \t\t\t-- Done\n',  datestr(now));

%% Closing

param.season_name = sprintf('%ssim',param.season_name);

% records file version
param.records.file.version_before_FullSim = param.records.file.version;
param.records.file.version = 1000;

% To use while loading simulated data
param.load_data.recs = [start_idx stop_idx];
param.load_data.imgs = param.sim.imgs;

exec_good = 1;

return;


end

% extra functions used here but not in cresis-toolbox/utility/
%% pos_vect

function pos_vect(A,varargin)

% function pos_vect(A,varargin)
% pos_vect(A) draws a line from origin to A
% pos_vect(A,varargin) draws a line from A to B
% plots distance near the midpoint
% Author: Hara Madhav Talasila

if nargin==1
  line([0;A(1)], [0;A(2)], [0;A(3)],'LineStyle','-','LineWidth',1,'Color',[152,251,152]/256);
  p2p_dist = norm(A);
  midP = A/2;
elseif nargin==2
  B = varargin{1};
  line([A(1);B(1)], [A(2);B(2)], [A(3);B(3)],'LineStyle','-','LineWidth',1,'Color',[152,251,152]/256);
  p2p_dist = norm(A-B);
  midP = (A+B)/2;
end

text(midP(1), midP(2), midP(3), sprintf('%f',p2p_dist));

end

%% unit_vect

function unit_vect(A, cap, dist)

% function unit_vect(A, cap, dist)
% unit_vect(A, cap, dist) draws a line from A to B
% B here is at a distance,'dist', in 'cap' direction from A
% Useful to visualize continuous local coordinate system
% Author: Hara Madhav Talasila

B = A + cap*dist;

for idx = 1:size(A,2)
  line([A(1,idx);B(1,idx)], [A(2,idx);B(2,idx)], [A(3,idx);B(3,idx)],'LineStyle','-.','LineWidth',0.5,'Color',[152,251,152]/256);
end

end
