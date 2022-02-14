% Script sim.input_full_sim
%
% Generates simulated data from an extracted flightline
% Save simulated raw_data, params, records, frames 
%
% Authors: John Paden, Hara Madhav Talasila
%
% See also sim.flightline_extract, run_load_data (example 7)

param=[];
raw_data = [];
frames = [];
records = [];

%% INPUTS to sim.flightline_extract.m
% =========================================================================
% REQUIRED

% param.sim.radar_name  = 'snow';
% param.sim.season_name = '2012_Greenland_P3'; %'2018_Antarctica_DC8';
% param.sim.day_seg     = '20120330_04'; % '20181010_01';
%
% param.sim.radar_name  = 'snow';
% param.sim.season_name = '2013_Greenland_P3';
% param.sim.day_seg     = '20130327_02';
%
% param.radar_name      = 'snow';
% param.season_name     = '2013_Greenland_P3';
% param.day_seg         = '20130327_02';
%
% param.sim.radar_name  = 'snow';
% param.sim.season_name = '2016_Greenland_Polar6';
% param.sim.day_seg     = '20160414_01';
%
% param.sim.radar_name  = 'snow';
% param.sim.season_name = '2018_Antarctica_DC8';
% param.sim.day_seg     = '20181010_01';

param.sim.radar_name  = 'rds';
param.sim.season_name = '2014_Greenland_P3';
param.sim.day_seg     = '20140325_05';

% Images to sim/load
param.sim.imgs = {[1 2]};% ,[2,1] {[ones([7 1]),(1:7).'], [2*ones([7 1]),(1:7).'], [3*ones([7 1]),(1:7).']};

% =========================================================================
% Optional

% Northward flight instead of actual trajectory
param.sim.north_along_track_en = 1;

% =========================================================================

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.sim.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% flightline_extract

[ param, frames, records, exec_good ] = sim.flightline_extract(param);

if ~exec_good; 
  fprintf('flightline_extract executed incompletely\n'); 
  return; 
else
  fprintf('(%s) Flightline extracted\n',  datestr(now));
  fprintf('=====================================================================\n');
end

clear exec_good;

%% Simulator

[c, WGS84] = physical_constants('c', 'WGS84');
debug_plot_en  = 0;
raw_data = [];

Nx = length(param.gps.elev); % rec_len
wfs = param.radar.wfs;

[output_dir,radar_type,radar_name] = ct_output_dir(param.sim.radar_name);

% Calculate data matrix cell(s)
for img = 1:length(param.sim.imgs)
  for wf_adc = 1:size(param.sim.imgs{img},1)
    
    wf = param.sim.imgs{img}(wf_adc,1);
    adc = param.sim.imgs{img}(wf_adc,2);
    
    % Create actual trajectory #### from data_load %% Create traj
    trajectory_param = struct('gps_source',records.gps_source, ...
      'season_name',param.season_name,'radar_name',param.radar_name, ...
      'rx_path', param.radar.wfs(wf).rx_paths(adc), ...
      'tx_weights', param.radar.wfs(wf).tx_weights, 'lever_arm_fh', param.radar.lever_arm_fh);
    [gps,lever_arm_val] = trajectory_with_leverarm(records,trajectory_param);
    
    [gps.x, gps.y, gps.z] = geodeticD2ecef(gps.lat, gps.lon, gps.elev, WGS84.ellipsoid);
    param.traj = gps;
    
    
    time_raw    = wfs(wf).time_raw;
    Tpd         = wfs(wf).Tpd;
    f0          = wfs(wf).f0;
    fc          = wfs(wf).fc;
    chirp_rate  = wfs(wf).chirp_rate;
    rx_tukey    = wfs(wf).tukey;
    fs          = param.radar.fs;
    
    % point target for now
    range = sqrt( (param.traj.x - param.target.x).^2 + (param.traj.y - param.target.y).^2 + (param.traj.z - param.target.z).^2 );
    twtt = 2*range/c;
    
    Nt = wfs(wf).Nt_raw;
    
    for rec = 1:Nx
      
      td    = twtt(rec);
      fb    = chirp_rate*td;
      time  = time_raw - td; % ################
      
      if strcmpi(radar_type,'deramp') % deramp
        
        raw_data{wf,adc}(:,rec) = tukeywin_cont(time/Tpd-0.5,rx_tukey) ...
          .* cos(2*pi*f0*td + pi*chirp_rate*(2*time_raw*td - td^2)); % Nt x Nrx;
        
      else % pulsed
        
        raw_data{img}(:,rec,wf_adc) = tukeywin_cont(time/Tpd-0.5,rx_tukey) ...
          * 0.5 .* exp(1i*(2*pi*f0*time + pi*chirp_rate*time.^2 )); % Nt x Nrx;
        
      end
      
    end %% for rec
    
  end %% for wf_adc
end %% for img

param.sim.range = range;
param.sim.twtt  = twtt;

param.records.file.version = 1000;

%% Output_Files

false_save_en = 0; % if 1; does not save the files

% RECORDS
param.sim.out_fn_dir_records = fullfile(gRadar.support_path,'records', ...
  sprintf('%s',param.sim.radar_name), ...
  sprintf('%s',param.season_name) );
param.sim.out_fn_records = fullfile(param.sim.out_fn_dir_records, ...
  sprintf('records_%s.mat',param.day_seg) );
fprintf('Saving records %s (%s)\n', param.sim.out_fn_records, datestr(now));
if ~false_save_en; ct_save(param.sim.out_fn_records,'-struct', 'records');
  if 0; fprintf('Creating auxiliary records files %s (%s)\n',param.sim.out_fn_records,datestr(now));
    create_records_aux_files(param.sim.out_fn_records);end; end;

% FRAMES
param.sim.out_fn_dir_frames = fullfile(gRadar.support_path,'frames', ...
  sprintf('%s',param.sim.radar_name), ...
  sprintf('%s',param.season_name) );
param.sim.out_fn_frames = fullfile(param.sim.out_fn_dir_frames, ...
  sprintf('frames_%s.mat',param.day_seg) );
fprintf('Saving frames %s (%s)\n', param.sim.out_fn_frames, datestr(now));
if ~false_save_en; ct_save(param.sim.out_fn_frames, 'frames'); end;

% RAW DATA
param.sim.out_fn_dir_raw_data = fullfile(gRadar.ct_tmp_path,'sim3D', ...
  sprintf('%s',param.sim.radar_name), ...
  sprintf('%s',param.season_name), ...
  sprintf('%s',param.day_seg(1:8)) );
for img = 1:length(param.sim.imgs)
  for wf_adc = 1:size(param.sim.imgs{img},1)
    wf = param.sim.imgs{img}(wf_adc,1);
    adc = param.sim.imgs{img}(wf_adc,2);
    param.sim.out_fn_raw_data = fullfile(param.sim.out_fn_dir_raw_data, ...
      sprintf('data_wfs_%02d_adc_%02d.mat',wf,adc) );
    fprintf('Saving raw_data %s (%s)\n', param.sim.out_fn_raw_data, datestr(now));
    if ~false_save_en; ct_save(param.sim.out_fn_raw_data, 'raw_data');end; 
    % 'hdr', 'dt', 't0', 'IF_filter_idx', 'NCO', 'DDC_filter_idx', 'dec'
  end
end

% PARAM
param.sim.out_fn_dir_param = param.sim.out_fn_dir_raw_data;
param.sim.out_fn_param = fullfile(param.sim.out_fn_dir_param, ...
  sprintf('param.mat') );
param.fn = param.sim.out_fn_param;
fprintf('Saving param %s (%s)\n', param.sim.out_fn_param, datestr(now));
if ~false_save_en; ct_save(param.fn, 'param');end;

if false_save_en; 
  fprintf('===== FALSE_SAVE_EN: none of the above files are written =====\n');
else
  fprintf('Files saved (%s)\n',  datestr(now));
  fprintf('=====================================================================\n');
end;

%% FIGURES (work best for single target case)

wf = 1;
t_ref = param.radar.wfs(wf).time_raw;
raw_data = raw_data{1};
tmp = 20*log10(abs(raw_data)) ;
y = 1:Nx;

call_sign = 'Simulated Data';

for compressing_this=1 %continue;
  
  fig_title = sprintf('%s_%s',mfilename, call_sign);
  fig_h = figure('Name',fig_title);
  
  subplot(121);
  imagesc(y, t_ref/1e-6, tmp-max(tmp(:)) ,[-30,0] );
  cb = colorbar; cb.Label.String = 'Relative Power, dB';
  grid on; hold on; axis tight;
  plot(y, twtt/1e-6, '--', 'Color', 'g');
  [td_closest, idx_closest] = min(twtt);
  plot(idx_closest, twtt(idx_closest)/1e-6, 's','LineWidth',2);
  xlabel('Along-track position, rlines'); ylabel('Fast-time, us');
  title('PhaseHistory Magnitude');
  
  subplot(122);
  imagesc(y, t_ref/1e-6, angle(raw_data) );
  cb = colorbar; cb.Label.String = 'Radians';
  grid on; hold on; axis tight;
  plot(y, twtt/1e-6, '--', 'Color', 'g');
  xlabel('Along-track position, rlines'); ylabel('Fast-time, us');
  title('PhaseHistory Phase');
  
  try
    sgtitle(fig_title,'FontWeight','Bold','FontSize',14,'Interpreter','None');
  end
  set(findobj(fig_h,'type','axes'),'FontWeight', 'Bold', 'FontSize',14);
  set(fig_h, 'Position', get(0, 'Screensize'));
  %   print(gcf, '-dpng', fig_title, '-r300');
  
  % Frequency check
  fig_title = 'Frequency check';
  fig_h = figure('Name', fig_title);
  hold on;
  
  f0    = param.radar.wfs(wf).f0;
  f1    = param.radar.wfs(wf).f1;
  fs    = param.radar.fs;
  IF_nz = [0:3];
  IF_cutoffs = [0 : 0.5 : 2]' * fs ; % Hz
  
  for idx = 1:length(IF_cutoffs)-1
    rectangle('Position',[IF_cutoffs(idx)/1e6 0 fs/2e6 1.5],'EdgeColor','k');
  end
  rectangle('Position',[f0/1e6 0 (f1-f0)/1e6 1],'EdgeColor','r','FaceColor', 'g');
  xticks( round(sort( [IF_cutoffs; f0; f1] ./1e6 ),2) );
  yticks([1 1.5]);
  grid on; axis tight;
  ylim([0 2]);

  xlabel('Frequency, MHz'); ylabel('Frequency bands');
  title('4 Nyquist zones and passband');

end
