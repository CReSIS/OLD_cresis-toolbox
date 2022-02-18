% Script sim.input_full_sim
%
% Generates simulated data from an extracted flightline
% Save simulated raw_data, params, records, frames
%
% Authors: John Paden, Hara Madhav Talasila
%
% See also sim.flightline_extract, run_load_data (example 7)

try; hara; end;

param=[];
data = [];
frames = [];
records = [];

%% INPUTS to sim.flightline_extract.m
% =========================================================================
% REQUIRED
% Images to sim/load  {[ones([7 1]),(1:7).'], [2*ones([7 1]),(1:7).'], [3*ones([7 1]),(1:7).']};

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

% param.sim.radar_name  = 'rds';
% param.sim.season_name = '2014_Greenland_P3';
% param.sim.day_seg     = '20140325_05';
% param.sim.imgs        = {[1 5]};
% param.sim.imgs        = {[1 5], [2 5] , [3 5]};

param.sim.radar_name  = 'rds';
param.sim.season_name = '2014_Greenland_P3';
param.sim.day_seg     = '20140410_01';
% param.sim.frame_idx   = 57;
param.sim.imgs        = {[1 5]};
param.sim.imgs        = {[1 5], [2 5] , [3 5]};

% param.sim.radar_name  = 'rds';
% param.sim.season_name = '2018_Greenland_P3';
% param.sim.day_seg     = '20180429_01';
% param.sim.imgs        = {[1 9]};
% param.sim.imgs        = {[1 9], [3 9], [5 9]};

% =========================================================================
% Optional

% Northward flight instead of actual trajectory
param.sim.north_along_track_en = 1;

% =========================================================================

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.sim.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% flightline_extract

[ param, records, frames, exec_good ] = sim.flightline_extract(param);

if ~exec_good;
  fprintf('flightline_extract executed incompletely\n');
  return;
else
  fprintf('(%s) FullSim Flightline EXTRACTED\n',  datestr(now));
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

if 1
  figure(369);
  plot(param.gps.lon, param.gps.lat, '.', 'LineWidth', 2);
  hold on; grid on;
  xlabel('Longitude, Degrees');
  ylabel('Latitude, Degrees');
  leg_str= {'Reference'};
  legend(leg_str);
  title('Trajectories');
  line_markers = ['x', 'd', 'h'];
end

% Calculate raw_data matrix cell(s)
for img = 1:length(param.sim.imgs)
  for wf_adc = 1:size(param.sim.imgs{img},1)
    
    wf = param.sim.imgs{img}(wf_adc,1);
    adc = param.sim.imgs{img}(wf_adc,2);
    
    % Create actual trajectory #### from data_load %% Create traj
    trajectory_param = struct('gps_source',records.gps_source, ...
      'season_name',param.season_name,'radar_name',param.radar_name, ...
      'rx_path', param.radar.wfs(wf).rx_paths(adc), ...
      'tx_weights', param.radar.wfs(wf).tx_weights, ...
      'lever_arm_fh', param.radar.lever_arm_fh);
    [gps,lever_arm_val] = trajectory_with_leverarm(records,trajectory_param);
    
    [gps.x, gps.y, gps.z] = geodeticD2ecef(gps.lat, gps.lon, gps.elev, WGS84.ellipsoid);
    
    if 0
      sum(abs(records.x - gps.x))
      sum(abs(records.y - gps.y))
      sum(abs(records.z - gps.z))
      sum(abs(records.lat - gps.lat))
      sum(abs(records.lon - gps.lon))
      sum(abs(records.elev - gps.elev))
    end
    
    if 0
      % override sim trajectory with gps from extracted flightline
      gps = param.gps;
    end
    
    if 1
      figure(369);
      plot(gps.lon, gps.lat, line_markers(wf), 'LineWidth', img);
      leg_str = [leg_str {sprintf('wfs %02d adc %02d',wf,adc)}];
      legend(leg_str);
    end
    
    time_raw    = wfs(wf).time_raw;
    Tpd         = wfs(wf).Tpd;
    f0          = wfs(wf).f0;
    fc          = wfs(wf).fc;
    chirp_rate  = wfs(wf).chirp_rate;
    rx_tukey    = wfs(wf).tukey;
    fs          = param.radar.fs;
    
    % point target for now
    range = sqrt( (gps.x - param.target.x).^2 + (gps.y - param.target.y).^2 + (gps.z - param.target.z).^2 );
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
    
    param.traj{img,wf_adc}      = gps;
    param.sim.range{img,wf_adc} = range;
    param.sim.twtt{img,wf_adc}  = twtt;
    
  end %% for wf_adc
end %% for img

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
param.sim.out_fn_raw_data = fullfile(param.sim.out_fn_dir_raw_data, ...
  sprintf('data_%s.mat', param.day_seg) );
fprintf('Saving raw_data %s (%s)\n', param.sim.out_fn_raw_data, datestr(now));
if ~false_save_en; ct_save(param.sim.out_fn_raw_data, 'raw_data');end;

% PARAM
param.sim.out_fn_dir_param = param.sim.out_fn_dir_raw_data;
param.sim.out_fn_param = fullfile(param.sim.out_fn_dir_param, ...
  sprintf('param.mat') );
param.fn = param.sim.out_fn_param;
fprintf('Saving param %s (%s)\n', param.sim.out_fn_param, datestr(now));
if ~false_save_en; ct_save(param.fn, 'param');end;

if false_save_en
  fprintf('===== FALSE_SAVE_EN: none of the above files are written =====\n');
else
  fprintf('Files saved (%s)\n',  datestr(now));
  fprintf('=====================================================================\n');
end

%% FIGURES (work best for single target case)

for compressing_this=1 %continue;
  
  for img = 1:length(param.sim.imgs)
    for wf_adc = 1:size(param.sim.imgs{img},1)
      wf = param.sim.imgs{img}(wf_adc,1);
      adc = param.sim.imgs{img}(wf_adc,2);
      
      data = raw_data{img}(:,:,wf_adc);
      
      t_ref = param.radar.wfs(wf).time_raw;
      tmp = 20*log10(abs(data)) ;
      x = 1:Nx;
      [td_closest, idx_closest] = min(twtt);
      
      call_sign = sprintf('Simulated Data wfs_%02d_adc_%02d',wf,adc);
      fig_title = sprintf('%s_%s',mfilename, call_sign);
      fig_h = figure('Name',fig_title);
      
      subplot(121);
      imagesc(x, t_ref/1e-6, tmp-max(tmp(:)) ,[-30,0] );
      cb = colorbar; cb.Label.String = 'Relative Power, dB';
      grid on; hold on; axis tight;
      plot(x, twtt/1e-6, '--', 'Color', 'g');
      plot(idx_closest, twtt(idx_closest)/1e-6, 's','LineWidth',2);
      xlabel('Along-track position, rlines'); ylabel('Fast-time, us');
      title('PhaseHistory Magnitude');
      
      subplot(122);
      imagesc(x, t_ref/1e-6, angle(data) );
      cb = colorbar; cb.Label.String = 'Radians';
      grid on; hold on; axis tight;
      plot(x, twtt/1e-6, '--', 'Color', 'g');
      plot(idx_closest, twtt(idx_closest)/1e-6, 's','LineWidth',2);
      xlabel('Along-track position, rlines'); ylabel('Fast-time, us');
      title('PhaseHistory Phase');
      
      try
        sgtitle(fig_title,'FontWeight','Bold','FontSize',14,'Interpreter','None');
      end
      set(findobj(fig_h,'type','axes'),'FontWeight', 'Bold', 'FontSize',14);
      set(fig_h, 'Position', get(0, 'Screensize'));
      %   print(gcf, '-dpng', fig_title, '-r300');
      
      if 0
        % Frequency check
        fig_title = sprintf('Frequency check wfs_%02d_adc_%02d',wf,adc);
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
      
    end %% for wf_adc
  end %% for img
  
end
