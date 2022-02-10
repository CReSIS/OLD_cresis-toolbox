try hm; end;

param=[];
data = [];
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
  fprintf('Flightline extracted (%s)\n',  datestr(now));
  fprintf('=====================================================================\n');
end

clear exec_good;

%% Simulator

c = physical_constants('c');
debug_plot_en  = 0;
data = [];

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
    
    time_raw    = wfs(wf).time_raw;
    Tpd         = wfs(wf).Tpd;
    f0          = wfs(wf).f0;
    fc          = wfs(wf).fc;
    chirp_rate  = wfs(wf).chirp_rate;
    rx_tukey    = wfs(wf).tukey;
    fs          = param.radar.fs;
    
    % point target for now
    range = sqrt( (param.gps.x - param.target.x).^2 + (param.gps.y - param.target.y).^2 + (param.gps.z - param.target.z).^2 );
    twtt = 2*range/c;
    
    Nt = wfs(wf).Nt_raw;
    
    for rec = 1:Nx
      
      td    = twtt(rec);
      fb    = chirp_rate*td;
      time  = time_raw - td;
      
      if strcmpi(radar_type,'deramp') % deramp
        
        data{wf,adc}(:,rec) = tukeywin_cont(time/Tpd-0.5,rx_tukey) .* cos(2*pi*f0*td + pi*chirp_rate*(2*time_raw*td - td^2)); % Nt x Nrx;
        
      else % pulsed
        
        data{img}(:,rec,wf_adc) = tukeywin_cont(time/Tpd-0.5,rx_tukey) * 0.5 .* exp(-1i*2*pi*f0*td + 1i*pi*chirp_rate*time.^2); % Nt x Nrx;
        
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
    raw_data = data;
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

%% PLOTS

fig_h = figure(123);clf(123);
subplot(221)
% try
imagesc(1:Nx,param.radar.wfs(1).time/1e-6,lp(data{1}));
% catch
%   imagesc(1:Nx,param.radar.wfs(1).time/1e-6,lp(data));
% end
title(param.sim.season_name,'Interpreter', 'none');
% hold on;
% plot(param.target.elev*2/3e8/1e-6,'x','LineWidth',2);
xlabel('rlines');
ylabel('fast-time, us');
subplot(222)
try
  plot(param.radar.wfs(1).time/1e-6,lp(data{1,1}(:,ceil(Nx/2))));
catch
  %   plot(param.radar.wfs(1).time/1e-6,lp(data(:,ceil(Nx/2))));
end
xlabel('fast-time, us');
ylabel('center rline');
title(param.sim.day_seg,'Interpreter', 'none');
axis tight;
subplot(223)
plot(1:Nx,range,'x');
xlabel('rlines');
ylabel('range, m');
axis tight;
subplot(224)
plot(1:Nx,twtt/1e-6,'x');
xlabel('rlines');
ylabel('twtt, us');
axis tight;

% figure(369);
% plot(1:Nx,twtt/1e-6);
% xlabel('rlines');
% ylabel('twtt, us');
% axis tight;

% fig_title = sprintf('%s %s',param.sim.season_name, param.sim.day_seg);
% print(fig_h, '-dpng', fig_title, '-r300');
% saveas(fig_h,fig_title);

% figure(124);clf(124);
% imagesc(1:Nx,param.radar.wfs.time/1e-6,lp(data_test{1,1}));
% hold on;
% plot(param.target.elev*2/3e8/1e-6,'x','LineWidth',2);
% xlabel('rlines');
% ylabel('fast-time, us');

