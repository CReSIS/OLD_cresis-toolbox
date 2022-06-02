% script create_settings_2022_Greenland_Polar5
%
% Creates NI radar depth sounder settings
%
% See link_budgets.xls

physical_constants; % c = speed of light

% Define waveforms
if ispc
  base_dir = 'C:\waveforms\';
  arena_base_dir = 'C:\Users\Administrator\Desktop\Arena_Shared\configs\';
  arena_base_dir = 'C:\arena_waveforms\';
else
  base_dir = '~/waveforms/';
  arena_base_dir = '~/rss_waveforms/';
end

% Specify a list of start/stop frequencies
if 1
  f0_list = [180e6];
  f1_list = [210e6];
else
  f0_list = [180e6 150e6];
  f1_list = [210e6 520e6];
end
%Cent.Freq: 335   195   335   300   525   195
DDC_select_list = [2]; % Which DDC mode to use

% cal_settings = [1 2 3 4 5 6];
cal_settings = [2];

prf = 10000;

% presums: Ensure lambda/4 sampling (fudge factor allows difference) and an
%   that presums are an even number.
velocity = [100]; % m/s 
presums = round(c./abs(f1_list+f0_list)/2 ./ velocity * prf / 2)*2

final_DDS_phase = [];
final_DDS_phase_no_time = [];
final_DDS_amp = [];
final_DDS_time = [];
if 0
  % Initial conditions (usually all zeros phase/time with max amplitude)
  for idx = 1:length(f0_list)
    final_DDS_phase{idx} = [0 0 0 0 0 0 0 0];
    final_DDS_phase_no_time = [0 0 0 0 0 0 0 0]; % not used usually
    final_DDS_amp{idx} = [4000 4000 4000 4000 4000 4000 4000 4000];
    final_DDS_time{idx} =  [0 0 2.5 2.5 3.125 3.125 0 0];
  end
else
  % COPY AND PASTE RESULTS FROM basic_tx_chan_equalization_SEASON_NAME.m
  % HERE:
  
  % 150-520 MHz
  final_DDS_phase{end+1} = [24.5	73.3	-18.4	0.0	-49.7	14.3	125.2	54.4];
  final_DDS_phase_no_time{end+1} = [0 0 0 0 0 0 0 0]; % not used usually
  final_DDS_amp{end+1} = [1167	3003	2568	4000	3911	2528	2408	1330];
  final_DDS_time{end+1} =  [0.37	0.58	2.15	2.50	2.39	2.42	0.78	0.43];
    
 % final_DDS_phase{end+1} = [70.7	99.2	-116.2	0.0	-13.4	-12.8	132.0	83.6];
 % final_DDS_phase_no_time{end+1} = [0 0 0 0 0 0 0 0]; % not used usually
 % final_DDS_amp{end+1} = [1285	2868	4000	4000	3702	2852	2145	1609];
 % final_DDS_time{end+1} =  [1.6942	1.7961	2.2041	2.5000	2.3630	2.3930	2.0502	1.9164];
  
  %    final_DDS_phase{end+1} = [0 0 -100 0 0 0 0 0];
  %  final_DDS_phase_no_time{end+1} = [0 0 0 0 0 0 0 0]; % not used usually
  %  final_DDS_amp{end+1} = [1140	2209	4000	4000	4000	3755	2248	1130];
  %  final_DDS_time{end+1} =  [2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5];

  % 180-210 MHz
  final_DDS_phase{end+1} = [64.5	61.3	163.6	0.0	-117.5	155.5	59.8	-64.0];
  final_DDS_phase_no_time{end+1} = [0 0 0 0 0 0 0 0]; % not used usually
  final_DDS_amp{end+1} = [930	2541	2931	4000	3432	3162	2166	1066];
  final_DDS_time{end+1} =  [-1.36	-2.65	2.51	2.50	0.42	2.59	-1.01	-2.00];
  
 % final_DDS_phase{end+1} = [142.3	-91.9	227.7	0.0	209.5	-53.3	-230.6	83.6];
 % final_DDS_phase_no_time{end+1} = [0 0 0 0 0 0 0 0]; % not used usually
 % final_DDS_amp{end+1} = [897	1830	4000	4000	3705	2976	1797	936];
 % final_DDS_time{end+1} =  [-0.4500	1.2410	-1.7874	2.5000	-1.2370	-1.8574	-0.3766	0.2870];
 
  % 320-350 MHz
  final_DDS_phase{end+1} = final_DDS_phase{1};
  final_DDS_phase_no_time{end+1} = final_DDS_phase_no_time{1};
  final_DDS_amp{end+1} = final_DDS_amp{1};
  final_DDS_time{end+1} =  final_DDS_time{1};
  
  % 285-315 MHz
  final_DDS_phase{end+1} = [0 0 0 0 0 0 0 0];
  final_DDS_phase_no_time{end+1} = [0 0 0 0 0 0 0 0]; % not used usually
  final_DDS_amp{end+1} = [4000 4000 4000 4000 4000 4000 4000 4000];
  final_DDS_time{end+1} =  [0 0 2.5 2.5 3.125 3.125 0 0];
    
  % 510-540 MHz
  final_DDS_phase{end+1} = [0 0 0 0 0 0 0 0];
  final_DDS_phase_no_time{end+1} = [0 0 0 0 0 0 0 0]; % not used usually
  final_DDS_amp{end+1} = [4000 4000 4000 4000 4000 4000 4000 4000];
  final_DDS_time{end+1} =  [0 0 2.5 2.5 3.125 3.125 0 0];
    
  % 150-240 MHz
  final_DDS_phase{end+1} = final_DDS_phase{2};
  final_DDS_phase_no_time{end+1} = final_DDS_phase_no_time{2};
  final_DDS_amp{end+1} = final_DDS_amp{2};
  final_DDS_time{end+1} =  final_DDS_time{2};


end


% Hwindow_orig: Desired window created during transmit calibration
%  This is used any time a window that is different from that used
%  during calibration is to be used.
Hwindow_orig = chebwin(8,30).';

%% SETUP
% =========================================================================

param = [];
param.season_name = '2022_Greenland_Polar5';
param.radar_name = 'rds';
param.gps_source = 'awi-final';
clear phase_centers;
for tx_chan = 1:8
  % Just enable the current antenna to determine its phase center
  tx_weights = zeros(1,8);
  tx_weights(tx_chan) = 1;
  rxchan = 4; % Fix the receiver (it should not matter which one you choose)
  % Determine phase center for the antenna
  phase_centers(:,tx_chan) = lever_arm(param, tx_weights, rxchan);
end
% Adjust phase centers to the mean phase center position
phase_centers = bsxfun(@minus,phase_centers,mean(phase_centers,2));

% Create waveforms directories if they do not exist
if ~exist(base_dir,'dir')
  mkdir(base_dir);
end
calval_dir = fullfile(base_dir, 'calval');
if ~exist(calval_dir,'dir')
  mkdir(calval_dir);
end

%% AWI MCoRDS Arena Parameters
[param_defaults,defaults] = default_radar_params_2022_Greenland_Polar5_rds;
arena = defaults{1}.arena;

%% Survey Mode + loopback, noise, and deconv modes
% <4000 m thick ice, 1200 +/- 700 ft AGL
ice_thickness = [3000];
for freq_idx = [1]
  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
  param = merge_structs(param,param_defaults);
  param.arena = arena;
  param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 755; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [1 2 3];
  param.tg.altitude_guard = 700*12*2.54/100;
  param.tg.Haltitude = 1200*12*2.54/100;
  param.tg.Hice_thick = ice_thickness(freq_idx);
  param.prf = prf;
  param.presums = [2 4 presums(freq_idx)-6];
  param.wfs(1).atten = 37;
  param.wfs(2).atten = 0;
  param.wfs(3).atten = 0;
  DDS_amp = final_DDS_amp{cal_settings(freq_idx)};
  param.tx_weights = DDS_amp;
  param.tukey = 0.08;
  param.wfs(1).Tpd = 1e-6;
  param.wfs(2).Tpd = 3e-6;
  param.wfs(3).Tpd = 10e-6;
  param.wfs(1).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(2).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(3).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = (param.f0+param.f1)/2;
  [param.wfs(1:3).tx_mask] = deal([0 0 0 0 0 0 0 0]);
  param.fn = fullfile(base_dir,sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));
  write_cresis_xml(param);
  if freq_idx == 1
    % Default Mode
    param.fn = fullfile(base_dir,'default.xml');
    write_cresis_xml(param);
  end
  % Loopback Mode without delay line
  param.tg.staged_recording = false;
  param.tg.altitude_guard = 1000*12*2.54/100;
  param.tg.Haltitude = 0e-6 * c/2;
  param.tg.Hice_thick = 0; % Long enough for 10 us delay line
  param.fn = fullfile(calval_dir,sprintf('survey_%.0f-%.0fMHz_%.0fus_LOOPBACK_NO_DELAY.xml',param.f0/1e6,param.f1/1e6,param.wfs(end).Tpd*1e6));
  write_cresis_xml(param);
  % Loopback Mode (10e-6 delay line)
  param.tg.staged_recording = false;
  param.tg.altitude_guard = 1000*12*2.54/100;
  param.tg.Haltitude = 10e-6 * c/2;
  param.tg.Hice_thick = 0; % Long enough for 10 us delay line
  param.fn = fullfile(calval_dir,sprintf('survey_%.0f-%.0fMHz_%.0fus_LOOPBACK.xml',param.f0/1e6,param.f1/1e6,param.wfs(end).Tpd*1e6));
  write_cresis_xml(param);
  % Deconvolution Mode (for over calm lake or sea ice lead)
  param.wfs(1).atten = 43;
  param.wfs(2).atten = 43;
  param.wfs(3).atten = 43;
  param.tg.staged_recording = false;
  param.tg.altitude_guard = 3000*12*2.54/100;
  param.tg.Haltitude = 4000*12*2.54/100;
  param.tg.Hice_thick = 0 * 12*2.54/100/sqrt(er_ice);
  param.fn = fullfile(calval_dir,sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_DECONV.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6));
  write_cresis_xml(param);
  if freq_idx == 1
    % Noise Mode
    param.tx_weights = [0 0 0 0 0 0 0 0];
    [param.wfs(1:3).tx_mask] = deal([1 1 1 1 1 1 1 1]);
    param.wfs(1).atten = 37;
    param.wfs(2).atten = 0;
    param.wfs(3).atten = 0;
    param.tg.staged_recording = [1 2 3];
    param.tg.altitude_guard = 500*12*2.54/100;
    param.tg.Haltitude = 1400*12*2.54/100;
    param.tg.Hice_thick = 3250;
    param.fn = fullfile(calval_dir,sprintf('survey_%.0f-%.0fMHz_%.0fus_NOISE.xml',param.f0/1e6,param.f1/1e6,param.wfs(end).Tpd*1e6));
    write_cresis_xml(param);
  end
end



%% Equalization (Using Ocean)
% Haltitude +/- 2000 ft
% For lower altitude, increase attenuation
% Use these settings over ocean or sea ice for fast-time equalization,
% transmit equalization, and receiver equalization.
% Creates one waveform for each of N DDS-transmitters plus a combined
% waveform with all transmitters going.
Haltitude = [2500 2500 0 3500 6000];
Tpd_list = [1e-6 1e-6 3e-6 3e-6 3e-6];
attenuation = [62 45 43 61 55];
fn_hint = {'WATER','ICE','NO_DELAY','WATER','WATER'};
for Tpd_idx = 1:length(Tpd_list)
  Tpd = Tpd_list(Tpd_idx);
  for freq_idx = [1]
    param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
  param = merge_structs(param,param_defaults);
  param.arena = arena;
    param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 700; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
    param.DDC_select = DDC_select_list(freq_idx);
    param.max_duty_cycle = 0.12;
    param.create_IQ = false;
    param.tg.staged_recording = false;
    param.tg.altitude_guard = 2000*12*2.54/100;
    param.tg.Haltitude = Haltitude(Tpd_idx)*12*2.54/100;
    param.tg.Hice_thick = 0;
    param.prf = prf;
    param.presums = [10 10 10 10 10 10 10 10 10];
    [param.wfs(1:8).atten] = deal(attenuation(Tpd_idx)-12);
    [param.wfs(9:9).atten] = deal(attenuation(Tpd_idx));
    param.tx_weights = final_DDS_amp{cal_settings(freq_idx)};
    param.tukey = 0.08;
    param.Tpd = Tpd;
    for wf=1:9
      param.wfs(wf).phase = final_DDS_phase{cal_settings(freq_idx)};
    end
    param.delay = final_DDS_time{cal_settings(freq_idx)};
    param.f0 = f0_list(freq_idx);
    param.f1 = f1_list(freq_idx);
    param.DDC_freq = (param.f0+param.f1)/2;
    for wf=1:8
      param.wfs(wf).tx_mask = ones(1,8);
      param.wfs(wf).tx_mask(9-wf) = 0;
    end
    for wf=9:9
      param.wfs(wf).tx_mask = [0 0 0 0 0 0 0 0];
    end
    param.fn = fullfile(calval_dir,sprintf('txequal_%.0f-%.0fMHz_%.0fft_%.0fft_%.0fus_%s.xml', ...
      param.f0/1e6, param.f1/1e6, (param.tg.Haltitude+[-param.tg.altitude_guard param.tg.altitude_guard])*100/12/2.54, ...
      param.Tpd*1e6,fn_hint{Tpd_idx}));
    write_cresis_xml(param);
  end
end


return

ice_thickness = [4000 4000 4000 4000 4000 4000];
for freq_idx = [1 2 3 4 5 6]
  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
  param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 755; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [1 2 3];
  param.tg.altitude_guard = 700*12*2.54/100;
  param.tg.Haltitude = 1200*12*2.54/100;
  param.tg.Hice_thick = ice_thickness(freq_idx);
  param.prf = prf;
  param.presums = [2 4 presums(freq_idx)-6];
  param.wfs(1).atten = 37;
  param.wfs(2).atten = 0;
  param.wfs(3).atten = 0;
  DDS_amp = final_DDS_amp{cal_settings(freq_idx)};
  param.tx_weights = DDS_amp;
  param.tukey = 0.08;
  param.wfs(1).Tpd = 1e-6;
  param.wfs(2).Tpd = 3e-6;
  param.wfs(3).Tpd = 10e-6;
  param.wfs(1).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(2).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(3).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = (param.f0+param.f1)/2;
  [param.wfs(1:3).tx_mask] = deal([0 0 0 0 0 0 0 0]);
  param.fn = fullfile(base_dir,sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));
  write_cresis_xml(param);
  if freq_idx == 1
    % Default Mode
    param.fn = fullfile(base_dir,'default.xml');
    write_cresis_xml(param);
  end
  % Loopback Mode without delay line
  param.tg.staged_recording = false;
  param.tg.altitude_guard = 1000*12*2.54/100;
  param.tg.Haltitude = 0e-6 * c/2;
  param.tg.Hice_thick = 0; % Long enough for 10 us delay line
  param.fn = fullfile(calval_dir,sprintf('survey_%.0f-%.0fMHz_%.0fus_LOOPBACK_NO_DELAY.xml',param.f0/1e6,param.f1/1e6,param.wfs(end).Tpd*1e6));
  write_cresis_xml(param);
  % Loopback Mode (10e-6 delay line)
  param.tg.staged_recording = false;
  param.tg.altitude_guard = 1000*12*2.54/100;
  param.tg.Haltitude = 10e-6 * c/2;
  param.tg.Hice_thick = 0; % Long enough for 10 us delay line
  param.fn = fullfile(calval_dir,sprintf('survey_%.0f-%.0fMHz_%.0fus_LOOPBACK.xml',param.f0/1e6,param.f1/1e6,param.wfs(end).Tpd*1e6));
  write_cresis_xml(param);
  % Deconvolution Mode (for over calm lake or sea ice lead)
  param.wfs(1).atten = 43;
  param.wfs(2).atten = 43;
  param.wfs(3).atten = 43;
  param.tg.staged_recording = false;
  param.tg.altitude_guard = 3000*12*2.54/100;
  param.tg.Haltitude = 4000*12*2.54/100;
  param.tg.Hice_thick = 0 * 12*2.54/100/sqrt(er_ice);
  param.fn = fullfile(calval_dir,sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_DECONV.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6));
  write_cresis_xml(param);
  if freq_idx == 1
    % Noise Mode
    param.tx_weights = [0 0 0 0 0 0 0 0];
    [param.wfs(1:3).tx_mask] = deal([1 1 1 1 1 1 1 1]);
    param.wfs(1).atten = 37;
    param.wfs(2).atten = 0;
    param.wfs(3).atten = 0;
    param.tg.staged_recording = [1 2 3];
    param.tg.altitude_guard = 500*12*2.54/100;
    param.tg.Haltitude = 1400*12*2.54/100;
    param.tg.Hice_thick = 3250;
    param.fn = fullfile(calval_dir,sprintf('survey_%.0f-%.0fMHz_%.0fus_NOISE.xml',param.f0/1e6,param.f1/1e6,param.wfs(end).Tpd*1e6));
    write_cresis_xml(param);
  end
end

%% Survey Mode
% <3250 m thick ice, 1200 +/- 700 ft AGL
ice_thickness = [3250 3250 3250 3250 3250 3250];
for freq_idx = [1 2 3 4 5 6]
  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
  param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 750; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [1 2 3];
  param.tg.altitude_guard = 700*12*2.54/100;
  param.tg.Haltitude = 1200*12*2.54/100;
  param.tg.Hice_thick = ice_thickness(freq_idx);
  param.prf = prf;
  param.presums = [2 4 presums(freq_idx)-6];
  param.wfs(1).atten = 37;
  param.wfs(2).atten = 0;
  param.wfs(3).atten = 0;
  DDS_amp = final_DDS_amp{cal_settings(freq_idx)};
  param.tx_weights = DDS_amp;
  param.tukey = 0.08;
  param.wfs(1).Tpd = 1e-6;
  param.wfs(2).Tpd = 3e-6;
  param.wfs(3).Tpd = 10e-6;
  param.wfs(1).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(2).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(3).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = (param.f0+param.f1)/2;
  [param.wfs(1:3).tx_mask] = deal([0 0 0 0 0 0 0 0]);
  param.fn = fullfile(base_dir,sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));
  write_cresis_xml(param);
end

%% Survey Mode
% <3250 m thick ice, 1200 +/- 700 ft AGL
ice_thickness = [3250 3250 3250 3250 3250 3250];
for freq_idx = [1 2 3 4 5 6]
  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
  param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 750; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [1 2 3];
  param.tg.altitude_guard = 700*12*2.54/100;
  param.tg.Haltitude = 1200*12*2.54/100;
  param.tg.Hice_thick = ice_thickness(freq_idx);
  param.prf = prf;
  param.presums = [2 4 presums(freq_idx)-6];
  param.wfs(1).atten = 37;
  param.wfs(2).atten = 0;
  param.wfs(3).atten = 0;
  DDS_amp = final_DDS_amp{cal_settings(freq_idx)};
  param.tx_weights = DDS_amp;
  param.tukey = 0.08;
  param.wfs(1).Tpd = 1e-6;
  param.wfs(2).Tpd = 3e-6;
  param.wfs(3).Tpd = 10e-6;
  param.wfs(1).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(2).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(3).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = (param.f0+param.f1)/2;
  [param.wfs(1:3).tx_mask] = deal([0 0 0 0 0 0 0 0]);
  param.fn = fullfile(base_dir,sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));
  write_cresis_xml(param);
end

%% Survey Mode
% <4000 m thick ice, 5000 +/- 1400 ft AGL
ice_thickness = [4000 4000 4000 4000 4000 4000];
for freq_idx = [1 2 3 4 5 6]
  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
  param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 750; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [1 2 3];
  param.tg.altitude_guard = 1400*12*2.54/100;
  param.tg.Haltitude = 5000*12*2.54/100;
  param.tg.Hice_thick = ice_thickness(freq_idx);
  param.prf = prf;
  param.presums = [2 4 presums(freq_idx)-6];
  param.wfs(1).atten = 37;
  param.wfs(2).atten = 0;
  param.wfs(3).atten = 0;
  DDS_amp = final_DDS_amp{cal_settings(freq_idx)};
  param.tx_weights = DDS_amp;
  param.tukey = 0.08;
  param.wfs(1).Tpd = 1e-6;
  param.wfs(2).Tpd = 3e-6;
  param.wfs(3).Tpd = 10e-6;
  param.wfs(1).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(2).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(3).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = (param.f0+param.f1)/2;
  [param.wfs(1:3).tx_mask] = deal([0 0 0 0 0 0 0 0]);
  param.fn = fullfile(base_dir,sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));
  write_cresis_xml(param);
end

%% Survey Mode for thin ice
% <2500 m thick ice, 1200 +/- 700 ft AGL
ice_thickness = [2500 2500 2500 2500 2500 2500];
for freq_idx = [1 2 3 4 5 6]
  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
  param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 750; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [1 2 3];
  param.tg.rg_stop_offset = [500*sqrt(3.15) 0 0]; % Keep waveform 1 on for 500 m of ice
  param.tg.altitude_guard = 700*12*2.54/100;
  param.tg.Haltitude = 1200*12*2.54/100;
  param.tg.Hice_thick = ice_thickness(freq_idx);
  param.prf = prf;
  param.presums = [8 2 presums(freq_idx)-10];
  param.wfs(1).atten = 35;
  param.wfs(2).atten = 0;
  param.wfs(3).atten = 0;
  DDS_amp = final_DDS_amp{cal_settings(freq_idx)};
  param.tx_weights = DDS_amp;
  param.tukey = 0.08;
  param.wfs(1).Tpd = 1e-6;
  param.wfs(2).Tpd = 1e-6;
  param.wfs(3).Tpd = 3e-6;
  param.wfs(1).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(2).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(3).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = (param.f0+param.f1)/2;
  [param.wfs(1:3).tx_mask] = deal([0 0 0 0 0 0 0 0]);
  param.fn = fullfile(base_dir,sprintf('thinice_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));
  write_cresis_xml(param);
end

%% Survey Mode for thin ice
% <2500 m thick ice, 1500 +/- 700 ft AGL
ice_thickness = [2500 2500 2500 2500 2500 2500];
for freq_idx = [1 2 3 4 5 6]
  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
  param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 750; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [1 2 3];
  param.tg.rg_stop_offset = [500*sqrt(3.15) 0 0]; % Keep waveform 1 on for 500 m of ice
  param.tg.altitude_guard = 800*12*2.54/100;
  param.tg.Haltitude = 1500*12*2.54/100;
  param.tg.Hice_thick = ice_thickness(freq_idx);
  param.prf = prf;
  param.presums = [8 2 presums(freq_idx)-10];
  param.wfs(1).atten = 35;
  param.wfs(2).atten = 0;
  param.wfs(3).atten = 0;
  DDS_amp = final_DDS_amp{cal_settings(freq_idx)};
  param.tx_weights = DDS_amp;
  param.tukey = 0.08;
  param.wfs(1).Tpd = 1e-6;
  param.wfs(2).Tpd = 1e-6;
  param.wfs(3).Tpd = 3e-6;
  param.wfs(1).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(2).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(3).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = (param.f0+param.f1)/2;
  [param.wfs(1:3).tx_mask] = deal([0 0 0 0 0 0 0 0]);
  param.fn = fullfile(base_dir,sprintf('thinice_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));
  write_cresis_xml(param);
end

%% Sea Ice
% 1200 +/- 1200 ft AGL
for freq_idx = [1]
  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
  param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 750; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [0 0];
  param.tg.rg_stop_offset = [0 0]; % Keep waveform 1 on for 500 m of ice
  param.tg.altitude_guard = 1200*12*2.54/100;
  param.tg.Haltitude = 1200*12*2.54/100;
  param.tg.Hice_thick = 0;
  param.prf = prf;
  param.presums = [round(presums(freq_idx)/2/2)*2 round(presums(freq_idx)/2/2)*2];
  param.wfs(1).atten = 13;
  param.wfs(2).atten = 13;
  DDS_amp = final_DDS_amp{cal_settings(freq_idx)};
  param.tx_weights = DDS_amp;
  param.tukey = 0.08;
  param.wfs(1).Tpd = 1e-6;
  param.wfs(2).Tpd = 1e-6;
  param.wfs(1).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(2).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = (param.f0+param.f1)/2;
  param.wfs(1).tx_mask = deal([0 1 1 1 1 1 1 1]);
  param.wfs(2).tx_mask = deal([1 1 1 1 1 1 1 0]);
  param.fn = fullfile(base_dir,sprintf('seaice_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));
  write_cresis_xml(param);
end

%% Image Mode (EGRIP Low Altitude, Thick Ice)
% Ice thickness "param.tg.Hice_thick_min" m to "param.tg.Hice_thick" m, "param.tg.Haltitude" +/- "param.tg.altitude_guard" ft AGL
for freq_idx = 2
  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
  param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 700; param.flight_hours = 4; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
  param.max_data_rate = 755;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [1 2 3 3];
  param.tg.altitude_guard = 700 * 12*2.54/100;
  param.tg.Haltitude = 1200 * 12*2.54/100;
  param.tg.Hice_thick = 3250;
  param.tg.look_angle_deg = [0 0 40 40];
  param.prf = prf;
  param.presums = [2 4 ceil((presums(freq_idx)-6)/4)*2 ceil((presums(freq_idx)-6)/4)*2];
  % Switch from tx calibration window to hanning window to broaden beam
  DDS_amp = final_DDS_amp{freq_idx} .* hanning(8).' ./ Hwindow_orig;
  % Renormalize the amplitudes
  [~,relative_max_idx] = max(DDS_amp./param.max_tx);
  DDS_amp = round(DDS_amp .* param.max_tx(relative_max_idx) / DDS_amp(relative_max_idx));
  param.tx_weights = DDS_amp;
  param.tukey = 0.08;
  param.wfs(1).Tpd = 1e-6;
  param.wfs(2).Tpd = 3e-6;
  param.wfs(3).Tpd = 10e-6;
  param.wfs(4).Tpd = 10e-6;
  param.wfs(1).phase = final_DDS_phase{freq_idx};
  param.wfs(2).phase = final_DDS_phase{freq_idx};
  param.wfs(3).phase = final_DDS_phase{freq_idx};
  param.wfs(4).phase = final_DDS_phase{freq_idx};
  param.wfs(1).name = 'nadir1';
  param.wfs(2).name = 'nadir3';
  param.wfs(3).name = 'left';
  param.wfs(4).name = 'right';
  % Add in time delays to each position, subtract out the nadir time delays since tx_equalization already took care of those
  beam_angle_deg = 0; % Nadir
  param.wfs(1).delay = final_DDS_time{freq_idx} ...
    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
  beam_angle_deg = 0; % Nadir
  param.wfs(2).delay = final_DDS_time{freq_idx} ...
    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
  beam_angle_deg = 20; % Positive to the left
  param.wfs(3).delay = final_DDS_time{freq_idx} ...
    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
  beam_angle_deg = -20; % Negative to the right
  param.wfs(4).delay = final_DDS_time{freq_idx} ...
    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = (param.f0+param.f1)/2;
  [param.wfs(1:4).tx_mask] = deal([0 0 0 0 0 0 0 0]);
  param.wfs(1).atten = 37;
  param.wfs(2).atten = 0;
  param.wfs(3).atten = 0;
  param.wfs(4).atten = 0;
  param.fn = fullfile(base_dir,sprintf('egrip_image_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml', ...
    param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/2.54/12,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));
  write_cresis_xml(param);
end

%% Image Mode (Low Altitude, Thick Ice)
% Ice thickness "param.tg.Hice_thick_min" m to "param.tg.Hice_thick" m, "param.tg.Haltitude" +/- "param.tg.altitude_guard" ft AGL
for freq_idx = [1 2]
  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
  param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 700; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
  param.max_data_rate = 755;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [0 0 0];
  param.tg.start_ref = {'bottom','surface','bottom'};
  param.tg.stop_ref = {'bottom','surface','bottom'};
  param.tg.altitude_guard = 700 * 12*2.54/100;
  param.tg.Haltitude = 1200 * 12*2.54/100;
  param.tg.Hice_thick_min = 1500;
  param.tg.Hice_thick = 3000;
  param.tg.look_angle_deg = [40 0 40];
  param.prf = prf;
  param.presums = [ceil(presums(freq_idx)/4)*2 2 ceil(presums(freq_idx)/4)*2];
  % Switch from tx calibration window to hanning window to broaden beam
  DDS_amp = final_DDS_amp{freq_idx} .* hanning(8).' ./ Hwindow_orig;
  % Renormalize the amplitudes
  [~,relative_max_idx] = max(DDS_amp./param.max_tx);
  DDS_amp = round(DDS_amp .* param.max_tx(relative_max_idx) / DDS_amp(relative_max_idx));
  param.tx_weights = DDS_amp;
  param.tukey = 0.08;
  param.wfs(1).Tpd = 10e-6;
  param.wfs(2).Tpd = 1e-6;
  param.wfs(3).Tpd = 10e-6;
  param.wfs(1).phase = final_DDS_phase{freq_idx};
  param.wfs(2).phase = final_DDS_phase{freq_idx};
  param.wfs(3).phase = final_DDS_phase{freq_idx};
  param.wfs(1).name = 'left';
  param.wfs(2).name = 'nadir';
  param.wfs(3).name = 'right';
  % Add in time delays to each position, subtract out the nadir time delays since tx_equalization already took care of those
  beam_angle_deg = 20; % Positive to the left
  param.wfs(1).delay = final_DDS_time{freq_idx} ...
    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
  beam_angle_deg = 0; % Nadir
  param.wfs(2).delay = final_DDS_time{freq_idx} ...
    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
  beam_angle_deg = -20; % Negative to the right
  param.wfs(3).delay = final_DDS_time{freq_idx} ...
    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = (param.f0+param.f1)/2;
  [param.wfs(1:3).tx_mask] = deal([0 0 0 0 0 0 0 0]);
  param.wfs(1).atten = 0;
  param.wfs(2).atten = 37;
  param.wfs(3).atten = 0;
  param.fn = fullfile(base_dir,sprintf('image_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml', ...
    param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/2.54/12,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));
  write_cresis_xml(param);
end

%% Image Mode Pattern Measurements
% 3500 ft +/- 1000 ft AGL
for freq_idx = [1 2]
  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
  param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 700; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
  param.max_data_rate = 755;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [0 0];
  param.tg.start_ref = {'surface','surface'};
  param.tg.altitude_guard = 1000 * 12*2.54/100;
  param.tg.Haltitude = 3500 * 12*2.54/100;
  param.tg.Hice_thick_min = 0;
  param.tg.Hice_thick = 0;
  param.tg.look_angle_deg = [0 0 0];
  param.prf = prf;
  param.presums = [ceil(presums(freq_idx)/4)*2 ceil(presums(freq_idx)/4)*2];
  % Switch from tx calibration window to hanning window to broaden beam
  DDS_amp = final_DDS_amp{freq_idx} .* hanning(8).' ./ Hwindow_orig;
  % Renormalize the amplitudes
  [~,relative_max_idx] = max(DDS_amp./param.max_tx);
  DDS_amp = round(DDS_amp .* param.max_tx(relative_max_idx) / DDS_amp(relative_max_idx));
  param.tx_weights = DDS_amp;
  param.tukey = 0.08;
  param.wfs(1).Tpd = 3e-6;
  param.wfs(2).Tpd = 3e-6;
  param.wfs(1).phase = final_DDS_phase{freq_idx};
  param.wfs(2).phase = final_DDS_phase{freq_idx};
  param.wfs(1).name = 'left';
  param.wfs(2).name = 'right';
  % Add in time delays to each position, subtract out the nadir time delays since tx_equalization already took care of those
  beam_angle_deg = 20; % Positive to the left
  param.wfs(1).delay = final_DDS_time{freq_idx} ...
    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
  beam_angle_deg = -20; % Negative to the right
  param.wfs(2).delay = final_DDS_time{freq_idx} ...
    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = (param.f0+param.f1)/2;
  [param.wfs(1:2).tx_mask] = deal([0 0 0 0 0 0 0 0]);
  [param.wfs(1:2).atten] = deal(43);
  param.fn = fullfile(calval_dir,sprintf('image_%.0f-%.0fMHz_%.0fft_%.0fus_PATTERN.xml', ...
    param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/2.54/12,param.wfs(end).Tpd*1e6));
  write_cresis_xml(param);
end

%% Image Mode (High Altitude, Thin Ice)
% Ice thickness "param.tg.Hice_thick_min" m to "param.tg.Hice_thick" m, "param.tg.Haltitude" +/- "param.tg.altitude_guard" ft AGL
freq_idx_WB = 1;
freq_idx_NB = 3;
param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 700; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
param.max_data_rate = 755;
param.DDC_select = DDC_select_list(freq_idx_WB);
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.staged_recording = [0 0 0];
param.tg.start_ref = {'bottom','surface','bottom'};
param.tg.stop_ref = {'bottom','bottom','bottom'};
param.tg.altitude_guard = 1000 * 12*2.54/100;
param.tg.Haltitude = 6000 * 12*2.54/100;
param.tg.Hice_thick_min = 0;
param.tg.Hice_thick = 1000;
param.tg.look_angle_deg = [40 0 40];
param.prf = prf;
param.presums = [ceil(presums(freq_idx_WB)/4)*2 4 ceil(presums(freq_idx_WB)/4)*2];
% Switch from tx calibration window to hanning window to broaden beam
DDS_amp = final_DDS_amp{freq_idx_WB} .* hanning(8).' ./ Hwindow_orig;
% Renormalize the amplitudes
[~,relative_max_idx] = max(DDS_amp./param.max_tx);
DDS_amp = round(DDS_amp .* param.max_tx(relative_max_idx) / DDS_amp(relative_max_idx));
param.tx_weights = DDS_amp;
param.tukey = 0.08;
param.wfs(1).Tpd = 3e-6;
param.wfs(2).Tpd = 1e-6;
param.wfs(3).Tpd = 3e-6;
param.wfs(1).phase = final_DDS_phase{freq_idx_NB};
param.wfs(2).phase = final_DDS_phase{freq_idx_WB};
param.wfs(3).phase = final_DDS_phase{freq_idx_NB};
param.wfs(1).name = 'left';
param.wfs(2).name = 'nadir';
param.wfs(3).name = 'right';
% Add in time delays to each position, subtract out the nadir time delays since tx_equalization already took care of those
beam_angle_deg = 20; % Positive to the left
param.wfs(1).delay = final_DDS_time{freq_idx_NB} ...
  - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
  - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
  + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
beam_angle_deg = 0; % Nadir
param.wfs(2).delay = final_DDS_time{freq_idx_WB} ...
  - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
  - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
  + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
beam_angle_deg = -20; % Negative to the right
param.wfs(3).delay = final_DDS_time{freq_idx_NB} ...
  - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
  - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
  + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
param.wfs(1).f0 = f0_list(freq_idx_NB);
param.wfs(2).f0 = f0_list(freq_idx_WB);
param.wfs(3).f0 = f0_list(freq_idx_NB);
param.wfs(1).f1 = f1_list(freq_idx_NB);
param.wfs(2).f1 = f1_list(freq_idx_WB);
param.wfs(3).f1 = f1_list(freq_idx_NB);
param.DDC_freq = (param.wfs(2).f0+param.wfs(2).f1)/2;
[param.wfs(1:3).tx_mask] = deal([0 0 0 0 0 0 0 0]);
param.wfs(1).atten = 0;
param.wfs(2).atten = 23;
param.wfs(3).atten = 0;
param.fn = fullfile(base_dir,sprintf('imagehighthin_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml', ...
  param.wfs(2).f0/1e6,param.wfs(2).f1/1e6,param.tg.Haltitude*100/2.54/12,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));
write_cresis_xml(param);

%% Image Mode (Low Altitude, Ice <2000 m thick)
% Ice thickness "param.tg.Hice_thick_min" m to "param.tg.Hice_thick" m, "param.tg.Haltitude" +/- "param.tg.altitude_guard" ft AGL
freq_idx_WB = 1;
freq_idx_NB = 2;
param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 700; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
param.max_data_rate = 755;
param.DDC_select = DDC_select_list(freq_idx_NB);
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.staged_recording = [0 0 0];
param.tg.start_ref = {'bottom','surface','bottom'};
param.tg.stop_ref = {'bottom','bottom','bottom'};
param.tg.altitude_guard = 500 * 12*2.54/100;
param.tg.Haltitude = 1200 * 12*2.54/100;
param.tg.Hice_thick_min = 0;
param.tg.Hice_thick = 2000;
param.tg.look_angle_deg = [40 0 40];
param.prf = prf;
param.presums = [8 8 8]; %[ceil(presums(freq_idx_WB)/4)*2 4 ceil(presums(freq_idx_WB)/4)*2];
% Switch from tx calibration window to hanning window to broaden beam
DDS_amp = final_DDS_amp{freq_idx_NB} .* hanning(8).' ./ Hwindow_orig;
% Renormalize the amplitudes
[~,relative_max_idx] = max(DDS_amp./param.max_tx);
DDS_amp = round(DDS_amp .* param.max_tx(relative_max_idx) / DDS_amp(relative_max_idx));
param.tx_weights = DDS_amp;
param.tukey = 0.08;
param.wfs(1).Tpd = 3e-6;
param.wfs(2).Tpd = 3e-6;
param.wfs(3).Tpd = 3e-6;
param.wfs(1).phase = final_DDS_phase{freq_idx_NB};
param.wfs(2).phase = final_DDS_phase{freq_idx_NB};
param.wfs(3).phase = final_DDS_phase{freq_idx_NB};
param.wfs(1).name = 'left';
param.wfs(2).name = 'nadir';
param.wfs(3).name = 'right';
% Add in time delays to each position, subtract out the nadir time delays since tx_equalization already took care of those
beam_angle_deg = 20; % Positive to the left
param.wfs(1).delay = final_DDS_time{freq_idx_NB} ...
  - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
  - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
  + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
beam_angle_deg = 0; % Nadir
param.wfs(2).delay = final_DDS_time{freq_idx_NB} ...
  - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
  - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
  + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
beam_angle_deg = -20; % Negative to the right
param.wfs(3).delay = final_DDS_time{freq_idx_NB} ...
  - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
  - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
  + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
param.wfs(1).f0 = f0_list(freq_idx_NB);
param.wfs(2).f0 = f0_list(freq_idx_NB);
param.wfs(3).f0 = f0_list(freq_idx_NB);
param.wfs(1).f1 = f1_list(freq_idx_NB);
param.wfs(2).f1 = f1_list(freq_idx_NB);
param.wfs(3).f1 = f1_list(freq_idx_NB);
param.DDC_freq = (param.wfs(2).f0+param.wfs(2).f1)/2;
[param.wfs(1:3).tx_mask] = deal([0 0 0 0 0 0 0 0]);
param.wfs(1).atten = 0;
param.wfs(2).atten = 0;
param.wfs(3).atten = 0;
param.fn = fullfile(base_dir,sprintf('image_angelika_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml', ...
  param.wfs(2).f0/1e6,param.wfs(2).f1/1e6,param.tg.Haltitude*100/2.54/12,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));
write_cresis_xml(param);

%% Image Mode (Low Altitude, Ice <1400 m thick)
% Ice thickness "param.tg.Hice_thick_min" m to "param.tg.Hice_thick" m, "param.tg.Haltitude" +/- "param.tg.altitude_guard" ft AGL
freq_idx_WB = 1;
freq_idx_NB = 2;
param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 700; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
param.max_data_rate = 755;
param.DDC_select = DDC_select_list(freq_idx_NB);
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.staged_recording = [0 0 0];
param.tg.start_ref = {'bottom','surface','bottom'};
param.tg.stop_ref = {'bottom','bottom','bottom'};
param.tg.altitude_guard = 500 * 12*2.54/100;
param.tg.Haltitude = 1200 * 12*2.54/100;
param.tg.Hice_thick_min = 0;
param.tg.Hice_thick = 1400;
param.tg.look_angle_deg = [40 0 40];
param.prf = prf;
param.presums = [8 8 8]; %[ceil(presums(freq_idx_WB)/4)*2 4 ceil(presums(freq_idx_WB)/4)*2];
% Switch from tx calibration window to hanning window to broaden beam
DDS_amp = final_DDS_amp{freq_idx_NB} .* hanning(8).' ./ Hwindow_orig;
% Renormalize the amplitudes
[~,relative_max_idx] = max(DDS_amp./param.max_tx);
DDS_amp = round(DDS_amp .* param.max_tx(relative_max_idx) / DDS_amp(relative_max_idx));
param.tx_weights = DDS_amp;
param.tukey = 0.08;
param.wfs(1).Tpd = 3e-6;
param.wfs(2).Tpd = 3e-6;
param.wfs(3).Tpd = 3e-6;
param.wfs(1).phase = final_DDS_phase{freq_idx_NB};
param.wfs(2).phase = final_DDS_phase{freq_idx_NB};
param.wfs(3).phase = final_DDS_phase{freq_idx_NB};
param.wfs(1).name = 'left';
param.wfs(2).name = 'nadir';
param.wfs(3).name = 'right';
% Add in time delays to each position, subtract out the nadir time delays since tx_equalization already took care of those
beam_angle_deg = 20; % Positive to the left
param.wfs(1).delay = final_DDS_time{freq_idx_NB} ...
  - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
  - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
  + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
beam_angle_deg = 0; % Nadir
param.wfs(2).delay = final_DDS_time{freq_idx_NB} ...
  - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
  - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
  + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
beam_angle_deg = -20; % Negative to the right
param.wfs(3).delay = final_DDS_time{freq_idx_NB} ...
  - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
  - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
  + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
param.wfs(1).f0 = f0_list(freq_idx_NB);
param.wfs(2).f0 = f0_list(freq_idx_NB);
param.wfs(3).f0 = f0_list(freq_idx_NB);
param.wfs(1).f1 = f1_list(freq_idx_NB);
param.wfs(2).f1 = f1_list(freq_idx_NB);
param.wfs(3).f1 = f1_list(freq_idx_NB);
param.DDC_freq = (param.wfs(2).f0+param.wfs(2).f1)/2;
[param.wfs(1:3).tx_mask] = deal([0 0 0 0 0 0 0 0]);
param.wfs(1).atten = 0;
param.wfs(2).atten = 0;
param.wfs(3).atten = 0;
param.fn = fullfile(base_dir,sprintf('image_angelika_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml', ...
  param.wfs(2).f0/1e6,param.wfs(2).f1/1e6,param.tg.Haltitude*100/2.54/12,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));
write_cresis_xml(param);


%% Image Mode (Low Altitude, Thin Ice)

% Ice thickness 0 m to "param.tg.Hice_thick" m, "param.tg.Haltitude" +/- "param.tg.altitude_guard" ft AGL

for freq_idx = 2

  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);

  param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 700; param.flight_hours = 4; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;

  param.max_data_rate = 755;

  param.DDC_select = DDC_select_list(freq_idx);

  param.max_duty_cycle = 0.12;

  param.create_IQ = false;

  param.tg.staged_recording = [1 1 2 2 3 3];

  param.tg.altitude_guard = 700 * 12*2.54/100;

  param.tg.Haltitude = 1200 * 12*2.54/100;

  param.tg.Hice_thick = 2000;

  param.tg.look_angle_deg = [40 40 40 40 40 40];

  param.prf = prf;

  param.presums = [2 2 4 4 ceil((presums(freq_idx)-12)/4)*2 ceil((presums(freq_idx)-12)/4)*2];

  % Switch from tx calibration window to hanning window to broaden beam

  DDS_amp = final_DDS_amp{freq_idx} .* hanning(8).' ./ Hwindow_orig;

  % Renormalize the amplitudes

  [~,relative_max_idx] = max(DDS_amp./param.max_tx);

  DDS_amp = round(DDS_amp .* param.max_tx(relative_max_idx) / DDS_amp(relative_max_idx));

  param.tx_weights = DDS_amp;

  param.tukey = 0.08;

  param.wfs(1).Tpd = 1e-6;

  param.wfs(2).Tpd = 1e-6;

  param.wfs(3).Tpd = 3e-6;

  param.wfs(4).Tpd = 3e-6;

  param.wfs(5).Tpd = 10e-6;

  param.wfs(6).Tpd = 10e-6;

  param.wfs(1).phase = final_DDS_phase{freq_idx};

  param.wfs(2).phase = final_DDS_phase{freq_idx};

  param.wfs(3).phase = final_DDS_phase{freq_idx};

  param.wfs(4).phase = final_DDS_phase{freq_idx};

  param.wfs(5).phase = final_DDS_phase{freq_idx};

  param.wfs(6).phase = final_DDS_phase{freq_idx};

  param.wfs(1).name = 'left1';

  param.wfs(2).name = 'right1';

  param.wfs(3).name = 'left3';

  param.wfs(4).name = 'right3';

  param.wfs(5).name = 'left10';

  param.wfs(6).name = 'right10';

  % Add in time delays to each position, subtract out the nadir time delays since tx_equalization already took care of those

  beam_angle_deg = 20; % Positive to the left

  param.wfs(1).delay = final_DDS_time{freq_idx} ...

    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...

    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...

    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;

  beam_angle_deg = -20; % Positive to the right

  param.wfs(2).delay = final_DDS_time{freq_idx} ...

    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...

    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...

    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;

  beam_angle_deg = 20; % Positive to the left

  param.wfs(3).delay = final_DDS_time{freq_idx} ...

    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...

    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...

    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;

  beam_angle_deg = -20; % Negative to the right

  param.wfs(4).delay = final_DDS_time{freq_idx} ...

    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...

    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...

    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;

  beam_angle_deg = 20; % Positive to the left

  param.wfs(5).delay = final_DDS_time{freq_idx} ...

    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...

    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...

    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;

  beam_angle_deg = -20; % Negative to the right

  param.wfs(6).delay = final_DDS_time{freq_idx} ...

    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...

    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...

    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;

  param.f0 = f0_list(freq_idx);

  param.f1 = f1_list(freq_idx);

  param.DDC_freq = (param.f0+param.f1)/2;

  [param.wfs(1:6).tx_mask] = deal([0 0 0 0 0 0 0 0]);

  param.wfs(1).atten = 37;

  param.wfs(2).atten = 37;

  param.wfs(3).atten = 0;

  param.wfs(4).atten = 0;

  param.wfs(5).atten = 0;

  param.wfs(6).atten = 0;

  param.fn = fullfile(base_dir,sprintf('image_lowthin_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml', ...
      param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/2.54/12,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));

  write_cresis_xml(param);

end

