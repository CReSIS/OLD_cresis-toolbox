% script create_settings_2017_Greenland_P3_Accum
%
% Creates NI radar depth sounder settings
%
% See link_budgets.xls

physical_constants; % c = speed of light

% Define waveforms
if ispc
  base_dir = 'c:\waveforms_accum\';
  rss_base_dir = 'c:\waveforms_accum_rss\';
else
  base_dir = '~/waveforms/';
  rss_base_dir = '~/rss_waveforms/';
end

f0_list = [600e6];
f1_list = [900e6];
DDC_select_list = [0]; % Which DDC mode to use
cal_settings = [1];
prf = 10000;

% presums: Ensure lambda/4 sampling (fudge factor allows difference) and an
%   that presums are an even number.
velocity = [125]; % m/s
presums = round(c./abs(f1_list+f0_list)/2 ./ velocity * prf / 4)*4;
presums = round(c./abs(f1_list+f0_list)/2 ./ velocity * prf / 4 * 128/40)*4

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
    final_DDS_time{idx} =  [0 0 0 0 0 0 0 0];
  end
else
  % COPY AND PASTE RESULTS FROM basic_tx_chan_equalization_SEASON_NAME.m
  % HERE:
  
  % These are from transmit calibration during 20160401 test flight  
  % NOTE: These values are valid for when DDS channels 5-8 come up one 1440
  % MHz clock cycle after channels 1-4.
  
  % 600-900 MHz
  final_DDS_phase{1} = [0 0 0 0 0 0 0 0];
  final_DDS_phase_no_time = [0 0 0 0 0 0 0 0]; % not used usually
  final_DDS_amp{1} = [4000 4000 4000 4000 4000 4000 4000 4000];
  final_DDS_time{1} =  [0 0 0 0 0 0 0 0];

end

% Hwindow_orig: Desired window created during transmit calibration
%  This is used any time a window that is different from that used
%  during calibration is to be used.
Hwindow_orig = boxcar(8).';

%% SETUP
% =========================================================================

param = [];
param.season_name = '2017_Greenland_P3';
param.radar_name = 'accum';
param.gps_source = 'atm-?';
clear phase_centers;
for tx_chan = 1
  % Just enable the current antenna to determine its phase center
  tx_weights = 1;
  rxchan = 2; % Fix the receiver (it should not matter which one you choose)
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

%% Survey Mode + loopback, noise, and deconv modes
% <4000 m thick ice, 1200 +/- 700 ft AGL
ice_thickness = [500];
for freq_idx = [1]
  param = struct('radar_name','mcords5','num_chan',4,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',150.0e6,'fs_dds',2400e6,'TTL_clock',2400e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'rss_base_dir',rss_base_dir);
  param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 96; param.flight_hours = 7; param.sys_delay = 0.75e-6; param.use_mcords4_names = true;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [1 2];
  param.tg.altitude_guard = 0*12*2.54/100;
  param.tg.Haltitude = 1500*12*2.54/100;
  param.tg.Hice_thick = ice_thickness(freq_idx);
  param.prf = prf;
  param.presums = [4 presums(freq_idx)-4];
  param.wfs(1).atten = 28;
  param.wfs(2).atten = 0;
  DDS_amp = final_DDS_amp{cal_settings(freq_idx)};
  param.tx_weights = DDS_amp;
  param.tukey = 0.05;
  param.wfs(1).Tpd = 2e-6;
  param.wfs(2).Tpd = 3e-6;
  param.wfs(1).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(2).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = (param.f0+param.f1)/2;
  [param.wfs(1:2).tx_mask] = deal([0 0 0 0 0 0 0 0]);
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
  param.wfs(1).atten = 31;
  param.wfs(2).atten = 31;
  param.tg.staged_recording = false;
  param.tg.altitude_guard = 0*12*2.54/100;
  param.tg.Haltitude = 20000*12*2.54/100;
  param.tg.Hice_thick = 0 * 12*2.54/100/sqrt(er_ice);
  param.fn = fullfile(calval_dir,sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_DECONV.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6));
  write_cresis_xml(param);
  if freq_idx == 1
    % Noise Mode
    param.tx_weights = [0 0 0 0 0 0 0 0];
    [param.wfs(1:2).tx_mask] = deal([1 1 1 1 1 1 1 1]);
    param.wfs(1).atten = 28;
    param.wfs(2).atten = 0;
    param.tg.staged_recording = [1 2];
    param.tg.altitude_guard = 500*12*2.54/100;
    param.tg.Haltitude = 1400*12*2.54/100;
    param.tg.Hice_thick = 500;
    param.fn = fullfile(calval_dir,sprintf('survey_%.0f-%.0fMHz_%.0fus_NOISE.xml',param.f0/1e6,param.f1/1e6,param.wfs(end).Tpd*1e6));
    write_cresis_xml(param);
  end
end


