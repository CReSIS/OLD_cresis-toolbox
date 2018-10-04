function [param,defaults] = default_radar_params_2018_Greenland_P3_accum
% [param,defaults] = default_radar_params_2018_Greenland_P3_accum
%
% Accum: 2018 Greenland P3
%
% Creates base "param" struct
% Creates defaults cell array for each type of radar setting
%
% Author: John Paden

%% Preprocess parameters
param.season_name = '2018_Greenland_P3';
param.radar_name = 'mcords5-accum';

param.preprocess.daq.type = 'cresis';
param.preprocess.daq.xml_version = 2.0;
param.preprocess.daq.header_load_func = @basic_load_mcords5;
param.preprocess.daq.board_map = {'chan1','chan2','chan3','chan4'};
param.preprocess.daq.clk = 1e9/8;

param.preprocess.wg.type = 'arena';
param.preprocess.wg.tx_map = {'awg0','awg1'};

param.preprocess.file.version = 103;
param.preprocess.file.prefix = param.radar_name;
param.preprocess.file.suffix = '.bin';
param.preprocess.max_time_gap = 10;
param.preprocess.min_seg_size = 1;

%% OIB Accum Arena Parameters
arena = [];
arena.awg = [];
arena.awg(end+1).awg = 2;
arena.awg(end).dacs = [4 5];
arena.awg(end).dacClk = [2400e6 2400e6];
arena.awg(end).desiredAlignMin = [13 8];
arena.awg(end).desiredAlignMax = [27 22];
arena.awg(end+1).awg = 3;
arena.awg(end).dacs = [6 7];
arena.awg(end).dacClk = [2400e6 2400e6];
arena.awg(end).desiredAlignMin = [9 24];
arena.awg(end).desiredAlignMax = [23 38];
arena.dacs = [4 5 6 7];
arena.dacs_sampFreq = [2400e6 2400e6 2400e6 2400e6];
arena.max_tx = 0.7*[1 1 1 1 1 1 1 1]; %for less power add "*.5" 
arena.zeropimods = [0 180 270 90];
arena.TTL_time = [0.1, 0.2, (0.2+3+0.3)];
arena.dacs_internal_delay = 1.8;
arena.dacs_start_delay = 3.2;

arena.TTL_names = {};
arena.TTL_names{end+1} = 'test';
arena.TTL_names{end+1} = 'EPRI';
arena.TTL_names{end+1} = 'PRI';
arena.TTL_names{end+1} = 'PA_ENA';
arena.TTL_states{1} = [
  0 1 0 0 % test
  1 0 0 0 % EPRI
  0 1 0 0 % PRI
  0 0 1 0 % PA ENA
  ];
arena.TTL_states{2} = [
  0 1 0 0 % test
  0 0 0 0 % EPRI
  0 1 0 0 % PRI
  0 0 1 0 % PA ENA
  ];

param.preprocess.arena = arena;

% %% Control parameters (not used in the parameter spreadsheet directly)
% default.xml_file_prefix = 'mcords5';
% default.data_file_prefix = 'mcords5';
% default.header_load_func = @basic_load_mcords5;
% default.header_load_params = struct('clk',1600e6,'presum_bug_fixed',true);
% default.xml_version = 2.0;
% 
% default.noise_50ohm = [0 0 0 0];
% 
% default.Pt = 4*400;
% default.Gt = 4*4;
% default.Ae = 2*0.2*0.2;
% 
% default.system_loss_dB = 10.^(-5.88/10);
% default.max_DDS_RAM = 4000;
% default.tx_voltage = sqrt(1000*50)*10^(-2/20);
% 
% default.iq_mode = 0;
% default.tx_DDS_mask = [1 1 1 1];
% 
% default.radar_worksheet_headers = {'Tpd','Tadc','Tadc_adjust','f0','f1','ref_fn','tukey','tx_weights','rx_paths','adc_gains','chan_equal_dB','chan_equal_deg','Tsys','DC_adjust','DDC_mode','DDC_freq'};
% default.radar_worksheet_headers_type = {'r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r'};
% 
% default.basic_surf_track_min_time = -inf; % Normally -inf for lab test, 2e-6 for flight test
% default.basic_surf_track_Tpd_factor = -inf; % Normally -inf for lab test, 1.1 for flight test
% default.adc_folder_name = 'chan%d';

if 1
  % Example 1: Normal configuration:
  %   Connect antenna N to WFG N for all N = 1 to 8
  ref_adc = 3;
  default.txequal.img = [(1:4).', ref_adc*ones(4,1)];
  default.txequal.ref_wf_adc = 3;
  default.txequal.wf_mapping = [1 2 3 4 ];
  %default.txequal.Hwindow_desired = [chebwin(4,30).' 0 0 0 0];
  default.txequal.Hwindow_desired = [0.7 1 1 0.7 0 0 0 0];
  default.txequal.max_DDS_amp = [4000 4000 4000 4000 4000 4000 4000 4000 ];
  default.txequal.time_delay_desired = [0 0 0 0 ];
  default.txequal.phase_desired = [0 0 0 0 ];
  default.txequal.time_validation = [0.4 0.4 0.4 0.4 0 0 0 0]*1e-9;
  default.txequal.amp_validation = [3 3 3 3 0 0 0 0];
  default.txequal.phase_validation = [35 35 35 35 0 0 0 0];
  default.txequal.remove_linear_phase_en = true;
end

%% Records worksheet
default.records.gps.time_offset = 1;
default.records.geotiff_fn = 'greenland/Landsat-7/Greenland_natural_150m';
default.records.file.adcs = 1;
default.records.file.adc_headers = 1;
default.records.gps.en = 1;
default.records.frame_mode = 0;
default.records.presum_bug_fixed = 1;
default.records.tmp_fn_uses_adc_folder_name = 1;

%% Qlook worksheet
default.qlook.out_path = '';
default.qlook.en = 1;
default.qlook.block_size = 5000;
default.qlook.frm_types = {0,[0 1],0,0,-1};
default.qlook.coh_noise_method = [];
default.qlook.coh_noise_arg = [];
default.qlook.ft_wind = @hanning;
default.qlook.ft_wind_time = false;
default.qlook.ft_dec = true;
default.qlook.pulse_comp = [];
default.qlook.pulse_rfi.en = [];
default.qlook.pulse_rfi.inc_ave= [];
default.qlook.pulse_rfi.thresh_scale = [];
default.qlook.roll_correction = 0;
default.qlook.lever_arm_fh = @lever_arm;
default.qlook.elev_correction = 0;
default.qlook.B_filter = ones(1,20)/20;
default.qlook.decimate_factor = 20;
default.qlook.inc_ave = 10;
default.qlook.surf.en = 1;
default.qlook.surf.method = 'threshold';
default.qlook.surf.noise_rng = [0 -50 10];
default.qlook.surf.min_bin = 2e-6;
default.qlook.surf.max_bin = [];
default.qlook.surf.threshold = 9;
default.qlook.surf.sidelobe = 15;
default.qlook.surf.medfilt = 3;
default.qlook.surf.search_rng = [0:2];

%% SAR worksheet
default.sar.out_path = '';
default.sar.imgs = {[1 1],[2 1]};
default.sar.frm_types = {0,[0 1],0,0,-1};
default.sar.chunk_len = 5000;
default.sar.chunk_overlap = 10;
default.sar.frm_overlap = 0;
default.sar.coh_noise_removal = 0;
default.sar.combine_rx = 0;
default.sar.time_of_full_support = inf;
default.sar.pulse_rfi.en = [];
default.sar.pulse_rfi.inc_ave= [];
default.sar.pulse_rfi.thresh_scale = [];
default.sar.trim_vals = [];
default.sar.pulse_comp = 1;
default.sar.ft_dec = 1;
default.sar.ft_wind = @hanning;
default.sar.ft_wind_time = 0;
default.sar.lever_arm_fh = @lever_arm;
default.sar.mocomp.en = 1;
default.sar.mocomp.type = 2;
default.sar.mocomp.filter = {@butter  [2]  [0.1000]};
default.sar.mocomp.uniform_en = 1;
default.sar.sar_type = 'f-k';
default.sar.sigma_x = 2.5;
default.sar.sub_aperture_steering = 0;
default.sar.st_wind = @hanning;
default.sar.start_eps = 3.15;

%% Array worksheet
default.array.in_path = '';
default.array.array_path = '';
default.array.out_path = '';
default.array.method = 'standard';
default.array.window = @hanning;
default.array.bin_rng = 0;
default.array.rline_rng = -5:5;
default.array.dbin = 1;
default.array.dline = 6;
default.array.DCM = [];
default.array.three_dim.en = 0;
default.array.three_dim.layer_fn = '';
default.array.Nsv = 1;
default.array.theta_rng = [0 0];
default.array.sv_fh = @array_proc_sv;
default.array.diag_load = 0;
default.array.Nsig = 2;

%% Radar worksheet
default.radar.fs = 1600e6;
default.radar.Tadc = []; % normally leave empty to use value in file header
default.radar.adc_bits = 12;
default.radar.adc_full_scale = 2;
default.radar.rx_paths = [1:4];
default.radar.noise_figure = 2;
default.radar.rx_gain = 45;
default.radar.adc_SNR_dB = 59;
default.radar.Tadc_adjust = 6.95E-06; % System time delay: leave this empty or set it to zero at first, determine this value later using data over surface with known height or from surface multiple

%% Post worksheet
default.post.data_dirs = {'qlook'};
default.post.layer_dir = 'layerData';
default.post.maps_en = 1;
default.post.echo_en = 1;
default.post.layers_en = 0;
default.post.data_en = 0;
default.post.csv_en = 1;
default.post.concat_en = 1;
default.post.pdf_en = 1;
default.post.map.location = 'Greenland';
default.post.map.type = 'combined';
default.post.echo.elev_comp = 2;
default.post.echo.depth = '[min(Surface_Depth)-20 max(Surface_Depth)+200]';
% default.post.echo.elev_comp = 3;
% default.post.echo.depth = '[min(Surface_Elev)-200 max(Surface_Elev)+20]';
default.post.echo.er_ice = 3.15;
default.post.ops.location = 'arctic';
  
%% Radar Settings

defaults = {};

% Wideband settings
default.radar.wfs(1).chan_equal_Tsys = [0.15	0.00	0.10	0.68]/1e9;
default.radar.wfs(1).chan_equal_dB = [0 0 0 0];
default.radar.wfs(1).chan_equal_deg = [-105.5	-0.0	55.2	128.8];

 % survey mode
default.qlook.img_comb = [1e-06 -inf 2e-06];
default.qlook.imgs = {[1 1],[2 1]};
default.array.imgs = default.qlook.imgs;
default.array.img_comb = default.qlook.img_comb;
default.radar.DC_adjust = {'','','','',''};
default.radar.ref_fn = '';
default.xml_regexp = 'survey_600-900MHz_.*.xml';
default.name = 'Survey Mode 600-900 MHz';
defaults{end+1} = default;

% Other settings
default.qlook.img_comb = [];
default.qlook.imgs = {[1 1],[2 1]};
default.sar.imgs = {[1 1],[2 1]};
default.array.imgs = default.qlook.imgs;
default.array.img_comb = default.qlook.img_comb;
default.radar.DC_adjust = {'','','','',''};
default.radar.ref_fn = '';

default.xml_regexp = '.*';
default.name = 'Default 600-900 MHz';
defaults{end+1} = default;

%% Add default settings

param.preprocess.defaults = defaults;
