function param = default_radar_params_2018_Antarctica_DC8_rds
% param = default_radar_params_2018_Antarctica_DC8_rds
%
% rds: 2018_Antarctica_DC8
%
% Creates base "param" struct
% Creates defaults cell array for each type of radar setting
%
% Author: John Paden

param.season_name = '2018_Antarctica_DC8';
param.radar_name = 'mcords3';

%% Control parameters
% default.xml_file_prefix = 'mcords3';
% default.data_file_prefix = 'mcords3';
% default.header_load_func = @basic_load_mcords3;
% default.header_load_params = struct('clk',1e9/9,'presum_bug_fixed',false);
% default.xml_version = 2.0;

% default.noise_50ohm = [-45.5	-44.3	-44.8	-44.9	-40.2	-40.8	-40.5	-41.4	-39.9	-40.2	-42.0	-43.6	-43.0	-44.1	-44.6];
% 
% default.Pt = 1000 * [1 1 1 1 1 1 1];
% default.Gt = 7*4;
% default.Ae = 2*0.468 * 0.468;
% 
% default.system_loss_dB = 10.^(-5.88/10);
% default.max_DDS_RAM = 40000;
% default.tx_voltage = sqrt(1000*50)*10^(-2/20);
% 
% default.iq_mode = 0;
% default.tx_DDS_mask = [1 1 1 1 1 1 1 0];
% 
% default.radar_worksheet_headers = {'Tpd','Tadc','Tadc_adjust','f0','f1','ref_fn','tukey','tx_weights','rx_paths','adc_gains','chan_equal_dB','chan_equal_deg','Tsys','DC_adjust','DDC_mode','DDC_freq'};
% default.radar_worksheet_headers_type = {'r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r'};
% 
% default.basic_surf_track_min_time = 2e-6; % Normally -inf for lab test, 2e-6 for flight test
% default.basic_surf_track_Tpd_factor = 1.1; % Normally -inf for lab test, 1.1 for flight test
% default.adc_folder_name = 'board%b';

if 1
  % Example 1: Normal configuration:
  %   Connect antenna N to WFG N for all N = 1 to 7
  ref_adc = 14;
  default.txequal.img = [(1:7).', ref_adc*ones(7,1)];
  default.txequal.ref_wf_adc = 3;
  default.txequal.wf_mapping = [1 2 3 4 5 6 7 0];
  default.txequal.Hwindow_desired = [1 1 1 1 1 1 1 0];
  default.txequal.max_DDS_amp = [40000 40000 40000 40000 40000 40000 40000 0];
  default.txequal.time_delay_desired = [0 0 0 0 0 0 0 0];
  default.txequal.phase_desired = [0 0 0 0 0 0 0 0];
  default.txequal.time_validation = [3 3 3 3 3 3 3 3]*1e-9;
  default.txequal.amp_validation = [3 3 3 3 3 3 3 3];
  default.txequal.phase_validation = [35 35 35 35 35 35 35 35];
  default.txequal.remove_linear_phase_en = true;
end

%% Vectors worksheet in parameter spreadsheet
default.vectors.gps.time_offset = 1;

%% Records worksheet in parameter spreadsheet
default.records.geotiff_fn = 'greenland\Landsat-7\mzl7geo_90m_lzw.tif';
default.records.file.adcs = [2:16];
default.records.file.adc_headers = [2:16];
default.records.gps.en = 1;
default.records.frame_mode = 0;
default.records.presum_bug_fixed = 0;
default.records.tmp_fn_uses_adc_folder_name = 1;

%% Get heights (quick-look) worksheet in parameter spreadsheet
default.qlook.out_path = '';
default.qlook.en = 1;
default.qlook.block_size = 10000;
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

%% CSARP worksheet in parameter spreadsheet
default.sar.out_path = '';
default.sar.imgs = {[1*ones(4,1),(2:5).'],[2*ones(4,1),(2:5).'],[3*ones(4,1),(2:5).']};
default.sar.frm_types = {0,[0 1],0,0,-1};
default.sar.chunk_len = 5000;
default.sar.chunk_overlap = 10;
default.sar.frm_overlap = 0;
default.sar.coh_noise_removal = 0;
default.sar.array_rx = 0;
default.sar.time_of_full_support = 3.5e-5;
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

%% Combine worksheet in parameter spreadsheet
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

%% Radar worksheet in parameter spreadsheet
default.radar.fs = 1e9/9;
default.radar.Tadc = []; % normally leave empty to use value in file header
default.radar.adc_bits = 14;
default.radar.adc_full_scale = 2;
default.radar.rx_paths = [1,1:15];
default.radar.rx_paths = [1 8 9 10 11 1 2 3 4 5 6 7 12 13 14 15];
default.radar.noise_figure = 2;
default.radar.rx_gain = 51.5;
default.radar.adc_SNR_dB = 70;
default.radar.Tadc_adjust = -1.4455e-06; % System time delay: leave this empty or set it to zero at first, determine this value later using data over surface with known height or from surface multiple

defaults = {};

%% Settings
default.radar.wfs(1).chan_equal_Tsys = [16.51	11.76	8.57	5.91	1.60	-0.13	0.00	-2.44	-4.83	-3.16	-2.76	-10.60	-6.05	-1.68	5.19]/1e9;
default.radar.wfs(1).chan_equal_dB = [1.1	1.0	1.0	3.0	3.9	1.0	0.0	-2.6	-0.7	0.9	1.3	3.7	2.0	0.2	3.3];
default.radar.wfs(1).chan_equal_deg = [97.1	25.2	-60.1	-126.1	-38.7	6.4	-0.0	-131.0	0.6	-25.1	-47.5	-86.3	-33.5	165.4	-106.6];

% Survey mode
default.get_heights.qlook.img_comb = [3e-06 -inf 1e-06 1e-05 -inf 3e-06];
default.get_heights.imgs = {[1*ones(6,1),(1:6).'],[2*ones(6,1),(1:6).'],[3*ones(6,1),(1:6).']};
default.sar.imgs = default.qlook.imgs;
default.array.imgs = default.qlook.imgs;
default.array.img_comb = default.qlook.img_comb;
default.radar.ref_fn = '';
default.xml_regexp = 'survey_.*thick.xml';
default.name = 'Survey Mode';
defaults{end+1} = default;

% High altitude mode
default.get_heights.qlook.img_comb = [3e-06 -inf 1e-06 1e-05 -inf 3e-06];
default.get_heights.imgs = {[1*ones(6,1),(1:6).'],[2*ones(6,1),(1:6).'],[3*ones(6,1),(1:6).']};
default.sar.imgs = default.qlook.imgs;
default.array.imgs = default.qlook.imgs;
default.array.img_comb = default.qlook.img_comb;
default.radar.ref_fn = '';
default.xml_regexp = 'survey_.*high_altitude.xml';
default.name = 'High Altitude Mode';
defaults{end+1} = default;

% Deconvolution mode
default.get_heights.qlook.img_comb = [];
default.get_heights.imgs = {[1*ones(6,1),(1:6).'],[2*ones(6,1),(1:6).'],[3*ones(6,1),(1:6).']};
default.sar.imgs = default.qlook.imgs;
default.array.imgs = default.qlook.imgs;
default.array.img_comb = default.qlook.img_comb;
default.radar.ref_fn = '';
default.xml_regexp = 'survey_.*DECONVOLUTION.xml';
default.name = 'Deconvolution Mode';
defaults{end+1} = default;

% Other settings
default.qlook.img_comb = [];
default.qlook.imgs = [];
default.array.imgs = default.qlook.imgs;
default.array.img_comb = default.qlook.img_comb;
default.radar.ref_fn = '';
default.xml_regexp = '.*';
default.name = 'Other Settings';
defaults{end+1} = default;
