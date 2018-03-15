function [param,defaults] = default_radar_params_2018_Greenland_P3_mcords
% [param,defaults] = default_radar_params_2018_Greenland_P3_mcords
%
% MCORDS 3: 2018 Greenland P3
%
% Creates base "param" struct
% Creates defaults cell array for each type of radar setting
%
% Author: John Paden

param.season_name = '2018_Greenland_P3';
param.radar_name = 'mcords3';

%% Control parameters (not used in the parameter spreadsheet directly)
default.xml_file_prefix = 'mcords3';
default.data_file_prefix = 'mcords3';
default.header_load_func = @basic_load_mcords3;
default.header_load_params = struct('clk',1e9/9,'presum_bug_fixed',false);
default.xml_version = 2.0;

default.noise_50ohm = [-45.5	-44.3	-44.8	-44.9	-40.2	-40.8	-40.5	-41.4	-39.9	-40.2	-42.0	-43.6	-43.0	-44.1	-44.6];

default.Pt = 500 * [1 1 1 1 1 1 1];
default.Gt = 7*4;
default.Ae = 2*0.468 * 0.468;

default.system_loss_dB = 10.^(-5.88/10);
default.max_DDS_RAM = 40000;
default.tx_voltage = sqrt(1000*50)*10^(-2/20);

default.iq_mode = 0;
default.tx_DDS_mask = [1 1 1 1 1 1 1 0];

default.radar_worksheet_headers = {'Tpd','Tadc','Tadc_adjust','f0','f1','ref_fn','tukey','tx_weights','rx_paths','adc_gains','chan_equal_dB','chan_equal_deg','Tsys','DC_adjust','DDC_mode','DDC_freq'};
default.radar_worksheet_headers_type = {'r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r'};

default.basic_surf_track_min_time = 2e-6; % Normally -inf for lab test, 2e-6 for flight test
default.basic_surf_track_Tpd_factor = 1.1; % Normally -inf for lab test, 1.1 for flight test
default.adc_folder_name = 'board%b';

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
default.get_heights.qlook.out_path = '';
default.get_heights.qlook.en = 1;
default.get_heights.block_size = 10000;
default.get_heights.frm_types = {0,[0 1],0,0,-1};
default.get_heights.coh_noise_method = [];
default.get_heights.coh_noise_arg = [];
default.get_heights.ft_wind = @hanning;
default.get_heights.ft_wind_time = false;
default.get_heights.ft_dec = true;
default.get_heights.pulse_comp = [];
default.get_heights.pulse_rfi.en = [];
default.get_heights.pulse_rfi.inc_ave= [];
default.get_heights.pulse_rfi.thresh_scale = [];
default.get_heights.roll_correction = 0;
default.get_heights.lever_arm_fh = @lever_arm;
default.get_heights.elev_correction = 0;
default.get_heights.B_filter = ones(1,20)/20;
default.get_heights.decimate_factor = 20;
default.get_heights.inc_ave = 10;
default.get_heights.surf.en = 1;
default.get_heights.surf.method = 'threshold';
default.get_heights.surf.noise_rng = [0 -50 10];
default.get_heights.surf.min_bin = 2e-6;
default.get_heights.surf.max_bin = [];
default.get_heights.surf.threshold = 9;
default.get_heights.surf.sidelobe = 15;
default.get_heights.surf.medfilt = 3;
default.get_heights.surf.search_rng = [0:2];

%% CSARP worksheet in parameter spreadsheet
default.csarp.out_path = '';
default.csarp.imgs = {[1*ones(4,1),(2:5).'],[2*ones(4,1),(2:5).'],[3*ones(4,1),(2:5).']};
default.csarp.frm_types = {0,[0 1],0,0,-1};
default.csarp.chunk_len = 5000;
default.csarp.chunk_overlap = 10;
default.csarp.frm_overlap = 0;
default.csarp.coh_noise_removal = 0;
default.csarp.combine_rx = 0;
default.csarp.time_of_full_support = 3.5e-5;
default.csarp.pulse_rfi.en = [];
default.csarp.pulse_rfi.inc_ave= [];
default.csarp.pulse_rfi.thresh_scale = [];
default.csarp.trim_vals = [];
default.csarp.pulse_comp = 1;
default.csarp.ft_dec = 1;
default.csarp.ft_wind = @hanning;
default.csarp.ft_wind_time = 0;
default.csarp.lever_arm_fh = @lever_arm;
default.csarp.mocomp.en = 1;
default.csarp.mocomp.type = 2;
default.csarp.mocomp.filter = {@butter  [2]  [0.1000]};
default.csarp.mocomp.uniform_en = 1;
default.csarp.sar_type = 'f-k';
default.csarp.sigma_x = 2.5;
default.csarp.sub_aperture_steering = 0;
default.csarp.st_wind = @hanning;
default.csarp.start_eps = 3.15;

%% Combine worksheet in parameter spreadsheet
default.combine.in_path = '';
default.combine.array_path = '';
default.combine.out_path = '';
default.combine.method = 'standard';
default.combine.window = @hanning;
default.combine.bin_rng = 0;
default.combine.rline_rng = -5:5;
default.combine.dbin = 1;
default.combine.dline = 6;
default.combine.DCM = [];
default.combine.three_dim.en = 0;
default.combine.three_dim.layer_fn = '';
default.combine.Nsv = 1;
default.combine.theta_rng = [0 0];
default.combine.sv_fh = @array_proc_sv;
default.combine.diag_load = 0;
default.combine.Nsig = 2;

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

% survey mode
default.get_heights.qlook.img_comb = [3e-06 -inf 1e-06 1e-05 -inf 3e-06];
default.get_heights.imgs = {[1*ones(4,1),(2:5).'],[2*ones(4,1),(2:5).'],[3*ones(4,1),(2:5).']};
default.csarp.imgs = default.get_heights.imgs;
default.combine.imgs = default.get_heights.imgs;
default.combine.img_comb = default.get_heights.qlook.img_comb;
default.radar.DC_adjust = {'','',''};
default.radar.ref_fn = '';
default.xml_regexp = 'survey_.*thick.xml';
default.name = 'Survey Mode';
defaults{end+1} = default;

% survey mode
default.get_heights.qlook.img_comb = [3e-06 -inf 1e-06];
default.get_heights.imgs = {[1*ones(4,1),(2:5).'],[2*ones(4,1),(2:5).']};
default.csarp.imgs = default.get_heights.imgs;
default.combine.imgs = default.get_heights.imgs;
default.combine.img_comb = default.get_heights.qlook.img_comb;
default.radar.DC_adjust = {'','',''};
default.radar.ref_fn = '';
default.xml_regexp = 'survey_.*thin_ice.xml';
default.name = 'Thin Ice Mode';
defaults{end+1} = default;

% high altitude mode
default.get_heights.qlook.img_comb = [1e-05 -inf 3e-06];
default.get_heights.imgs = {[1*ones(4,1),(2:5).'],[2*ones(4,1),(2:5).']};
default.csarp.imgs = default.get_heights.imgs;
default.combine.imgs = default.get_heights.imgs;
default.combine.img_comb = default.get_heights.qlook.img_comb;
default.radar.DC_adjust = {'','',''};
default.radar.ref_fn = '';
default.xml_regexp = 'survey_.*high_altitude.xml';
default.name = 'High Altitude Mode';
defaults{end+1} = default;

% deconvolution mode
default.get_heights.qlook.img_comb = [];
default.get_heights.imgs = {[1*ones(4,1),(2:5).'],[2*ones(4,1),(2:5).'],[3*ones(4,1),(2:5).']};
default.csarp.imgs = default.get_heights.imgs;
default.combine.imgs = default.get_heights.imgs;
default.combine.img_comb = default.get_heights.qlook.img_comb;
default.radar.DC_adjust = {'','',''};
default.radar.ref_fn = '';
default.xml_regexp = 'survey_.*DECONVOLUTION.xml';
default.name = 'Deconvolution Mode';
defaults{end+1} = default;

%% Other settings

default.get_heights.qlook.img_comb = [];
default.get_heights.imgs = [];
default.combine.imgs = default.get_heights.imgs;
default.combine.img_comb = default.get_heights.qlook.img_comb;

default.radar.DC_adjust = {'DC_20160413_04_wf1.mat','DC_20160413_04_wf2.mat','DC_20160413_04_wf2.mat'};
default.radar.ref_fn = 'deconv_wf_%w_adc_%a_20160413_06';
default.xml_regexp = 'survey_180-210MHz_.*DECONV.xml';
default.name = 'Deconv 180-210 MHz';
defaults{end+1} = default;

default.radar.DC_adjust = [];
default.radar.ref_fn = '';
default.xml_regexp = '.*';
default.name = 'Other Settings';
defaults{end+1} = default;

return;
