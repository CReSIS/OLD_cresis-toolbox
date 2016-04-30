function [param,defaults] = default_radar_params_2016_Greenland_Polar6_mcords
% [param,defaults] = default_radar_params_2016_Greenland_Polar6_mcords
%
% MCORDS 5: 2016 Greenland Polar6
%
% Creates base "param" struct
% Creates defaults cell array for each type of radar setting
%
% Author: John Paden

param.season_name = '2016_Greenland_Polar6';
param.radar_name = 'mcords5';

%% Control parameters (not used in the parameter spreadsheet directly)
default.xml_file_prefix = 'mcords5';
default.data_file_prefix = 'mcords5';
default.header_load_func = @basic_load_mcords5;
default.header_load_params = struct('clk',1600e6,'presum_bug_fixed',true);
default.xml_version = 2.0;

% default.noise_50ohm = [-41.6	-42.2	-42.4	-41.9	-42.5	-42.9	-41.7	-43.0	-44.1	-44.7	-43.1	-44.1	-41.8	-42.6	-41.4	-42.6	-41.8	-43.1	-42.0	-42.7	-41.1	-43.4	-42.1	-41.9];
default.noise_50ohm = [-45.4	-45.7	-45.5	-45.6	-46.2	-46.7	-44.8	-46.1	-47.8	-48.8	-45.6	-47.5	-45.8	-46.1	-44.6	-45.1	-46.1	-46.7	-44.3	-45.8	-44.5	-46.5	-45.5	-44.9	];

default.Pt = (4*1000 + 4*500) * sum(chebwin(8,30).^2)/8;
default.Gt = 8*4;
default.Ae = 2*0.468 * 0.468;

default.system_loss_dB = 10.^(-5.88/10);
default.max_DDS_RAM = 4000;
default.tx_voltage = sqrt(1000*50)*10^(-2/20);

default.iq_mode = 0;
default.tx_DDS_mask = [1 1 1 1 1 1 1 1];

default.radar_worksheet_headers = {'Tpd','Tadc','Tadc_adjust','f0','f1','ref_fn','tukey','tx_weights','rx_paths','adc_gains','chan_equal_dB','chan_equal_deg','Tsys','DC_adjust','DDC_mode','DDC_freq'};
default.radar_worksheet_headers_type = {'r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r'};

%% Vectors worksheet in parameter spreadsheet
default.vectors.gps.time_offset = 1;

%% Records worksheet in parameter spreadsheet
default.records.geotiff_fn = 'greenland/Landsat-7/Greenland_natural_150m';
default.records.file.adcs = [1:24];
default.records.file.adc_headers = [1:24];
default.records.gps.en = 1;
default.records.frame_mode = 0;
default.records.presum_bug_fixed = 1;
default.records.tmp_fn_uses_adc_folder_name = 1;

%% Get heights (quick-look) worksheet in parameter spreadsheet
default.get_heights.qlook.out_path = '';
default.get_heights.qlook.en = 1;
default.get_heights.block_size = 5000;
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
default.csarp.imgs = {[1*ones(24,1),(1:24).'],[2*ones(24,1),(1:24).'],[3*ones(24,1),(1:24).']};
default.csarp.frm_types = {0,[0 1],0,0,-1};
default.csarp.chunk_len = 3500;
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
default.radar.fs = 1600e6;
default.radar.Tadc = []; % normally leave empty to use value in file header
default.radar.fs = 1600e6;
default.radar.adc_bits = 12;
default.radar.adc_full_scale = 2;
default.radar.rx_paths = [1:22,24,23];
default.radar.noise_figure = 2;
default.radar.rx_gain = 10^(48/20);
default.radar.adc_SNR_dB = 59;
default.radar.Tadc_adjust = 0.000010179163; % System time delay: leave this empty or set it to zero at first, determine this value later using data over surface with known height or from surface multiple

defaults = {};

%% Wideband settings
default.radar.wfs(1).chan_equal_dB = [-0.7 -0.7 -0.7 -0.2 -0.1 -2.3 0.2 -0.1 -4.2 -4.5 -3.1 -4.6 -1.1 -1.5 -0.9 -1.7 0 -0.6 0.2 3.4 0.7 -2.1 0.2 0.9];
default.radar.wfs(1).chan_equal_deg = [-168.6 -114.1 -5.7 9 30 24.1 -144.3 -137.7 113.1 64.8 124.3 133.7 108 138.1 71.6 102.9 -95.8 -127.2 -143.1 -139.6 -128.4 -172.3 -158.5 -152.4];
default.radar.wfs(1).chan_equal_Tsys = [0.3 0.7 0 0.2 0.1 0.2 0.2 0.3 -31.8 -32.1 -31.7 -31.6 -31.2 -30.8 -31.5 -31.1 -4.3 -4.5 -4.6 -4.5 -4.5 -4.6 -4.5 -4.5]/1e9;

 % survey mode
default.get_heights.qlook.img_comb = [3e-06 -inf 1e-06 1e-05 -inf 3e-06];
default.get_heights.imgs = {[1*ones(8,1),(9:16).'],[2*ones(8,1),(9:16).'],[3*ones(8,1),(9:16).']};
default.combine.imgs = default.get_heights.imgs;
default.combine.img_comb = default.get_heights.qlook.img_comb;
default.radar.DC_adjust = {'DC_20160413_04_wf1.mat','DC_20160413_04_wf2.mat','DC_20160413_04_wf3.mat'};
default.radar.ref_fn = '';
default.xml_regexp = 'survey_150-520MHz_.*thick.xml';
default.name = 'Survey Mode 150-520 MHz';
defaults{end+1} = default;

 % thin ice mode
default.get_heights.qlook.img_comb = [1e-06 -inf 1e-06 3e-06 -inf 1e-06];
default.get_heights.imgs = {[1*ones(8,1),(9:16).'],[2*ones(8,1),(9:16).'],[3*ones(8,1),(9:16).']};
default.combine.imgs = default.get_heights.imgs;
default.combine.img_comb = default.get_heights.qlook.img_comb;
default.radar.DC_adjust = {'DC_20160413_04_wf1.mat','DC_20160413_04_wf1.mat','DC_20160413_04_wf2.mat'};
default.radar.ref_fn = '';
default.xml_regexp = 'thinice_150-520MHz_.*thick.xml';
default.name = 'Thin Ice Mode 150-520 MHz';
defaults{end+1} = default;

 % 2 beam imaging mode
default.get_heights.qlook.img_comb = [];
default.get_heights.imgs = {[1*ones(8,1),(9:16).']};
default.combine.imgs = default.get_heights.imgs;
default.combine.img_comb = default.get_heights.qlook.img_comb;
default.radar.DC_adjust = {'DC_20160413_04_wf3.mat','DC_20160413_04_wf3.mat'};
default.radar.ref_fn = '';
default.xml_regexp = 'image_150-520MHz_.*thick.xml';
default.name = '2 Beam Image Mode 150-520 MHz';
defaults{end+1} = default;

 % 3 beam imaging mode
default.get_heights.qlook.img_comb = [];
default.get_heights.imgs = {[2*ones(8,1),(9:16).']};
default.combine.imgs = default.get_heights.imgs;
default.combine.img_comb = default.get_heights.qlook.img_comb;
default.radar.DC_adjust = {'DC_20160413_04_wf1.mat','DC_20160413_04_wf2.mat','DC_20160413_04_wf2.mat'};
default.radar.ref_fn = '';
default.xml_regexp = 'image3_150-520MHz_.*thick.xml';
default.name = '3 Beam Image Mode 150-520 MHz';
defaults{end+1} = default;

%% Narrowband settings
default.radar.wfs(1).chan_equal_dB = zeros(1,24);
default.radar.wfs(1).chan_equal_deg = zeros(1,24);
default.radar.wfs(1).chan_equal_Tsys = zeros(1,24);

% survey mode
default.get_heights.qlook.img_comb = [3e-06 -inf 1e-06 1e-05 -inf 3e-06];
default.get_heights.imgs = {[1*ones(8,1),(9:16).'],[2*ones(8,1),(9:16).'],[3*ones(8,1),(9:16).']};
default.combine.imgs = default.get_heights.imgs;
default.combine.img_comb = default.get_heights.qlook.img_comb;
default.radar.DC_adjust = {'DC_20160413_04_wf1.mat','DC_20160413_04_wf2.mat','DC_20160413_04_wf3.mat'};
default.radar.ref_fn = 'deconv_wf_%w_adc_%a_20160413_06';
default.xml_regexp = 'survey_180-210MHz_.*thick.xml';
default.name = 'Survey Mode 180-210 MHz';
defaults{end+1} = default;

% thin ice mode
default.get_heights.qlook.img_comb = [1e-06 -inf 1e-06 3e-06 -inf 1e-06];
default.get_heights.imgs = {[1*ones(8,1),(9:16).'],[2*ones(8,1),(9:16).'],[3*ones(8,1),(9:16).']};
default.combine.imgs = default.get_heights.imgs;
default.combine.img_comb = default.get_heights.qlook.img_comb;
default.radar.DC_adjust = {'DC_20160413_04_wf1.mat','DC_20160413_04_wf1.mat','DC_20160413_04_wf2.mat'};
default.radar.ref_fn = 'deconv_wf_%w_adc_%a_20160413_06_thin';
default.xml_regexp = 'thinice_180-210MHz_.*thick.xml';
default.name = 'Thin Ice Mode 180-210 MHz';
defaults{end+1} = default;

 % 2 beam imaging mode
default.get_heights.qlook.img_comb = [];
default.get_heights.imgs = {[1*ones(8,1),(9:16).']};
default.combine.imgs = default.get_heights.imgs;
default.combine.img_comb = default.get_heights.qlook.img_comb;
default.radar.DC_adjust = {'DC_20160413_04_wf3.mat','DC_20160413_04_wf3.mat'};
default.radar.ref_fn = 'deconv_wf_%w_adc_%a_20160413_06_image';
default.xml_regexp = 'image_180-210MHz_.*thick.xml';
default.name = '2 Beam Image Mode 180-210 MHz';
defaults{end+1} = default;

% 3 beam imaging mode
default.get_heights.qlook.img_comb = [];
default.get_heights.imgs = {[2*ones(8,1),(9:16).']};
default.combine.imgs = default.get_heights.imgs;
default.combine.img_comb = default.get_heights.qlook.img_comb;
default.radar.DC_adjust = {'DC_20160413_04_wf1.mat','DC_20160413_04_wf2.mat','DC_20160413_04_wf2.mat'};
default.radar.ref_fn = 'deconv_wf_%w_adc_%a_20160413_06_image3';
default.xml_regexp = 'image3_180-210MHz_.*thick.xml';
default.name = '3 Beam Image Mode 180-210 MHz';
defaults{end+1} = default;

%% Other settings

default.get_heights.qlook.img_comb = [];
default.get_heights.imgs = [];
default.combine.imgs = default.get_heights.imgs;
default.combine.img_comb = default.get_heights.qlook.img_comb;
default.radar.DC_adjust = [];
default.radar.ref_fn = '';
default.xml_regexp = '.*';
default.name = 'Other Settings';
defaults{end+1} = default;

return;
