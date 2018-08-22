function [param,defaults] = default_radar_params_2018_Antarctica_TObas
% [param,defaults] = default_radar_params_2018_Antarctica_TObas
%
% Accum3: 2018_Antarctica_TObas
%
% Creates base "param" struct
% Creates defaults cell array for each type of radar setting
%
% Author: John Paden

param.season_name = '2018_Antarctica_TObas';
param.radar_name = 'accum3';

%% Control parameters (not used in the parameter spreadsheet directly)
default.header_load_func = @basic_load_arena;
default.header_load_params = struct('clk',1600e6,'presum_bug_fixed',true);
default.xml_version = 100.0;

default.noise_50ohm = [0 0 0 0];

default.Pt = 400;
default.Gt = 4*2;
default.Ae = default.Gt*(3e8/750e6)^2;

default.system_loss_dB = 10.^(-5.88/10);
default.tx_voltage = sqrt(1000*50)*10^(-2/20);

default.iq_mode = 0;
default.tx_DDS_mask = [1 1 1 1];

default.basic_surf_track_min_time = 2e-6;
default.basic_surf_track_Tpd_factor = 1.1; % Normally -inf for lab test, 1.1 for flight test
default.records.file.board_folder_name = 'chan%b';
default.records.file.boards = [1 2];
default.records.file.version = 103;

if 1
  % Example 1: Normal configuration:
  %   Connect antenna N to WFG N for all N = 1 to 8
  ref_adc = 1;
  default.txequal.img = [(1:1).', ref_adc*ones(1,1)];
  default.txequal.ref_wf_adc = 1;
  default.txequal.wf_mapping = [1];
  default.txequal.Hwindow_desired = chebwin(1,30).';
  default.txequal.max_DDS_amp = [4000];
  default.txequal.time_delay_desired = [0 0 0 0 ];
  default.txequal.phase_desired = [0 0 0 0 ];
  default.txequal.time_validation = [0.4 0.4 0.4 0.4 ]*1e-9;
  default.txequal.amp_validation = [3 3 3 3 ];
  default.txequal.phase_validation = [35 35 35 35 ];
  default.txequal.remove_linear_phase_en = true;
end

%% Records worksheet in parameter spreadsheet
default.records.gps.time_offset = 0;
default.records.frames.geotiff_fn = 'antarctica/Landsat-7/Antarctica_LIMA_480m.tif';
default.records.frames.mode = 0;

%% Quick Look worksheet in parameter spreadsheet
default.qlook.out_path = '';
default.qlook.block_size = 5000;
default.qlook.motion_comp = 0;
default.qlook.dec = 20;
default.qlook.inc_dec = 10;
default.qlook.surf.en = 1;
default.qlook.surf.min_bin = 2e-6;
default.qlook.surf.method = 'threshold';
default.qlook.surf.threshold = 17;
default.qlook.surf.filter_len = 7;
default.qlook.surf.sidelobe = 17;
default.qlook.surf.noise_rng = [0 -50 10];
default.qlook.surf.search_rng = [0:2];

%% SAR worksheet in parameter spreadsheet
default.sar.out_path = '';
default.sar.imgs = {[1*ones(4,1),(1:4).'],[2*ones(4,1),(1:4).']};
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
default.radar.adc_bits = 14;
default.radar.Vpp_scale = 2;
default.radar.noise_figure = 2;
default.radar.adc_SNR_dB = 59;
default.radar.Tadc_adjust = 8.3042e-06; % System time delay: leave this empty or set it to zero at first, determine this value later using data over surface with known height or from surface multiple
default.radar.lever_arm_fh = @lever_arm;
default.radar.adc_gains = 10.^([45 27]/20);
default.radar.rx_paths = [1 1];

defaults = {};

% Deconvolution Mode
default.records.arena.data_map = {[0 0 1 1],[0 0 2 1]};
default.qlook.qlook.img_comb = [];
default.qlook.imgs = {[1*ones(1,1),(1:1).'],[2*ones(1,1),(1:1).']};
default.sar.imgs = default.qlook.imgs;
default.array.imgs = default.qlook.imgs;
default.array.img_comb = default.qlook.qlook.img_comb;
default.radar.ref_fn = '';
default.radar.wfs(1).chan_equal_Tsys = [0]/1e9;
default.radar.wfs(1).chan_equal_dB = [0];
default.radar.wfs(1).chan_equal_deg = [0];

default.xml_regexp = '.*deconv.xml';
default.name = 'Deconv Mode 600-900 MHz';
defaults{end+1} = default;

% Survey Mode
default.records.arena.data_map = {[0 0 1 1],[0 0 2 1]};
default.qlook.qlook.img_comb = [2e-06 -inf 2e-06];
default.qlook.imgs = {[1*ones(1,1),(1:1).'],[2*ones(1,1),(1:1).']};
default.sar.imgs = default.qlook.imgs;
default.array.imgs = default.qlook.imgs;
default.array.img_comb = default.qlook.qlook.img_comb;
default.radar.ref_fn = '';
default.radar.wfs(1).chan_equal_Tsys = [0]/1e9;
default.radar.wfs(1).chan_equal_dB = [0];
default.radar.wfs(1).chan_equal_deg = [0];

default.xml_regexp = 'survey_600-900MHz_.*.xml';
default.name = 'Survey Mode 600-900 MHz';
defaults{end+1} = default;

%% Other settings
default.records.arena.data_map = {[0 0 1 1],[0 0 2 1]};
default.qlook.qlook.img_comb = [];
default.qlook.imgs = {[1*ones(1,1),(1:1).'],[2*ones(1,1),(1:1).']};
default.sar.imgs = default.qlook.imgs;
default.array.imgs = default.qlook.imgs;
default.array.img_comb = default.qlook.qlook.img_comb;
default.radar.ref_fn = '';
default.radar.wfs(1).chan_equal_Tsys = [0]/1e9;
default.radar.wfs(1).chan_equal_dB = [0];
default.radar.wfs(1).chan_equal_deg = [0];

default.xml_regexp = '.*';
default.name = 'Default 600-900 MHz';
defaults{end+1} = default;
