function [param,defaults] = default_radar_params_2018_Antarctica_Ground_mcords6
% [param,defaults] = default_radar_params_2018_Antarctica_Ground_mcords6
%
% MCoRDS6: 2018_Antarctica_Ground
%
% Creates base "param" struct
% Creates defaults cell array for each type of radar setting
%
% Author: John Paden

param.season_name = '2018_Antarctica_Ground';
param.radar_name = 'mcords6';

%% Control parameters (not used in the parameter spreadsheet directly)
default.header_load_func = @basic_load_arena;

default.noise_50ohm = [0 0 0 0 0 0 0 0];

default.Pt = 400; % Transmit power at the transmit antenna
default.Gt = 4*2; % Transmit antenna gain
default.Ae = default.Gt*(3e8/195e6)^2; % Receiver antenna effective area
default.system_loss_dB = 10.^(-5.88/10); % Losses from the receive antenna to before the first LNA
default.noise_figure = 2; % Noise figure of receiver starting at the first LNA
default.adc_SNR_dB = 59; % ADC full scale signal SNR (relative to quantization noise)
default.fs = 640e6;
default.fs_dac = 1280e6;
default.max_duty_cycle = 0.1;
default.max_data_rate = 60;
default.max_tx = [1 1 1 1];
default.PRI_guard = 1e-6;
default.PRI_guard_percentage = 450e6/500e6;

default.tx_enable = [1 1 1 1];

default.basic_surf_track_min_time = 2e-6;
default.basic_surf_track_Tpd_factor = 1.1; % Normally -inf for lab test, 1.1 for flight test
default.board_map = {'digrx0','digrx1','digrx2','digrx3'};
default.tx_map = {'awg0','awg1','awg2','awg3'};
default.records.file.boards = [1 2 3 4];
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
  default.txequal.time_delay_desired = [0];
  default.txequal.phase_desired = [0];
  default.txequal.time_validation = [0.4]*1e-9;
  default.txequal.amp_validation = [3];
  default.txequal.phase_validation = [35];
  default.txequal.remove_linear_phase_en = true;
end

%% MCoRDS6 Arena Parameters
arena = [];
subsystem_idx = 0;
subsystem_idx = subsystem_idx + 1;
arena.subsystem(subsystem_idx).name = 'ARENA-CTU';
arena.subsystem(subsystem_idx).subSystem{1} = 'ctu';
subsystem_idx = subsystem_idx + 1;
arena.subsystem(subsystem_idx).name = 'ARENA0';
arena.subsystem(subsystem_idx).subSystem{1} = 'awg0';
arena.subsystem(subsystem_idx).subSystem{2} = 'digrx0';
subsystem_idx = subsystem_idx + 1;
arena.subsystem(subsystem_idx).name = 'ARENA1';
arena.subsystem(subsystem_idx).subSystem{1} = 'awg1';
arena.subsystem(subsystem_idx).subSystem{2} = 'digrx1';
subsystem_idx = subsystem_idx + 1;
arena.subsystem(subsystem_idx).name = 'ARENA2';
arena.subsystem(subsystem_idx).subSystem{1} = 'awg2';
arena.subsystem(subsystem_idx).subSystem{2} = 'digrx2';
subsystem_idx = subsystem_idx + 1;
arena.subsystem(subsystem_idx).name = 'ARENA3';
arena.subsystem(subsystem_idx).subSystem{1} = 'awg3';
arena.subsystem(subsystem_idx).subSystem{2} = 'digrx3';
subsystem_idx = subsystem_idx + 1;
arena.subsystem(subsystem_idx).name = 'Data Server';

dac_idx = 0;
dac_idx = dac_idx + 1;
arena.dac(dac_idx).name = 'awg0';
arena.dac(dac_idx).type = 'dac-ad9129_0012';
arena.dac(dac_idx).dacClk = default.fs_dac;
arena.dac(dac_idx).desiredAlignMin = -24;
arena.dac(dac_idx).desiredAlignMax = -10;
arena.dac(dac_idx).dcoPhase = 80;
dac_idx = dac_idx + 1;
arena.dac(dac_idx).name = 'awg1';
arena.dac(dac_idx).type = 'dac-ad9129_0012';
arena.dac(dac_idx).dacClk = default.fs_dac;
arena.dac(dac_idx).desiredAlignMin = -14;
arena.dac(dac_idx).desiredAlignMax = 0;
arena.dac(dac_idx).dcoPhase = 80;
dac_idx = dac_idx + 1;
arena.dac(dac_idx).name = 'awg2';
arena.dac(dac_idx).type = 'dac-ad9129_0012';
arena.dac(dac_idx).dacClk = default.fs_dac;
arena.dac(dac_idx).desiredAlignMin = -17;
arena.dac(dac_idx).desiredAlignMax = -3;
arena.dac(dac_idx).dcoPhase = 80;
dac_idx = dac_idx + 1;
arena.dac(dac_idx).name = 'awg3';
arena.dac(dac_idx).type = 'dac-ad9129_0012';
arena.dac(dac_idx).dacClk = default.fs_dac;
arena.dac(dac_idx).desiredAlignMin = 0;
arena.dac(dac_idx).desiredAlignMax = 14;
arena.dac(dac_idx).dcoPhase = 80;

adc_idx = 0;
adc_idx = adc_idx + 1;
arena.adc(adc_idx).name = 'digrx0';
arena.adc(adc_idx).type = 'adc-ad9680_0017';
arena.adc(adc_idx).sampFreq = default.fs;
arena.adc(adc_idx).adcMode = 2;
arena.adc(adc_idx).desiredAlignMin = 0;
arena.adc(adc_idx).desiredAlignMax = 14;
arena.adc(adc_idx).ip = '10.0.0.100';
arena.adc(adc_idx).outputSelect = 0;
adc_idx = adc_idx + 1;
arena.adc(adc_idx).name = 'digrx1';
arena.adc(adc_idx).type = 'adc-ad9680_0017';
arena.adc(adc_idx).sampFreq = default.fs;
arena.adc(adc_idx).adcMode = 2;
arena.adc(adc_idx).desiredAlignMin = -7;
arena.adc(adc_idx).desiredAlignMax = 10;
arena.adc(adc_idx).ip = '10.0.0.100';
arena.adc(adc_idx).outputSelect = 0;
adc_idx = adc_idx + 1;
arena.adc(adc_idx).name = 'digrx2';
arena.adc(adc_idx).type = 'adc-ad9680_0017';
arena.adc(adc_idx).sampFreq = default.fs;
arena.adc(adc_idx).adcMode = 2;
arena.adc(adc_idx).desiredAlignMin = -10;
arena.adc(adc_idx).desiredAlignMax = 4;
arena.adc(adc_idx).ip = '10.0.0.100';
arena.adc(adc_idx).outputSelect = 0;
adc_idx = adc_idx + 1;
arena.adc(adc_idx).name = 'digrx3';
arena.adc(adc_idx).type = 'adc-ad9680_0017';
arena.adc(adc_idx).sampFreq = default.fs;
arena.adc(adc_idx).adcMode = 2;
arena.adc(adc_idx).desiredAlignMin = 13;
arena.adc(adc_idx).desiredAlignMax = 27;
arena.adc(adc_idx).ip = '10.0.0.100';
arena.adc(adc_idx).outputSelect = 0;

daq_idx = 0;
daq_idx = daq_idx + 1;
arena.daq(daq_idx).name = 'daq0';
arena.daq(daq_idx).type = 'daq_0001';
arena.daq(daq_idx).auxDir = '/data/';
arena.daq(daq_idx).fileStripe = '/data/%b/';
arena.daq(daq_idx).fileName = 'mcords';

arena.system.name = 'ku0001';
arena.param.board_map = default.board_map;
arena.param.tx_map = default.tx_map;

arena.param.tx_max = [1 1 1 1];
arena.param.PA_setup_time = 2e-6; % Time required to enable PA before transmit
arena.param.TTL_time_delay = 0.0; % TTL time delay relative to transmit start
arena.param.ADC_time_delay = 0.0; % ADC time delay relative to transmit start
   
% mode 0, subchannel 0, board_idx 1 is wf-adc 1-1
% mode 0, subchannel 0, board_idx 2 is wf-adc 2-1

arena.psc.type = 'psc_0003';

arena.daq.type = 'daq_0001';

arena.ctu.name = 'ctu';
arena.ctu.type = 'ctu_001D';
arena.ctu.nmea = 31;
arena.ctu.pps = 10;
arena.ctu.pps_polarity = 1;
idx = 0;
idx = idx + 1;
arena.ctu.out.bit_group(idx).name = 'EPRI';
arena.ctu.out.bit_group(idx).bits = 0;
arena.ctu.out.bit_group(idx).epri = [1 0];
arena.ctu.out.bit_group(idx).pri = [0 0];
idx = idx + 1;
arena.ctu.out.bit_group(idx).name = 'PRI';
arena.ctu.out.bit_group(idx).bits = 1;
arena.ctu.out.bit_group(idx).epri = [1 0];
arena.ctu.out.bit_group(idx).pri = [1 0];
idx = idx + 1;
arena.ctu.out.bit_group(idx).name = 'TR';
arena.ctu.out.bit_group(idx).bits = 2;
arena.ctu.out.bit_group(idx).epri = [1 0];
arena.ctu.out.bit_group(idx).pri = [1 0];
idx = idx + 1;
arena.ctu.out.bit_group(idx).name = 'CalSwitch';
arena.ctu.out.bit_group(idx).bits = 3;
arena.ctu.out.bit_group(idx).epri = [1 0];
arena.ctu.out.bit_group(idx).pri = [1 0];
idx = idx + 1;
arena.ctu.out.bit_group(idx).name = 'Atten';
arena.ctu.out.bit_group(idx).bits = 4:8;
arena.ctu.out.bit_group(idx).epri = [0 0];
arena.ctu.out.bit_group(idx).pri = [0 0];
%   idx = idx + 1;
%   arena.ctu.out.bit_group(idx).name = 'AttenFirst18dB';
%   arena.ctu.out.bit_group(idx).bits = 4;
%   arena.ctu.out.bit_group(idx).epri = [1 1];
%   arena.ctu.out.bit_group(idx).pri = [1 1];
%   idx = idx + 1;
%   arena.ctu.out.bit_group(idx).name = 'AttenSecond7dB';
%   arena.ctu.out.bit_group(idx).bits = 5;
%   arena.ctu.out.bit_group(idx).epri = [0 0];
%   arena.ctu.out.bit_group(idx).pri = [0 0];

arena.ctu.out.time_cmd = {'2e-6+param.wfs(wf).Tpd+0.5e-6' '2/param.prf'};

default.arena = arena;

%% Records worksheet in parameter spreadsheet
default.records.gps.time_offset = 0;
default.records.frames.geotiff_fn = 'antarctica/Landsat-7/Antarctica_LIMA_480m.tif';
default.records.frames.mode = 1;

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
default.radar.Tadc_adjust = 8.3042e-06; % System time delay: leave this empty or set it to zero at first, determine this value later using data over surface with known height or from surface multiple
default.radar.lever_arm_fh = @lever_arm;
default.radar.adc_gains_dB = [45 45 45 45 45 45 45 45]; % Gain from the first LNA to the ADC
default.radar.rx_paths = [1 2 3 4 5 6 7 8];
chan_equal_Tsys = [0 0 0 0 0 0 0 0]/1e9;
chan_equal_dB = [0 0 0 0 0 0 0 0];
chan_equal_deg = [0 0 0 0 0 0 0 0];

defaults = {};

% Deconvolution Mode
default.records.data_map = {[2 0 1 1;2 1 1 2;5 0 2 1;5 1 2 2],[2 0 1 3;2 1 1 4;5 0 2 3;5 1 2 4],[2 0 1 5;2 1 1 6;5 0 2 5;5 1 2 6],[2 0 1 7;2 1 1 8;5 0 2 7;5 1 2 8]};
default.qlook.qlook.img_comb = [];
default.qlook.imgs = {[1*ones(8,1),(1:8).'],[2*ones(8,1),(1:8).']};
default.sar.imgs = default.qlook.imgs;
default.array.imgs = default.qlook.imgs;
default.array.img_comb = default.qlook.qlook.img_comb;
default.radar.ref_fn = '';
for wf = 1:2
  default.radar.wfs(wf).chan_equal_Tsys = chan_equal_Tsys;
  default.radar.wfs(wf).chan_equal_dB = chan_equal_dB;
  default.radar.wfs(wf).chan_equal_deg = chan_equal_deg;
  default.radar.wfs(wf).adcs = [1 2 3 4 5 6 7 8];
end

default.config_regexp = '.*deconv.*';
default.name = 'Deconv Mode 160-230 MHz';
defaults{end+1} = default;

% Survey Mode
default.records.data_map = {[2 0 1 1;2 1 1 2;5 0 2 1;5 1 2 2],[2 0 1 3;2 1 1 4;5 0 2 3;5 1 2 4],[2 0 1 5;2 1 1 6;5 0 2 5;5 1 2 6],[2 0 1 7;2 1 1 8;5 0 2 7;5 1 2 8]};
default.qlook.img_comb = [10e-06 -inf 3e-06];
default.qlook.imgs = {[1*ones(8,1),(1:8).'],[2*ones(8,1),(1:8).']};
default.sar.imgs = default.qlook.imgs;
default.array.imgs = default.qlook.imgs;
default.array.img_comb = default.qlook.qlook.img_comb;
default.radar.ref_fn = '';
for wf = 1:2
  default.radar.wfs(wf).chan_equal_Tsys = chan_equal_Tsys;
  default.radar.wfs(wf).chan_equal_dB = chan_equal_dB;
  default.radar.wfs(wf).chan_equal_deg = chan_equal_deg;
  default.radar.wfs(wf).adcs = [1 2 3 4 5 6 7 8];
end

default.config_regexp = '.*survey.*';
default.name = 'Survey Mode 160-230 MHz';
defaults{end+1} = default;

%% Other settings
default.records.data_map = {[2 0 1 1;2 1 1 2;5 0 2 1;5 1 2 2],[2 0 1 3;2 1 1 4;5 0 2 3;5 1 2 4],[2 0 1 5;2 1 1 6;5 0 2 5;5 1 2 6],[2 0 1 7;2 1 1 8;5 0 2 7;5 1 2 8]};
default.qlook.img_comb = [];
default.qlook.imgs = {[1*ones(8,1),(1:8).'],[2*ones(8,1),(1:8).']};
default.sar.imgs = default.qlook.imgs;
default.array.imgs = default.qlook.imgs;
default.array.img_comb = default.qlook.qlook.img_comb;
default.radar.ref_fn = '';
for wf = 1:2
  default.radar.wfs(wf).chan_equal_Tsys = chan_equal_Tsys;
  default.radar.wfs(wf).chan_equal_dB = chan_equal_dB;
  default.radar.wfs(wf).chan_equal_deg = chan_equal_deg;
  default.radar.wfs(wf).adcs = [1 2 3 4 5 6 7 8];
end

default.config_regexp = '.*';
default.name = 'Default 160-230 MHz';
defaults{end+1} = default;