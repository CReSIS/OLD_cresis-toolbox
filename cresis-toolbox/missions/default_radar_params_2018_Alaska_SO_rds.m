function param = default_radar_params_2018_Alaska_SO_rds
% param = default_radar_params_2018_Alaska_SO_rds
%
% RDS: 2018_Alaska_SO
%
% Creates base "param" struct
% Creates defaults cell array for each type of radar setting
%
% Author: John Paden

%% Preprocess parameters
param.season_name = '2018_Alaska_SO';
param.radar_name = 'rds';

param.config.field_time_gap = 'gps_time';
param.config.file.version = 413;
param.config.file.prefix = '';
param.config.file.suffix = '.tdms';
param.config.max_time_gap = 10;
param.config.min_seg_size = 1;

param.config.daq_type = 'utua';
param.config.wg_type = 'utua';
param.config.header_load_func = [];
param.config.board_map = {''};
param.config.tx_map = {'','','','','','','',''};

param.config.max_data_rate = 100;
param.config.max_duty_cycle = 0.12;
param.config.prf_multiple = []; % Power supply sync signal that PRF must be a factor of these numbers
param.config.PRI_guard = 10e-6;
param.config.PRI_guard_percentage = 1;
param.config.tx_enable = [1];
param.config.max_tx = 1;
param.config.max_tx_voltage = sqrt(250*50)*10^(-2/20); % voltage at max_tx

%% UTUA parameters
param.config.utua.clk = 1e9/9;
param.config.utua.rx_gain_dB = 51.5;

%% Command worksheet
default.cmd.records = 1;
default.cmd.qlook = 1;
default.cmd.generic = 1;

%% Records worksheet
default.records.gps.time_offset = 1;
default.records.file.adcs = [1 2];
default.records.frames.mode = 1;
default.records.frames.geotiff_fn = fullfile('alaska','Landsat-7','Alaska_90m.tif');

%% Qlook worksheet
default.qlook.out_path = '';
default.qlook.block_size = 50000;
default.qlook.dec = 20;
default.qlook.inc_dec = 5;
default.qlook.surf.en = 1;
default.qlook.surf.profile = 'RDS_OIB';

%% SAR worksheet
default.sar.out_path = '';
default.sar.imgs = {[1 1]};
default.sar.frm_types = {0,[0 1],0,0,-1};
default.sar.chunk_len = 5000;
default.sar.combine_rx = 0;
default.sar.time_of_full_support = 3.5e-5;
default.sar.mocomp.en = 1;
default.sar.mocomp.type = 2;
default.sar.mocomp.filter = {@butter  [2]  [0.1000]};
default.sar.mocomp.uniform_en = 1;
default.sar.sar_type = 'fk';
default.sar.sigma_x = 60;
default.sar.sub_aperture_steering = 0;
default.sar.st_wind = @hanning;
default.sar.start_eps = 3.15;

%% Array worksheet
default.array.in_path = '';
default.array.out_path = '';
default.array.img_comb = [];
default.array.imgs = {[1 1]};
default.array.method = 'standard';
default.array.window = @hanning;
default.array.bin_rng = 0;
default.array.line_rng = -1:1;
default.array.dbin = 1;
default.array.dline = 1;

%% Radar worksheet
default.radar.fs = 50e6;
default.radar.Tadc = []; % normally leave empty to use value in file header
default.radar.adc_bits = 14;
default.radar.Vpp_scale = 2;

default.radar.wfs.rx_paths = [1 1];
default.radar.wfs.noise_figure = 2;
default.radar.wfs.Tadc_adjust = 1.2e-06; % System time delay: leave this empty or set it to zero at first, determine this value later using data over surface with known height or from surface multiple

%% Post worksheet
default.post.data_dirs = {'qlook'};
default.post.layer_dir = 'layerData';
default.post.maps_en = 1;
default.post.echo_en = 1;
default.post.layers_en = 1;
default.post.data_en = 0;
default.post.csv_en = 1;
default.post.concat_en = 1;
default.post.pdf_en = 1;
default.post.map.location = 'Custom';
default.post.map.type = 'combined';
default.post.echo.elev_comp = 3;
default.post.echo.depth = '[publish_echogram_switch(Bbad,0.25,Surface_Elev,-3500,DBottom,-100),max(Surface_Elev+100)]';
default.post.echo.er_ice = 3.15;
default.post.ops.en = 0;
default.post.ops.location = 'arctic';
default.post.ops.layers = {'bottom','surface'};
default.post.ops.gaps_dist = [300 60];


%% Radar Settings
defaults = {};

default.radar.wfs(1).Tsys = [0]/1e9;
default.radar.wfs(1).chan_equal_dB = [0];
default.radar.wfs(1).chan_equal_deg = [0];

% survey mode
default.qlook.img_comb = [];
default.qlook.imgs = {[1 1]};
default.sar.imgs = default.qlook.imgs;
default.array.imgs = default.qlook.imgs;
default.array.img_comb = default.qlook.img_comb;
default.radar.ref_fn = '';
default.config_regexp = '';
default.name = '2.5 MHz Survey Mode';
defaults{end+1} = default;

% survey mode
default.qlook.img_comb = [];
default.qlook.imgs = {[1 1]};
default.sar.imgs = default.qlook.imgs;
default.array.imgs = default.qlook.imgs;
default.array.img_comb = default.qlook.img_comb;
default.config_regexp = '';
default.name = '5 MHz Survey Mode';
defaults{end+1} = default;

% Other settings

default.qlook.img_comb = [];
default.qlook.imgs = [];
default.sar.imgs = default.qlook.imgs;
default.array.imgs = default.qlook.imgs;
default.array.img_comb = default.qlook.img_comb;
default.config_regexp = '';
default.name = 'Other Settings';
defaults{end+1} = default;

%% Add default settings

param.config.defaults = defaults;
