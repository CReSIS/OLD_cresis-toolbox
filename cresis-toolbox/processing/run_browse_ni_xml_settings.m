% script run_browse_ni_xml_settings
%
% Example script for running browse_ni_xml_settings.m
%
% Author: John Paden

clear % REMOVE THIS LINE!!!

base_dir = 'Y:\paden\awi\';
base_dir_in_param = '/cresis/snfs1/scratch/paden/awi/';

adc_folder_names = {''};

param.radar_name = 'mcords5';
param.season_name = '2016_Greenland_Polar6';

header_load_date = datenum(2016,3,11);
default.cmd.mission_name = 'Lab Test';

param_file.write_en = true;
param_file.path = ''; % Set to override the default path

%% MCORDS 5: 2015 Greenland Polar6
param.radar_name = 'mcords5';
xml_file_prefix = 'mcords5_20160311';
data_file_prefix = 'mcords5';
board_sub_directory = 'chan1';

default.records.geotiff_fn = 'greenland/Landsat-7/Greenland_natural_150m'
default.records.file.adcs = [1:24];
default.records.file.adc_headers = [1:24];
default.records.gps.en = 1;
default.records.frame_mode = 0;
default.records.presum_bug_fixed = 1;
default.records.tmp_fn_uses_adc_folder_name = 1;

%default.get_heights.qlook.out_path = '';
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
default.get_heights.surf.en = 0;
default.get_heights.surf.method = 'threshold';
default.get_heights.surf.noise_rng = [0 -50 10];
default.get_heights.surf.min_bin = 2e-6;
default.get_heights.surf.max_bin = [];
default.get_heights.surf.threshold = 9;
default.get_heights.surf.sidelobe = 15;
default.get_heights.surf.medfilt = 3;
default.get_heights.surf.search_rng = [0:2];

default.csarp.out_path = '';
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

header_load_func = @basic_load_mcords5;
header_load_params.clk = 1600e6;
adc_bits = 12;
adc_full_scale = 2; 
rx_paths = [1:24];
chan_equal_dB = '[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
chan_equal_deg = '[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
chan_equal_Tsys = '[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]/1e9';
Tadc = []; % normally leave empty to use value in file header
Tadc_adjust = 0.000010179163; % leave this empty or set it to zero at first, determine this value later using data over surface with known height or from surface multiple
max_DDS_RAM = 4000;
tx_voltage = sqrt(1000*50)*10^(-2/20);
max_adc_gain_dB = 52;
iq_mode = 0;
fs = 1600e6;
fs_sync = 1.6e9/8;
tx_DDS_mask = [1 1 1 1 1 1 1 1];
imgs_adcs = 1:8;
radar_worksheet_headers = {'Tpd','Tadc','Tadc_adjust','f0','f1','ref_fn','tukey','tx_weights','rx_paths','adc_gains','chan_equal_dB','chan_equal_deg','Tsys','DDC_mode','DDC_freq'};
radar_worksheet_headers_type = {'r','r','r','r','r','r','r','r','r','r','r','r','r','r','r'};

%% MCORDS 5: 2015 Greenland C130
% param.radar_name = 'mcords5';
% xml_file_prefix = 'mcords5';
% data_file_prefix = 'mcords5';
% board_sub_directory = 'chan1';
% header_load_func = @basic_load_mcords5;
% header_load_params.clk = 1600e6;
% adc_bits = 12;
% adc_full_scale = 2;
% rx_paths = [1:2];
% chan_equal_dB = '[0 0]';
% chan_equal_deg = '[0 0]';
% chan_equal_Tsys = '[0 0]/1e9';
% Tadc = []; % normally leave empty to use value in file header
% Tadc_adjust = 9.92e-6; % leave this empty or set it to zero at first, determine this value later using data over surface with known height or from surface multiple
% max_DDS_RAM = 3000;
% tx_voltage = sqrt(1000*50)*10^(-2/20);
% iq_mode = 0;
% fs = 1600e6;
% fs_sync = 1.6e9/8;
% tx_DDS_mask = [1 1];
% imgs_adcs = 1:2;
% radar_worksheet_headers = {'Tpd','Tadc','Tadc_adjust','f0','f1','ref_fn','tukey','tx_weights','rx_paths','adc_gains','chan_equal_dB','chan_equal_deg','Tsys','DDC_mode','DDC_freq'};
% radar_worksheet_headers_type = {'r','r','r','r','r','r','r','r','r','r','r','r','r','r','r'};

%% MCORDS 4: 2013 Antarctica Basler
% param.radar_name = 'mcords4';
% xml_file_prefix = 'mcords4';
% data_file_prefix = 'mcords4';
% board_sub_directory = 'chan2';
% header_load_func = @basic_load_mcords4;
% header_load_params.clk = 1e9/2;
% adc_bits = 12;
% adc_full_scale = 2;
% rx_paths = [1 3 5 7 2 4 6 8];
% chan_equal_dB = '[-0.9 0 0.3 0.4 0.2 0.4 0.2 -1]';
% chan_equal_deg = '[12.6 19.3 27.3 35.2 28.5 28 9.4 6.3]';
% chan_equal_Tsys = '[82.06 82.00 81.88 82.10 82.08 81.97 81.90 81.77]/1e9';
% Tadc_adjust = 0; % leave this empty or set it to zero first, you can later determine this value from surface multiple.
% max_DDS_RAM = 60000;
% tx_voltage = sqrt(250*50)*10^(-2/20);
% iq_mode = -j;
% fs = 1e9/2;
% fs_sync = 1e9/8;
% tx_DDS_mask = [1 1 1 1 1 1 1 1];
% imgs_adcs = 1:8;
% % radar_worksheet_headers = {'Tpd','Tadc','Tadc_adjust','f0','f1','ref_fn','tukey','tx_weights','rx_paths','adc_gains','chan_equal_dB','chan_equal_deg','Tsys'};
% % radar_worksheet_headers_type = {'r','r','r','r','r','r','r','r','r','r','r','r','r'};
% radar_worksheet_headers = {'Tpd','Tadc','f0','f1','ref_fn','tukey','tx_weights','rx_paths','adc_gains','chan_equal_dB','chan_equal_deg','Tsys'};
% radar_worksheet_headers_type = {'r','r','r','r','r','r','r','r','r','r','r','r'};

%% MCORDS 3: 2014 Greenland P3
% param.radar_name = 'mcords3';
% data_file_prefix = 'mcords3';
% board_sub_directory = 'board0';
% header_load_func = @basic_load_mcords3;
% header_load_params.clk = 1e9/9;
% adc_bits = 14;
% adc_full_scale = 2;
% rx_paths = [1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
% chan_equal_dB = '[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
% chan_equal_deg = '[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
% chan_equal_Tsys = '([0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]+1640)/1e9';
% Tadc = []; % normally leave empty to use value in file header
% Tadc_adjust = []; % leave this empty or set it to zero at first, determine this value later using data over surface with known height or from surface multiple
% max_DDS_RAM = 30000;
% tx_voltage = sqrt(300*50)*10^(-2/20);
% iq_mode = 0;
% fs = 1e9/9;
% fs_sync = 1e9/18;
% tx_DDS_mask = [1 1 1 1 1 1 1 0];
% imgs_adcs = 2:15;
% radar_worksheet_headers = {'Tpd','Tadc','Tadc_adjust','f0','f1','ref_fn','tukey','tx_weights','rx_paths','adc_gains','chan_equal_dB','chan_equal_deg','Tsys'};
% radar_worksheet_headers_type = {'r','r','r','r','r','r','r','r','r','r','r','r','r'};

%% MCORDS 3: 2014 Antarctica DC8
% param.radar_name = 'mcords3';
% data_file_prefix = 'mcords3';
% board_sub_directory = 'board0';
% header_load_func = @basic_load_mcords3;
% header_load_params.clk = 300e6/2;
% adc_bits = 14;
% adc_full_scale = 2;
% %rx_paths = [1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
% rx_paths = [1 2 3 4 5 6 1 1];
% chan_equal_dB = '[-0.3 -3.7 -1.9  0 -4.7  0.30]';
% chan_equal_deg = '[-11.4 -33.9 -44.6 0 -39.8 -48.5]';
% chan_equal_Tsys = '[0 0 0 0 0 0]/1e9';
% Tadc_adjust = 1.6e-6; % leave this empty or set it to zero first, you can later determine this value from surface multiple.
% max_DDS_RAM = 30000;
% tx_voltage = sqrt(300*50)*10^(-2/20);
% iq_mode = 0;
% %fs = 1e9/9;
% %fs_sync = 1e9/18;
% fs = 300e6/2;
% fs_sync = 300e6/4;
% tx_DDS_mask = [1 1 1 1 1 1];
% imgs_adcs = [1 2 3 4 5 6];
% radar_worksheet_headers = {'Tpd','Tadc','Tadc_adjust','f0','f1','ref_fn','tukey','tx_weights','rx_paths','adc_gains','chan_equal_dB','chan_equal_deg','Tsys'};
% radar_worksheet_headers_type = {'r','r','r','r','r','r','r','r','r','r','r','r','r'};

%% MCORDS 2
% param.radar_name = 'mcords2';
% data_file_prefix = 'mcords2';
% board_sub_directory = 'board0';
% header_load_func = @basic_load_mcords2;
% header_load_params.clk = 1e9/9;
% adc_bits = 14;
% adc_full_scale = 2;
% rx_paths = [1 1 2 3 4 5 6 7];
% chan_equal_dB = '[0 0 0 0 0 0 0 0]';
% chan_equal_deg = '[0 0 0 0 0 0 0 0]';
% chan_equal_Tsys = '[790 790 790 790 790 790 790 790]/1e9';
% max_DDS_RAM = 30000;
% tx_voltage = sqrt(300*50)*10^(-2/20);
% iq_mode = 0;
% fs = 1e9/9;
% tx_DDS_mask = [1 1 1 1 1 1 1 0];
% imgs_adcs = 2:15;
% radar_worksheet_headers = {'Tpd','Tadc','Tadc_adjust','f0','f1','ref_fn','tukey','tx_weights','rx_paths','adc_gains','chan_equal_dB','chan_equal_deg','Tsys'};
% radar_worksheet_headers_type = {'r','r','r','r','r','r','r','r','r','r','r','r','r'};

%% 
% xml_version = There are three versions of XML files
%   mcords3: pre-2014 Greenland P3 is 1
%   mcords3: post-2014 Greenland P3 is 3
%   mcords4: 2
%   mcords5: 2 (2016 Polar 6) or 4 (e.g. 2015 LC130)
%   (THIS OFTEN NEEDS TO BE SET)
xml_version = 2;

% gps_fn: Set to plot GPS locations on a map (leave empty to disable) and
%   also to plot roll data
gps_fn = '/cresis/snfs1/dataproducts/csarp_support/gps/2013_Antarctica_Basler/gps_20131216.mat';
gps_fn = '';

if ~isempty(gps_fn)
  gps = load(gps_fn);
end

MIN_FILES_IN_SEGMENT = 3; % Since the first file is usually corrupted and short segments are often not useful for processing/data-interpretation
SKIP_FILES_START = 1; % Often the first file is corrupt, so we leave it out
SKIP_FILES_END = 0;

% geotiff_fn: Set to a geotiff file to plot GPS location on (leave empty to
%   disable)
geotiff_fn = ct_filename_gis(param,'antarctica/Landsat-7/Antarctica_LIMA_480m.tif');
geotiff_fn = ct_filename_gis(param,'greenland/Landsat-7/mzl7geo_90m_lzw.tif');
geotiff_fn = '';

% =======================================================================
% =======================================================================
%% Automated Section
% =======================================================================
% =======================================================================

browse_ni_xml_settings;
