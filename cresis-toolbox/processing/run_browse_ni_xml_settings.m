% script run_browse_ni_xml_settings
%
% Example script for running browse_ni_xml_settings.m
%
% Author: John Paden

clear;
setup = 'mcords3_p3_setup';

%% MCoRDS5 Polar 6 Setup
% =======================================================================
if strcmpi(setup, 'mcords5_polar6_setup')

  base_dir = '/cresis/snfs1/scratch/2016_Germany_AWI_tests/AWI_ICE_bak/test_flight/';
  
  adc_folder_names = {'chan%d'};
  
  header_load_date = datenum(2016,8,30);
  
  mission_name = 'Testflight';

  [param,defaults] = default_radar_params_2017_Antarctica_Polar6_mcords;

  base_dir_in_param = base_dir; % Only modify if you want params to have a different path than the current path
    
  param_file.write_en = false; % Windows only and experimental still (allows direct writes to param spreadsheets)
  param_file.path = ''; % Set to override the default path to the parameter file
    
  % gps_fn: Set to plot GPS locations on a map (leave empty to disable) and
  %   also to plot roll data
  gps_fn = ct_filename_support(setfield(param,'day_seg',datestr(header_load_date,'YYYYmmDD')),'','gps',true);
  %gps_fn = '';
    
  MIN_FILES_IN_SEGMENT = 2; % Since the first file is usually corrupted and short segments are often not useful for processing/data-interpretation
  SKIP_FILES_START = 0; % Often the first file is corrupt, leave it out by setting this to 1
  SKIP_FILES_END = 0;
  
  % geotiff_fn: Set to a geotiff file to plot GPS location on (leave empty to
  %   disable)
  geotiff_fn = ct_filename_gis(param,'antarctica/Landsat-7/Antarctica_LIMA_480m.tif');
  geotiff_fn = ct_filename_gis(param,'greenland/Landsat-7/mzl7geo_90m_lzw.tif');
  %geotiff_fn = '';
  
  manual_enable = true;
  
  browse_ni_xml_settings;
  return;
end

%% MCoRDS3 P3 Setup
% =======================================================================
if strcmpi(setup, 'mcords3_p3_setup')

  base_dir = '/process-archive/20170226/mcords';
  
  adc_folder_names = {'board%b'};
  
  header_load_date = datenum(2017,2,26);
  
  mission_name = 'Testflight';

  [param,defaults] = default_radar_params_2017_Greenland_P3_mcords;

  base_dir_in_param = base_dir; % Only modify if you want params to have a different path than the current path
    
  param_file.write_en = false; % Windows only and experimental still (allows direct writes to param spreadsheets)
  param_file.path = ''; % Set to override the default path to the parameter file
    
  % gps_fn: Set to plot GPS locations on a map (leave empty to disable) and
  %   also to plot roll data
  gps_fn = ct_filename_support(setfield(param,'day_seg',datestr(header_load_date,'YYYYmmDD')),'','gps',true);
  %gps_fn = '';
    
  MIN_FILES_IN_SEGMENT = 2; % Since the first file is usually corrupted and short segments are often not useful for processing/data-interpretation
  SKIP_FILES_START = 0; % Often the first file is corrupt, leave it out by setting this to 1
  SKIP_FILES_END = 0;
  
  % geotiff_fn: Set to a geotiff file to plot GPS location on (leave empty to
  %   disable)
%   geotiff_fn = ct_filename_gis(param,'antarctica/Landsat-7/Antarctica_LIMA_480m.tif');
  geotiff_fn = ct_filename_gis(param,'greenland/Landsat-7/mzl7geo_90m_lzw.tif');
  %geotiff_fn = '';
  
  manual_enable = true;
  
  browse_ni_xml_settings;
  return;
end

%% MCoRDS3 DC8 Setup
% =======================================================================
if strcmpi(setup, 'mcords3_dc8_setup')

  base_dir = '/process/mcords';
  
  adc_folder_names = {'board%b'};
  
  header_load_date = datenum(2016,10,20);
  
  mission_name = 'Getz A';

  [param,defaults] = default_radar_params_2016_Antarctica_DC8_mcords;

  base_dir_in_param = base_dir; % Only modify if you want params to have a different path than the current path
    
  param_file.write_en = false; % Windows only and experimental still (allows direct writes to param spreadsheets)
  param_file.path = ''; % Set to override the default path to the parameter file
    
  % gps_fn: Set to plot GPS locations on a map (leave empty to disable) and
  %   also to plot roll data
  gps_fn = ct_filename_support(setfield(param,'day_seg',datestr(header_load_date,'YYYYmmDD')),'','gps',true);
  %gps_fn = '';
    
  MIN_FILES_IN_SEGMENT = 2; % Since the first file is usually corrupted and short segments are often not useful for processing/data-interpretation
  SKIP_FILES_START = 0; % Often the first file is corrupt, leave it out by setting this to 1
  SKIP_FILES_END = 0;
  
  % geotiff_fn: Set to a geotiff file to plot GPS location on (leave empty to
  %   disable)
  geotiff_fn = ct_filename_gis(param,'antarctica/Landsat-7/Antarctica_LIMA_480m.tif');
%   geotiff_fn = ct_filename_gis(param,'greenland/Landsat-7/mzl7geo_90m_lzw.tif');
  %geotiff_fn = '';
  
  manual_enable = true;
  
  browse_ni_xml_settings;
  return;
end

if strcmpi(setup, 'mcords4_setup')
  % =======================================================================
  % =======================================================================
  %% MCoRDS4 Setup
  % =======================================================================
  % =======================================================================
  
  if 0
    base_dir = '/cresis/snfs1/data/MCoRDS/2013_Antarctica_Basler/';
    base_dir_in_param = '/cresis/snfs1/data/MCoRDS/2013_Antarctica_Basler/';
    adc_folder_names = {'20140102'};
    param.season_name = '2013_Antarctica_Basler';
    header_load_date = datenum(2014,1,2);
    default.cmd.mission_name = 'Mission 1';
    xml_file_prefix = 'mcords4';
    board_sub_directory = 'chan1';
    xml_version = 1.5;
    param.radar_name = 'mcords4';
    data_file_prefix = 'mcords4';
    
    header_load_func = @basic_load_mcords4;
    header_load_params.clk = 500e6;
    iq_mode = 1;
    fs = 500e6;
    fs_sync = 1e9/16;
    radar_worksheet_headers = {'Tpd','Tadc','Tadc_adjust','f0','f1','ref_fn','tukey','tx_weights','rx_paths','adc_gains','chan_equal_dB','chan_equal_deg','Tsys'};
    radar_worksheet_headers_type = {'r','r','r','r','r','r','r','r','r','r','r','r','r'};
  else
    base_dir = '/cresis/snfs1/data/MCoRDS/2015_Greenland_C130/';
    base_dir_in_param = '/cresis/snfs1/data/MCoRDS/2015_Greenland_C130/';
    adc_folder_names = {'20150420'};
    param.season_name = '2015_Greenland_C130';
    header_load_date = datenum(2015,4,20);
    default.cmd.mission_name = 'Mission 1';
    xml_file_prefix = 'mcords5';
    board_sub_directory = 'chan1';
    xml_version = 1.6;
    param.radar_name = 'mcords5';
    data_file_prefix = 'mcords5';
    
    header_load_func = @basic_load_mcords5;
    header_load_params.clk = 1600e6;
    iq_mode = 0;
    fs = 1600e6;
    fs_sync = 1e9/16;
    radar_worksheet_headers = {'Tpd','Tadc','Tadc_adjust','f0','f1','ref_fn','tukey','tx_weights','rx_paths','adc_gains','chan_equal_dB','chan_equal_deg','Tsys','DDC_mode','DDC_freq'};
    radar_worksheet_headers_type = {'r','r','r','r','r','r','r','r','r','r','r','r','r','r','r'};
  end
  
  param_file.write_en = false;
  param_file.path = ''; % Set to override the default path
    
  if ~isempty(regexp(param.season_name,'Greenland'))
    default.records.geotiff_fn = 'greenland/Landsat-7/Greenland_natural_150m'
  elseif ~isempty(regexp(param.season_name,'Antarctica'))
    default.records.geotiff_fn = 'greenland/Landsat-7/Antarctica_LIMA_480m'
  end
  default.records.file.adcs = [1:16];
  default.records.file.adc_headers = [1:16];
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
  
  adc_bits = 12;
  adc_full_scale = 2;
  rx_paths = [1,1:15];
  chan_equal_dB = '[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]';
  chan_equal_deg = '[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
  chan_equal_Tsys = '[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]/1e9';
  Tadc = []; % normally leave empty to use value in file header
  Tadc_adjust = 0.000010179163; % leave this empty or set it to zero at first, determine this value later using data over surface with known height or from surface multiple
  max_DDS_RAM = 30000;
  tx_voltage = sqrt(300*50)*10^(-2/20);
  max_adc_gain_dB = 52;
  TTL_prog_delay = 650;
  tx_DDS_mask = [1 1 1 1 1 1 1 0];
  imgs_adcs = 2:15;
  
  % gps_fn: Set to plot GPS locations on a map (leave empty to disable) and
  %   also to plot roll data
  gps_fn = ct_filename_support(setfield(param,'day_seg',datestr(header_load_date,'YYYYmmDD')),'','gps',true);
  gps_fn = '';
  
  if ~isempty(gps_fn)
    gps = load(gps_fn);
  end
  
  MIN_FILES_IN_SEGMENT = 2; % Since the first file is usually corrupted and short segments are often not useful for processing/data-interpretation
  SKIP_FILES_START = 1; % Often the first file is corrupt, so we leave it out
  SKIP_FILES_END = 0;
  
  % geotiff_fn: Set to a geotiff file to plot GPS location on (leave empty to
  %   disable)
  geotiff_fn = ct_filename_gis(param,'antarctica/Landsat-7/Antarctica_LIMA_480m.tif');
  geotiff_fn = ct_filename_gis(param,'greenland/Landsat-7/Greenland_natural_150m.tif');
  geotiff_fn = '';

  browse_ni_xml_settings;
  return;

end

if strcmpi(setup, 'mcords2_setup')
  % =======================================================================
  % =======================================================================
  %% MCoRDS2 Setup
  % =======================================================================
  % =======================================================================
  
  if 0
    base_dir = '/cresis/snfs1/data/MCoRDS/2011_Greenland_P3/';
    base_dir_in_param = '/cresis/snfs1/data/MCoRDS/2011_Greenland_P3/';
    adc_folder_names = {'20110406'};
    param.season_name = '2011_Greenland_P3';
    header_load_date = datenum(2011,4,6);
    default.cmd.mission_name = 'Mission 1';
    xml_file_prefix = 'DDS_20110406';
    board_sub_directory = 'board0/seg_01';
    xml_version = 1.0;
    param.radar_name = 'mcords2';
    data_file_prefix = 'mcords2';
  elseif 0
    base_dir = '/cresis/snfs1/data/MCoRDS/2011_Antarctica_TO/';
    base_dir_in_param = '/cresis/snfs1/data/MCoRDS/2011_Antarctica_TO/';
    adc_folder_names = {'20111202'};
    param.season_name = '2011_Antarctica_TO';
    header_load_date = datenum(2011,12,2);
    default.cmd.mission_name = 'Mission 1';
    xml_file_prefix = 'DDS_20111202';
    board_sub_directory = 'board0';
    xml_version = 1.1;
    param.radar_name = 'mcords2';
    data_file_prefix = 'mcords2';
  elseif 0
    base_dir = '/cresis/snfs1/data/MCoRDS/2012_Greenland_P3/';
    base_dir_in_param = '/cresis/snfs1/data/MCoRDS/2012_Greenland_P3/';
    adc_folder_names = {'20120410'};
    param.season_name = '2012_Greenland_P3';
    header_load_date = datenum(2012,4,10);
    default.cmd.mission_name = 'Mission 1';
    xml_file_prefix = 'DDS_20120410_095222';
    board_sub_directory = 'board0/seg_01';
    xml_version = 1.1;
    param.radar_name = 'mcords2';
    data_file_prefix = 'mcords2';
  elseif 0
    base_dir = '/cresis/snfs1/data/MCoRDS/2012_Antarctica_DC8/';
    base_dir_in_param = '/cresis/snfs1/data/MCoRDS/2012_Antarctica_DC8/';
    adc_folder_names = {'20121015'};
    param.season_name = '2012_Antarctica_DC8';
    header_load_date = datenum(2012,10,15);
    default.cmd.mission_name = 'Mission 1';
    xml_file_prefix = 'DDS_20121015_160701';
    board_sub_directory = 'board0';
    xml_version = 1.2;
    param.radar_name = 'mcords2';
    data_file_prefix = 'mcords2';
  elseif 0
    base_dir = '/cresis/snfs1/data/MCoRDS/2013_Greenland_P3/';
    base_dir_in_param = '/cresis/snfs1/data/MCoRDS/2013_Greenland_P3/';
    adc_folder_names = {'20130426'};
    param.season_name = '2013_Greenland_P3';
    header_load_date = datenum(2013,4,26);
    default.cmd.mission_name = 'Mission 1';
    xml_file_prefix = 'DDS_20130426';
    board_sub_directory = 'board0';
    xml_version = 1.3;
    param.radar_name = 'mcords3';
    data_file_prefix = 'mcords3';
  elseif 0
    base_dir = '/cresis/snfs1/data/MCoRDS/2013_Antarctica_P3/';
    base_dir_in_param = '/cresis/snfs1/data/MCoRDS/2013_Antarctica_P3/';
    adc_folder_names = {'20131127'};
    param.season_name = '2013_Antarctica_P3';
    header_load_date = datenum(2013,11,27);
    default.cmd.mission_name = 'Mission 1';
    xml_file_prefix = 'DDS_2013112';
    board_sub_directory = 'board0';
    xml_version = 1.3;
    param.radar_name = 'mcords3';
    data_file_prefix = 'mcords3';
  elseif 1
    base_dir = '/cresis/snfs1/data/MCoRDS/2014_Greenland_P3/';
    base_dir_in_param = '/cresis/snfs1/data/MCoRDS/2014_Greenland_P3/';
    adc_folder_names = {'20140325'};
    param.season_name = '2014_Greenland_P3';
    header_load_date = datenum(2014,3,25);
    default.cmd.mission_name = 'Mission 1';
    xml_file_prefix = 'radar_20140325';
    board_sub_directory = 'board0';
    xml_version = 1.4;
    param.radar_name = 'mcords3';
    data_file_prefix = 'mcords3';
  elseif 0
    base_dir = '/cresis/snfs1/data/MCoRDS/2014_Antarctica_DC8/';
    base_dir_in_param = '/cresis/snfs1/data/MCoRDS/2014_Antarctica_DC8/';
    adc_folder_names = {'20141023'};
    param.season_name = '2014_Antarctica_DC8';
    header_load_date = datenum(2014,10,24);
    default.cmd.mission_name = 'Mission 1';
    xml_file_prefix = 'radar_20141024';
    board_sub_directory = 'board0';
    xml_version = 1.4;
    param.radar_name = 'mcords3';
    data_file_prefix = 'mcords3';
  end
  
  param_file.write_en = false;
  param_file.path = ''; % Set to override the default path
    
  if ~isempty(regexp(param.season_name,'Greenland'))
    default.records.geotiff_fn = 'greenland/Landsat-7/Greenland_natural_150m'
  elseif ~isempty(regexp(param.season_name,'Antarctica'))
    default.records.geotiff_fn = 'greenland/Landsat-7/Antarctica_LIMA_480m'
  end
  default.records.file.adcs = [1:16];
  default.records.file.adc_headers = [1:16];
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
  
  header_load_func = @basic_load_mcords2;
  header_load_params.clk = 1e9/9;
  adc_bits = 14;
  adc_full_scale = 2;
  rx_paths = [1,1:15];
  chan_equal_dB = '[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]';
  chan_equal_deg = '[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
  chan_equal_Tsys = '[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]/1e9';
  Tadc = []; % normally leave empty to use value in file header
  Tadc_adjust = 0.000010179163; % leave this empty or set it to zero at first, determine this value later using data over surface with known height or from surface multiple
  max_DDS_RAM = 30000;
  tx_voltage = sqrt(300*50)*10^(-2/20);
  max_adc_gain_dB = 52;
  iq_mode = 0;
  fs = 1e9/9;
  fs_sync = 1e9/18;
  TTL_prog_delay = 650;
  tx_DDS_mask = [1 1 1 1 1 1 1 0];
  imgs_adcs = 2:15;
  radar_worksheet_headers = {'Tpd','Tadc','Tadc_adjust','f0','f1','ref_fn','tukey','tx_weights','rx_paths','adc_gains','chan_equal_dB','chan_equal_deg','Tsys'};
  radar_worksheet_headers_type = {'r','r','r','r','r','r','r','r','r','r','r','r','r'};
  
  % gps_fn: Set to plot GPS locations on a map (leave empty to disable) and
  %   also to plot roll data
  gps_fn = ct_filename_support(setfield(param,'day_seg',datestr(header_load_date,'YYYYmmDD')),'','gps',true);
  gps_fn = '';
  
  if ~isempty(gps_fn)
    gps = load(gps_fn);
  end
  
  MIN_FILES_IN_SEGMENT = 2; % Since the first file is usually corrupted and short segments are often not useful for processing/data-interpretation
  SKIP_FILES_START = 1; % Often the first file is corrupt, so we leave it out
  SKIP_FILES_END = 0;
  
  % geotiff_fn: Set to a geotiff file to plot GPS location on (leave empty to
  %   disable)
  geotiff_fn = ct_filename_gis(param,'antarctica/Landsat-7/Antarctica_LIMA_480m.tif');
  geotiff_fn = ct_filename_gis(param,'greenland/Landsat-7/Greenland_natural_150m.tif');
  geotiff_fn = '';

  browse_ni_xml_settings;
  return;

end
