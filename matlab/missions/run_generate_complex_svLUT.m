% script run_generate_complex_svLUT.m
%
% Example script for running generate_complex_svLUT.m
%
% Author: John Paden

if 1
  %% 2018_Greenland_P3
  % fn = output from coh_noise surf tracker
  param = read_param_xls(ct_filename_param('rds_param_2018_Greenland_P3.xls'),'20180315_09',{'analysis_equal' 'analysis'});
  output_fn = '/scratch/rds/2018_Greenland_P3/CSARP_analysis/sv_table_2018_Greenland_P3.mat';
  
  param.collate_equal.in_dir = 'analysis';
  param.collate_equal.motion_comp_en = true;
  param.collate_equal.cmd_idx = 1;
  param.collate_equal.chan_eq_en = true;
  param.collate_equal.zero_surf_bin = [];
  param.collate_equal.rlines = [];
  
  param.collate_equal.imgs = {[1]};
  param.collate_equal.wf_adc_idxs = {{[1:7]}};
  for wf = 1:8
    params = ct_set_params(params,sprintf('radar.wfs(%d).Tsys',wf),[0.46 -4.66 0.14 -1.77 0 -2.63 -3.38 -69.66 -75.57 -75.45 -80.42 -80.49 -75.71 -77.69 -70.53]/1e9);
    params = ct_set_params(params,sprintf('radar.wfs(%d).chan_equal_dB',wf),[6.8 -0.6 3 0.1 0 3.5 3.9 7 3.3 4.8 6.1 6.2 4.6 3.1 6.2]);
    params = ct_set_params(params,sprintf('radar.wfs(%d).chan_equal_deg',wf),[-166.2 -142.7 177 -95.9 0 -25.9 -86.5 -27.4 128.1 41.6 -46.8 43 90.7 121.3 31.6]);
  end
  
  % Transmit windowing
  Hchan = boxcar(7).';
  
  analysis.surf.motion_comp.en = true;
  analysis.surf.chan_eq.en = true;
  tx_ant = [1:7]; % wf-adc index for each antenna
  analysis.surf.rlines = [];
  
  ref_ant = 3; % Antenna to use as phase reference channel (usually a center element)
  
  good_mask_min_samples = 10;
  good_mask_min_angle = -55;
  good_mask_max_angle = 55;
  
  plot_min_angle = -41;
  plot_max_angle = 41;
  
  fit_method = 'filter'; % 'filter' or 'sgolay'
  degrees_of_freedom = 3; % order of spatial filter fitted to the data
  
  debug_level = 1;
  
elseif 0
  %% 2017_Greenland_P3
  % fn = output from coh_noise surf tracker
  fn = '/scratch/rds/2017_Greenland_P3/CSARP_noise/surf_20170226_01_img_01.mat';
  output_fn = '/scratch/rds/2017_Greenland_P3/CSARP_noise/sv_table_2017_Greenland_P3.mat';
  
  % Transmit windowing
  Hchan = boxcar(7).';
  
  analysis.surf.motion_comp.en = true;
  analysis.surf.chan_eq.en = true;
  tx_ant = 1:7; % wf-adc index for each antenna
  analysis.surf.rlines = [];
  
  ref_ant = 3; % Antenna to use as phase reference channel (usually a center element)
  
  good_mask_min_samples = 10;
  good_mask_min_angle = -45;
  good_mask_max_angle = 45;
  
  plot_min_angle = -41;
  plot_max_angle = 41;
  
  fit_method = 'filter'; % 'filter' or 'sgolay'
  degrees_of_freedom = 3; % order of spatial filter fitted to the data
  
  debug_level = 1;
  
elseif 0
  %% 2016_Antarctica_DC8
  % fn = output from coh_noise surf tracker
  fn = '/scratch/rds/2016_Antarctica_DC8/CSARP_noise/surf_20161004_08_img_01.mat';
  output_fn = '/scratch/rds/2016_Antarctica_DC8/CSARP_noise/sv_table_2016_Antarctica_DC8.mat';
  
  % Transmit windowing
  Hchan = boxcar(6).';
  
  analysis.surf.motion_comp.en = true;
  analysis.surf.chan_eq.en = true;
  tx_ant = 1:6; % wf-adc index for each antenna
  analysis.surf.rlines = [];
  
  ref_ant = 5; % Antenna to use as phase reference channel (usually a center element)
  
  good_mask_min_samples = 10;
  good_mask_min_angle = -45;
  good_mask_max_angle = 45;
  
  plot_min_angle = -41;
  plot_max_angle = 41;
  
  fit_method = 'filter'; % 'filter' or 'sgolay'
  degrees_of_freedom = 3; % order of spatial filter fitted to the data
  
  debug_level = 1;
  
elseif 0
  %% 2015_Greenland_C130
  % fn = output from coh_noise surf tracker
  fn = 'F:\rds\2015_Greenland_C130\CSARP_noise\surf_20150313_14_img_01.mat';
  output_fn = 'F:/rds/2015_Greenland_C130/CSARP_noise/sv_table_2015_Greenland_C130.mat';
  
  % Transmit windowing
  Hchan = [1 1];
  
  analysis.surf.motion_comp.en = true;
  analysis.surf.chan_eq.en = true;
  tx_ant = 1:2; % wf-adc index for each antenna
  analysis.surf.rlines = [];
  
  ref_ant = 1; % Antenna to use as phase reference channel (usually a center element)
  
  good_mask_min_samples = 10;
  good_mask_min_angle = -50;
  good_mask_max_angle = 40;
  
  plot_min_angle = -50;
  plot_max_angle = 40;
  
  fit_method = 'filter'; % 'filter' or 'sgolay'
  degrees_of_freedom = 3; % order of spatial filter fitted to the data
  
  debug_level = 1;
  
  
elseif 0
  %% 2013_Antarctica_Basler
  % fn = output from coh_noise surf tracker
  fn = 'F:\rds\2013_Antarctica_Basler\CSARP_noise\surf_20131216_05_img_01.mat';
  output_fn = 'F:/rds/2013_Antarctica_Basler/CSARP_noise/sv_table_2013_Antarctica_Basler.mat';
  
  % Transmit windowing
  Hchan = chebwin(8,30).'; % NOT SURE
  
  analysis.surf.motion_comp.en = true;
  analysis.surf.chan_eq.en = true;
  tx_ant = 1:8; % wf-adc index for each antenna
  analysis.surf.rlines = [];
  
  ref_ant = 4; % Antenna to use as phase reference channel (usually a center element)
  
  good_mask_min_samples = 10;
  good_mask_min_angle = -57;
  good_mask_max_angle = 51;
  
  plot_min_angle = -57;
  plot_max_angle = 51;
  
  fit_method = 'filter'; % 'filter' or 'sgolay'
  degrees_of_freedom = 3; % order of spatial filter fitted to the data
  
  debug_level = 1;
  
elseif 0
  %% 2015_Greenland_Polar6
  % fn = output from coh_noise surf tracker
  fn = 'F:\rds\2015_Greenland_Polar6\CSARP_noise\surf_20150911_15_img_01.mat';
  output_fn = 'F:/rds/2015_Greenland_Polar6/CSARP_noise/sv_table_2015_Greenland_Polar6.mat';
  
  % Transmit windowing
  Hchan = chebwin(8,30).';
  
  analysis.surf.motion_comp.en = true;
  analysis.surf.chan_eq.en = true;
  tx_ant = 1:8; % wf-adc index for each antenna
  analysis.surf.rlines = [];
  
  ref_ant = 4; % Antenna to use as phase reference channel (usually a center element)
  
  good_mask_min_samples = 8;
  good_mask_min_angle = -40;
  good_mask_max_angle = 40;
  
  plot_min_angle = -35;
  plot_max_angle = 35;
  
  fit_method = 'filter'; % 'filter' or 'sgolay'
  degrees_of_freedom = 3; % order of spatial filter fitted to the data
  
  debug_level = 1;
  
elseif 0
  %% 2016_Greenland_Polar6
  % fn = output from coh_noise surf tracker
  fn = 'F:\rds\2016_Greenland_Polar6\CSARP_noise\surf_20160401_13_img_01.mat';
  output_fn = 'F:/rds/2016_Greenland_Polar6/CSARP_noise/sv_table_2016_Greenland_Polar6.mat';
  
  % Transmit windowing
  Hchan = chebwin(8,30).';
  
  analysis.surf.motion_comp.en = true;
  analysis.surf.chan_eq.en = true;
  tx_ant = 1:8; % wf-adc index for each antenna
  analysis.surf.rlines = 1500:6500;
  
  ref_ant = 4; % Antenna to use as phase reference channel (usually a center element)
  
  good_mask_min_samples = 10;
  good_mask_min_angle = -30;
  good_mask_max_angle = 30;
  
  plot_min_angle = -25;
  plot_max_angle = 25;
  
  fit_method = 'filter'; % 'filter' or 'sgolay'
  degrees_of_freedom = 3; % order of spatial filter fitted to the data
  
  debug_level = 1;
  
elseif 0
  %% 2017_Antarctica_Polar6
  % fn = output from coh_noise surf tracker
  fn = '/cresis/snfs1/dataproducts/ct_data/rds/2017_Antarctica_Polar6/CSARP_noise/surf_20160830_03_img_01.mat';
  output_fn = '/cresis/snfs1/dataproducts/ct_data/rds/2017_Antarctica_Polar6/CSARP_noise/sv_table_2017_Antarctica_Polar6_1us.mat';
  
  % Transmit windowing
  Hchan = chebwin(8,30).';
  
  analysis.surf.motion_comp.en = true;
  analysis.surf.chan_eq.en = true;
  tx_ant = 1:8; % wf-adc index for each antenna
  analysis.surf.rlines = [];
  
  ref_ant = 4; % Antenna to use as phase reference channel (usually a center element)
  
  good_mask_min_samples = 10;
  good_mask_min_angle = -30;
  good_mask_max_angle = 30;
  
  plot_min_angle = -25;
  plot_max_angle = 25;
  
  fit_method = 'filter'; % 'filter' or 'sgolay'
  degrees_of_freedom = 3; % order of spatial filter fitted to the data
  
  debug_level = 1;
end

generate_complex_svLUT;

return;
