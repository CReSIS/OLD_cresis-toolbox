% script run_generate_complex_svLUT.m
%
% Example script for running generate_complex_svLUT.m
%
% Author: John Paden

if 0
  %% 2014_Antarctica_DC8
  % fn = output from coh_noise surf tracker
  
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
  
  good_mask_min_samples = 10;
  good_mask_min_angle = -40;
  good_mask_max_angle = 40;
  
  plot_min_angle = -35;
  plot_max_angle = 35;
  
  fit_method = 'filter'; % 'filter' or 'sgolay'
  degrees_of_freedom = 3; % order of spatial filter fitted to the data
  
  debug_level = 1;
  
elseif 1
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
end

generate_complex_svLUT;

return;
