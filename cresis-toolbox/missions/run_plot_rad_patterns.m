% script run_plot_rad_patterns.m
%
% Example script for running plot_rad_patterns.m
%
% Author: John Paden

if 1
  %% 2016_Antarctica_DC8
  
  fn = '/scratch/rds/2016_Antarctica_DC8/CSARP_noise/surf_20161004_08_img_01.mat';
  elements_fn = '/scratch/rds/2016_Antarctica_DC8/CSARP_noise/sv_table_2016_Antarctica_DC8.mat';
  
  analysis.surf.motion_comp.en = true;
  analysis.surf.chan_eq.en = true;
  analysis.surf.rlines = [];
  
  good_mask_min_samples = 10;
  good_mask_min_angle = -49;
  good_mask_max_angle = 49; 
  
  plot_min_angle = -45;
  plot_max_angle = 45;
  
  fit_method = 'filter'; % 'filter' or 'sgolay'
  
  debug_level = 1;
  f0 = 175e6;
  f1 = 205e6;
    
  % Combined Pattern
%   rad_patterns = [13]; % wf-adc indexes for patterns you want to generate
%   ref_pattern = [1];
%   sv_ant_ref = [5];
%   retrack_en = false;
%   equalize_angle = 0;
%   output_fn = fullfile('/scratch/rds/2016_Antarctica_DC8/CSARP_noise/', ...
%     sprintf('combined_pattern_2016_Antarctica_DC8.mat'));
%   degrees_of_freedom = 3; % spatial filter fitting (usually around the # of antennas)
%   
%   plot_rad_patterns;
%   figure(6);
%   title(sprintf('%3.0f - %3.0f MHz',f0/1e6,f1/1e6));
%   legend off;
%   saveas(6,fullfile('/scratch/rds/2016_Antarctica_DC8/CSARP_noise/', ...
%     sprintf('combined_pattern_2016_Antarctica_DC8.fig')));
  
  % Individual Elements
  rad_patterns = [7:12]; % wf-adc indexes for patterns you want to generate
  ref_pattern = 5;
  sv_LUT_ref = 5*ones(1,6);
  retrack_en = false;
  equalize_angle = 0;
  output_fn = '/scratch/rds/2016_Antarctica_DC8/CSARP_noise/sv_table_all_2016_Antarctica_DC8.mat';
  degrees_of_freedom = 3*ones(1,6); % spatial filter fitting (usually around the # of antennas)
  plot_rad_patterns;
  
elseif 0
  %% 2015_Greenland_C130
  fn = 'D:\rds\2015_Greenland_C130\CSARP_noise\surf_20150313_14_img_01.mat';
  elements_fn = 'D:\rds\2015_Greenland_C130\CSARP_noise\sv_table_2015_Greenland_C130.mat';
  
  analysis.surf.motion_comp.en = true;
  analysis.surf.chan_eq.en = true;
  analysis.surf.rlines = [];
  
  good_mask_min_samples = 10;
  good_mask_min_angle = -50;
  good_mask_max_angle = 40; 
  
  plot_min_angle = -50;
  plot_max_angle = 40;
  
  fit_method = 'filter'; % 'filter' or 'sgolay'
  
  debug_level = 1;
  
  error('Need more setup variables')
  
  plot_rad_patterns;

elseif 0
  %% 2015_Greenland_Polar6
  fn = 'J:\rds\2015_Greenland_Polar6\CSARP_noise\surf_20150911_17_img_01.mat';
  elements_fn = 'J:\rds\2015_Greenland_Polar6\CSARP_noise\sv_table_2015_Greenland_Polar6.mat';
  
  analysis.surf.motion_comp.en = true;
  analysis.surf.chan_eq.en = true;
  analysis.surf.rlines = [];
  
  good_mask_min_samples = 10;
  good_mask_min_angle = -35;
  good_mask_max_angle = 35; 
  
  plot_min_angle = -30;
  plot_max_angle = 30;
  
  fit_method = 'filter'; % 'filter' or 'sgolay'
  
  debug_level = 1;
  
  error('Need more setup variables')
  
  plot_rad_patterns;

elseif 0
  %% 2016_Greenland_Polar6
  fn = 'F:\rds\2016_Greenland_Polar6\CSARP_noise\surf_20160401_13_img_01.mat';
  fn = '/home/administrator/Scratch/rds/2016_Greenland_Polar6/CSARP_noise/surf_20160414_13_img_01.mat';
  elements_fn = '/home/administrator/Scratch/rds/2016_Greenland_Polar6/CSARP_noise/sv_table_2016_Greenland_Polar6.mat';
  
  analysis.surf.motion_comp.en = true;
  analysis.surf.chan_eq.en = true;
  analysis.surf.rlines = [];
  
  good_mask_min_samples = 10;
  good_mask_min_angle = -49;
  good_mask_max_angle = 49; 
  
  plot_min_angle = -45;
  plot_max_angle = 45;
  
  fit_method = 'filter'; % 'filter' or 'sgolay'
  
  debug_level = 1;
  
  for f0 = [168.5e6 : 50e6 : 450e6]
    f1 = f0 + 50e6;
    
    % Combined Pattern
    % rad_patterns = [33]; % wf-adc indexes for patterns you want to generate
    % ref_pattern = [1];
    % sv_ant_ref = [4];
    % retrack_en = false;
    % equalize_angle = 0;
    % output_fn = fullfile('F:\rds\2016_Greenland_Polar6\CSARP_noise\', ...
    %   sprintf('combined_pattern_2016_Greenland_Polar6_%3.0f_%3.0fMHz.mat',f0/1e6,f1/1e6));
    % degrees_of_freedom = round(16*(f0+f1)/450e6/2) % spatial filter fitting (usually around the # of antennas)
    % 
    % plot_rad_patterns;
    % figure(6);
    % title(sprintf('%3.0f - %3.0f MHz',f0/1e6,f1/1e6));
    % legend off;
    % saveas(6,fullfile('F:\rds\2016_Greenland_Polar6\CSARP_noise\', ...
    %   sprintf('combined_pattern_2016_Greenland_Polar6_%3.0f_%3.0fMHz.fig',f0/1e6,f1/1e6)));
    
    % Left Pattern
    rad_patterns = [1]; % wf-adc indexes for patterns you want to generate
    ref_pattern = [1];
    sv_ant_ref = [4];
    retrack_en = true;
    equalize_angle = -20;
    output_fn = fullfile('/home/administrator/Scratch/rds/2016_Greenland_Polar6/CSARP_noise/', ...
      sprintf('left_pattern_2016_Greenland_Polar6_%3.0f_%3.0fMHz.mat',f0/1e6,f1/1e6));
    degrees_of_freedom = round(16*(f0+f1)/450e6/2) % spatial filter fitting (usually around the # of antennas)

    plot_rad_patterns;
    figure(6);
    title(sprintf('%3.0f - %3.0f MHz',f0/1e6,f1/1e6));
    legend off;
    saveas(6,fullfile('/home/administrator/Scratch/rds/2016_Greenland_Polar6/CSARP_noise/', ...
      sprintf('left_pattern_2016_Greenland_Polar6_%3.0f_%3.0fMHz.fig',f0/1e6,f1/1e6)));
    
    % Right Pattern
%     rad_patterns = [3]; % wf-adc indexes for patterns you want to generate
%     ref_pattern = [1];
%     sv_ant_ref = [4];
%     equalize_angle = 20;
%     retrack_en = true;
%     output_fn = fullfile('/home/administrator/Scratch/rds/2016_Greenland_Polar6/CSARP_noise/', ...
%       sprintf('right_pattern_2016_Greenland_Polar6_%3.0f_%3.0fMHz.mat',f0/1e6,f1/1e6));
%     degrees_of_freedom = round(16*(f0+f1)/450e6/2) % spatial filter fitting (usually around the # of antennas)
%     
%     plot_rad_patterns;
%     figure(6);
%     title(sprintf('%3.0f - %3.0f MHz',f0/1e6,f1/1e6));
%     legend off;
%     saveas(6,fullfile('/home/administrator/Scratch/rds/2016_Greenland_Polar6/CSARP_noise/', ...
%       sprintf('right_pattern_2016_Greenland_Polar6_%3.0f_%3.0fMHz.fig',f0/1e6,f1/1e6)));
    
    % Individual Elements (Center Elements Only)
    %   rad_patterns = [17:24]; % wf-adc indexes for patterns you want to generate
    %   %ref_pattern = [repmat(5,[1 8]) repmat(13,[1 8]) repmat(21,[1 8])];
    %   ref_pattern = repmat(4,[1 8]);
    %   sv_ant_ref = repmat(4,[1 8]);
    % retrack_en = false;
    % equalize_angle = 0;
    %   output_fn = 'F:\rds\2016_Greenland_Polar6\CSARP_noise\sv_table_all_2016_Greenland_Polar6.mat';
    %   degrees_of_freedom = repmat(3,[1 8]); % spatial filter fitting (usually around the # of antennas)
    
    % Individual Elements
    %   rad_patterns = [9:30, 32, 31]; % wf-adc indexes for patterns you want to generate
    %   %ref_pattern = [repmat(5,[1 8]) repmat(13,[1 8]) repmat(21,[1 8])];
    %   ref_pattern = repmat(12,[1 24]);
    %   sv_ant_ref = repmat(4,[1 24]);
    % retrack_en = false;
    % equalize_angle = 0;
    %   output_fn = 'F:\rds\2016_Greenland_Polar6\CSARP_noise\sv_table_all_2016_Greenland_Polar6.mat';
    %   degrees_of_freedom = repmat(3,[1 24]); % spatial filter fitting (usually around the # of antennas)
    
  end
  

elseif 1
  %% 2017_Antarctica_Polar6
  fn = '/cresis/snfs1/dataproducts/ct_data/rds/2017_Antarctica_Polar6/CSARP_noise/surf_20160830_03_img_01.mat';
  elements_fn = '/cresis/snfs1/dataproducts/ct_data/rds/2017_Antarctica_Polar6/CSARP_noise/sv_table_2017_Antarctica_Polar6.mat';
  
  analysis.surf.motion_comp.en = true;
  analysis.surf.chan_eq.en = true;
  analysis.surf.rlines = [];
  
  good_mask_min_samples = 10;
  good_mask_min_angle = -30;
  good_mask_max_angle = 30; 
  
  plot_min_angle = -25;
  plot_max_angle = 25;
  
  fit_method = 'filter'; % 'filter' or 'sgolay'
  
  debug_level = 1;
  
  for f0 = 150e6%[168.5e6 : 50e6 : 450e6]
    f1 = f0 + 50e6;
    f1 = 520e6;
    
    % Combined Pattern
    % rad_patterns = [33]; % wf-adc indexes for patterns you want to generate
    % ref_pattern = [1];
    % sv_ant_ref = [4];
    % retrack_en = false;
    % equalize_angle = 0;
    % output_fn = fullfile('F:\rds\2017_Antarctica_Polar6\CSARP_noise\', ...
    %   sprintf('combined_pattern_2017_Antarctica_Polar6_%3.0f_%3.0fMHz.mat',f0/1e6,f1/1e6));
    % degrees_of_freedom = round(16*(f0+f1)/450e6/2) % spatial filter fitting (usually around the # of antennas)
    % 
    % plot_rad_patterns;
    % figure(6);
    % title(sprintf('%3.0f - %3.0f MHz',f0/1e6,f1/1e6));
    % legend off;
    % saveas(6,fullfile('F:\rds\2017_Antarctica_Polar6\CSARP_noise\', ...
    %   sprintf('combined_pattern_2017_Antarctica_Polar6_%3.0f_%3.0fMHz.fig',f0/1e6,f1/1e6)));
    
    % Left Pattern
%     rad_patterns = [1]; % wf-adc indexes for patterns you want to generate
%     ref_pattern = [1];
%     sv_ant_ref = [4];
%     retrack_en = true;
%     equalize_angle = -20;
%     output_fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2017_Antarctica_Polar6/CSARP_noise/', ...
%       sprintf('left_pattern_2017_Antarctica_Polar6_%3.0f_%3.0fMHz.mat',f0/1e6,f1/1e6));
%     degrees_of_freedom = round(16*(f0+f1)/450e6/2) % spatial filter fitting (usually around the # of antennas)
% 
%     plot_rad_patterns;
%     figure(6);
%     title(sprintf('%3.0f - %3.0f MHz',f0/1e6,f1/1e6));
%     legend off;
%     saveas(6,fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2017_Antarctica_Polar6/CSARP_noise/', ...
%       sprintf('left_pattern_2017_Antarctica_Polar6_%3.0f_%3.0fMHz.fig',f0/1e6,f1/1e6)));
    
    % Right Pattern
%     rad_patterns = [3]; % wf-adc indexes for patterns you want to generate
%     ref_pattern = [1];
%     sv_ant_ref = [4];
%     equalize_angle = 20;
%     retrack_en = true;
%     output_fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2017_Antarctica_Polar6/CSARP_noise/', ...
%       sprintf('right_pattern_2017_Antarctica_Polar6_%3.0f_%3.0fMHz.mat',f0/1e6,f1/1e6));
%     degrees_of_freedom = round(16*(f0+f1)/450e6/2) % spatial filter fitting (usually around the # of antennas)
%     
%     plot_rad_patterns;
%     figure(6);
%     title(sprintf('%3.0f - %3.0f MHz',f0/1e6,f1/1e6));
%     legend off;
%     saveas(6,fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2017_Antarctica_Polar6/CSARP_noise/', ...
%       sprintf('right_pattern_2017_Antarctica_Polar6_%3.0f_%3.0fMHz.fig',f0/1e6,f1/1e6)));
    
    % Individual Elements (Center Elements Only)
    %   rad_patterns = [17:24]; % wf-adc indexes for patterns you want to generate
    %   %ref_pattern = [repmat(5,[1 8]) repmat(13,[1 8]) repmat(21,[1 8])];
    %   ref_pattern = repmat(4,[1 8]);
    %   sv_ant_ref = repmat(4,[1 8]);
    % retrack_en = false;
    % equalize_angle = 0;
    %   output_fn = '/cresis/snfs1/dataproducts/ct_data/rds/2017_Antarctica_Polar6/CSARP_noise/sv_table_all_2017_Antarctica_Polar6.mat';
    %   degrees_of_freedom = repmat(3,[1 8]); % spatial filter fitting (usually around the # of antennas)
    
    % Individual Elements
    rad_patterns = [9:30, 32, 31]; % wf-adc indexes for patterns you want to generate
    ref_pattern = [repmat(5,[1 8]) repmat(13,[1 8]) repmat(21,[1 8])];
%     ref_pattern = repmat(12,[1 24]);
    sv_ant_ref = repmat(4,[1 24]);
    
%     rad_patterns = [17:24]; % wf-adc indexes for patterns you want to generate
%     ref_pattern = repmat(4,[1 8]);
%     sv_ant_ref = repmat(4,[1 8]);
    retrack_en = true;
    equalize_angle = 0;
    output_fn = '/cresis/snfs1/dataproducts/ct_data/rds/2017_Antarctica_Polar6/CSARP_noise/sv_table_all_2017_Antarctica_Polar6.mat';
    degrees_of_freedom = repmat(3,[1 24]); % spatial filter fitting (usually around the # of antennas)
    
    plot_rad_patterns;
    figure(6);
    title(sprintf('%3.0f - %3.0f MHz',f0/1e6,f1/1e6));
    legend off;
    saveas(6,fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2017_Antarctica_Polar6/CSARP_noise/', ...
      sprintf('individual_2017_Antarctica_Polar6_%3.0f_%3.0fMHz.fig',f0/1e6,f1/1e6)));
    
  end
  
end

return;

