clearvars -except gRadar
param = [];

param = [];

%% Petermann Line 1 2014
% if ispc
%   fn = fullfile('X:\ct_data\rds\2014_Greenland_P3\CSARP_multipass\',sprintf('Petermann_line1_2014.mat'));
% else
%   fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_multipass/',sprintf('Petermann_line1_2014'));
% end

%% Petermann Line 1 2011, 2014, 2018
% if ispc
%   param.multipass.fn = fullfile('X:\ct_data\rds\2014_Greenland_P3\CSARP_multipass\',sprintf('Petermann_line1_2011_2014_2018.mat'));
% else
%   param.multipass.fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_multipass/',sprintf('Petermann_line1_2011_2014_2018'));
% end
% 
% param.multipass.rbins = [];
% 
% param.multipass.baseline_master_idx = 2;
% param.multipass.master_idx = 2;
% 
% param.multipass.pass_en_mask = [];
% param.multipass.output_fn_midfix = [];
% param.multipass.coregistration_time_shift = [0 0 -2];
% param.multipass.comp_mode = 2;


%% Petermann Line 2 2013, 2014
% if ispc
%   param.multipass.fn = fullfile('X:\ct_data\rds\2014_Greenland_P3\CSARP_multipass\',sprintf('Petermann_line2_2013_2014.mat'));
% else
%   param.multipass.fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_multipass/',sprintf('Petermann_line2_2013_2014'));
% end
% 
% param.multipass.rbins = [];
% 
% param.multipass.baseline_master_idx = 2;
% param.multipass.master_idx = 2;
% 
% param.multipass.pass_en_mask = [];
% param.multipass.output_fn_midfix = [];
% param.multipass.coregistration_time_shift = [-0.5 0];
% param.multipass.comp_mode = 2;

%% Petermann Line 4 2010, 2011, 2013, 2014
% if ispc
%   param.multipass.fn = fullfile('X:\ct_data\rds\2014_Greenland_P3\CSARP_multipass\',sprintf('Petermann_line4_2010_2011_2013_2014.mat'));
% else
%   param.multipass.fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_multipass/',sprintf('Petermann_line4_2010_2011_2013_2014'));
% end
% 
% param.multipass.rbins = [];
% 
% param.multipass.baseline_master_idx = 2;
% param.multipass.master_idx = 2;
% 
% param.multipass.pass_en_mask = [];
% param.multipass.output_fn_midfix = [];
% param.multipass.coregistration_time_shift = [0 -0.5 0 0];
% param.multipass.comp_mode = 2;


%% 79N Line 1 2010, 2014, 2016, 2018
% if ispc
%   param.multipass.fn = fullfile('X:\ct_data\rds\2014_Greenland_P3\CSARP_multipass\',sprintf('79N_line1_2010_2014_2016_2018.mat'));
% else
%   param.multipass.fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_multipass/',sprintf('79N_line1_2010_2014_2016_2018'));
% end
% 
% param.multipass.rbins = [];
% 
% param.multipass.baseline_master_idx = 2;
% param.multipass.master_idx = 2;
% 
% param.multipass.pass_en_mask = [];
% param.multipass.output_fn_midfix = [];
% param.multipass.coregistration_time_shift = [1 0 0 -2];
% param.multipass.time_gate = [2e-6 13e-6];
% 
% comp_mode = 2;
%% 2011 to 2012 Greenland P3
% radartype = 'rds';
% passname = 'rds_thule_2011_2012';
% param.multipass.fn = fullfile(gRadar.out_path,radartype,'2014_Greenland_P3','CSARP_multipass',passname);
% 
% param.multipass.rbins = 220:420;
% 
% param.multipass.baseline_master_idx = 8;
% param.multipass.master_idx = 8;
% 
% param.multipass.pass_en_mask = [];
% param.multipass.output_fn_midfix = [];
% param.multipass.coregistration_time_shift = [];
% param.multipass.comp_mode = 1;
% param.multipass.time_gate = [];

%% 2014 Same Day Greenland P3
if 1
  radartype = 'rds';
  passname = 'rds_thule_2014_SameDay_allwf';
  param.multipass.fn = fullfile(gRadar.out_path,radartype,'2014_Greenland_P3','CSARP_multipass',passname);
  
  param.multipass.rbins = [];
  
  param.multipass.baseline_master_idx = 8;
  param.multipass.master_idx = 8;
  
  param.multipass.pass_en_mask = [];
  param.multipass.output_fn_midfix = [];
  param.multipass.coregistration_time_shift = [];
  param.multipass.time_gate = [];
  
  %Load equalization vector (sar specific)
  eqvec1 = [122.5 121.6 126.7 105.5 130.1 120.1 128.7 -0.0 -134.2 121.1 36.4 125.8 -171.5 -1.0 128.4];
  eqvec2 = eqvec1;
  neweq = [42.8 45.5 48.2 50.9 53.6 56.3 58.9 0.0 2.7 5.4 8.1 93.6 96.3 99.0 101.7]/2;
  neweq = [neweq neweq];

  if 1 && exist('neweq','var') && ~isempty(neweq)
    eqvec1 = eqvec1-neweq(1:length(eqvec1));
    eqvec2 = eqvec2-neweq(length(eqvec1)+1:end);
  end

  equalization1 = 10.^(zeros(1,15)/20) .* exp(1i*(eqvec1)/180*pi);
  equalization2 = 10.^(zeros(1,15)/20) .* exp(1i*(eqvec2)/180*pi);
  param.multipass.equalization = [equalization1 equalization2];
  
  param.multipass.debug_plots = {'NA'};
  
  comp_mode = 1:4;
end

%% Summit Camp: 2012-2014
if 0
  radartype = 'rds';
  passname = 'summit_2012_2014_allwf';
  param.multipass.fn = fullfile(gRadar.out_path,radartype,'2014_Greenland_P3','CSARP_multipass',passname);
  
  param.multipass.rbins = [];
  
  param.multipass.baseline_master_idx = 8;
  param.multipass.master_idx = 8;
  
  param.multipass.pass_en_mask = [];
  param.multipass.output_fn_midfix = [];
  param.multipass.coregistration_time_shift = [];
  param.multipass.time_gate = [];
  
%   %Load equalization vector (sar specific)
%   eqvec1 = [122.5 121.6 126.7 105.5 130.1 120.1 128.7 -0.0 -134.2 121.1 36.4 125.8 -171.5 -1.0 128.4];
%   eqvec2 = eqvec1;
%   neweq = [42.8 45.5 48.2 50.9 53.6 56.3 58.9 0.0 2.7 5.4 8.1 93.6 96.3 99.0 101.7]/2;
%   neweq = [neweq neweq];
% 
%   if 1 && exist('neweq','var') && ~isempty(neweq)
%     eqvec1 = eqvec1-neweq(1:length(eqvec1));
%     eqvec2 = eqvec2-neweq(length(eqvec1)+1:end);
%   end
% 
%   equalization1 = 10.^(zeros(1,15)/20) .* exp(1i*(eqvec1)/180*pi);
%   equalization2 = 10.^(zeros(1,15)/20) .* exp(1i*(eqvec2)/180*pi);
%   param.multipass.equalization = [equalization1 equalization2];
  
  param.multipass.debug_plots = {'debug'};
  
  comp_mode = 1:4;
end

%% Automated section
% =========================================================================
%Load param variable
[fn_dir, param.multipass.pass_name] = fileparts(param.multipass.fn);
% Input checking
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

% Run multipass
for mode_id = 1:length(comp_mode)
  param.multipass.comp_mode = comp_mode(mode_id);
  if comp_mode(mode_id)==3 %&& ~exist(fullfile(fn_dir,[param.multipass.pass_name '_slope.mat']),'file')
    %Check if slope file exists
    tomo.cross_track_slope_est(param.multipass.fn)
  end
  multipass.multipass
end