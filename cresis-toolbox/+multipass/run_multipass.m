% script run_multipass
%
% Script for running multipass.m
%
% Authors: Cody Barnett, Bailey Miller, John Paden
%
% See also: multipass.combine_passes.m, multipass.run_combine_passes.m,
% multipass.multipass.m, multipass.run_multipass.m

%% User Settings
% =========================================================================

param_override = [];
param = [];

example_str = 'egig_2018_allwf';

% param.multipass.layer: Choose which layers to include in the coregistered output.
% Default is surface and bottom from CSARP_layer.
%param_override.multipass.layer = struct('name',{'surface' 'bottom_qc'},'source','ops');

if strcmpi(example_str,'Thwaites_201902_201912_202001')
  %% Thwaites Line 1 20190201_01, 20191225_01, 20200127_01
  param.multipass.fn = fullfile(gRadar.out_path,'accum','2018_Antarctica_TObas','CSARP_multipass','Thwaites_201902_201912_202001.mat');
  
  param.multipass.rbins = [];
  
  param.multipass.baseline_master_idx = 1;
  param.multipass.master_idx = 1;
  
  param.multipass.pass_en_mask = [];
  param.multipass.output_fn_midfix = [];
  param.multipass.coregistration_time_shift = [0 0 0];
  param.multipass.comp_mode = 2;
end

if strcmpi(example_str,'Petermann_line1_2011_2014_2015_2017_2018_2019')
  %% Petermann Line 1 2011, 2014, 2015, 2017, 2018, 2019
  param.multipass.fn = fullfile(gRadar.out_path,'rds','2014_Greenland_P3','CSARP_multipass','Petermann_line1_2011_2014_2015_2017_2018_2019');
  
  param.multipass.rbins = [];
  
  param.multipass.baseline_master_idx = 2;
  param.multipass.master_idx = 2;
  
  param.multipass.pass_en_mask = [];
  param.multipass.output_fn_midfix = [];
  param.multipass.coregistration_time_shift = [0 0 0 0 -2 0];
  param.multipass.comp_mode = 2;
end

if strcmpi(example_str,'Petermann_line2_2013_2014_2017')
  %% Petermann Line 2 2013, 2014, 2017
  param.multipass.fn = fullfile(gRadar.out_path,'rds','2014_Greenland_P3','CSARP_multipass','Petermann_line2_2013_2014_2017');
  
  param.multipass.rbins = [];
  
  param.multipass.baseline_master_idx = 2;
  param.multipass.master_idx = 2;
  
  param.multipass.pass_en_mask = [];
  param.multipass.output_fn_midfix = [];
  param.multipass.coregistration_time_shift = [-0.5 0 0];
  param.multipass.comp_mode = 2;
end

if strcmpi(example_str,'Petermann_line4_2010_2011_2013_2014_2017')
  %% Petermann Line 4 2010, 2011, 2013, 2014, 2017
  param.multipass.fn = fullfile(gRadar.out_path,'rds','2014_Greenland_P3','CSARP_multipass','Petermann_line4_2010_2011_2013_2014_2017');
  
  param.multipass.rbins = [];
  
  param.multipass.baseline_master_idx = 2;
  param.multipass.master_idx = 2;
  
  param.multipass.pass_en_mask = [];
  param.multipass.output_fn_midfix = [];
  param.multipass.coregistration_time_shift = [0 -0.5 0 0 0];
  param.multipass.comp_mode = 2;
  
end

if strcmpi(example_str,'79N_line1_2010_2014_2016_2018_2019')
  %% 79N Line 1 2010, 2014, 2016, 2018, 2019
  param.multipass.fn = fullfile(gRadar.out_path,'rds','2014_Greenland_P3','CSARP_multipass','79N_line1_2010_2014_2016_2018_2019');
  
  param.multipass.rbins = [];
  
  param.multipass.baseline_master_idx = 2;
  param.multipass.master_idx = 2;
  
  param.multipass.pass_en_mask = [];
  param.multipass.output_fn_midfix = [];
  param.multipass.coregistration_time_shift = [-1 0 0 -2 0];
  param.multipass.comp_mode = 2;
  param.multipass.time_gate = [2e-6 13e-6];
end

if strcmpi(example_str,'Humboldt_line1_2012_2013_2014_2017')
  %% Humboldt Line 1 2012, 2013, 2014, 2017
  param.multipass.fn = fullfile(gRadar.out_path,'rds','2014_Greenland_P3','CSARP_multipass','Humboldt_line1_2012_2013_2014_2017');
  
  param.multipass.rbins = [];
  
  param.multipass.baseline_master_idx = 3;
  param.multipass.master_idx = 3;
  
  param.multipass.pass_en_mask = [];
  param.multipass.output_fn_midfix = [];
  param.multipass.coregistration_time_shift = [0 -0.5 0 0];
  param.multipass.comp_mode = 2;
  param.multipass.time_gate = [2e-6 13e-6];
end

if strcmpi(example_str,'Ryder_line1_2011_2013_2015_2019')
  %% Ryder Line 1 2011, 2013, 2015, 2019
  param.multipass.fn = fullfile(gRadar.out_path,'rds','2014_Greenland_P3','CSARP_multipass','Ryder_line1_2011_2013_2015_2019');
  
  param.multipass.rbins = [];
  
  param.multipass.baseline_master_idx = 2;
  param.multipass.master_idx = 2;
  
  param.multipass.pass_en_mask = [];
  param.multipass.output_fn_midfix = [];
  param.multipass.coregistration_time_shift = [0 0 0 0];
  param.multipass.comp_mode = 2;
  param.multipass.time_gate = [2e-6 13e-6];
end

if strcmpi(example_str,'Steensby_line1_2011_2013_2015_2019')
  %% Steensby Line 1 2011, 2013, 2015, 2019
  param.multipass.fn = fullfile(gRadar.out_path,'rds','2014_Greenland_P3','CSARP_multipass','Steensby_line1_2011_2013_2015_2019');
  
  param.multipass.rbins = [];
  
  param.multipass.baseline_master_idx = 2;
  param.multipass.master_idx = 2;
  
  param.multipass.pass_en_mask = [];
  param.multipass.output_fn_midfix = [];
  param.multipass.coregistration_time_shift = [0 0 0 0];
  param.multipass.comp_mode = 2;
  param.multipass.time_gate = [2e-6 13e-6];
end

if strcmpi(example_str,'ZI_line1_2010_2014_2016_2017_2018_2019')
  %% Zachariae Isstrom Line 1 2010, 2014, 2016, 2017, 2018, 2019
  param.multipass.fn = fullfile(gRadar.out_path,'rds','2014_Greenland_P3','CSARP_multipass','ZI_line1_2010_2014_2016_2017_2018_2019');
  
  param.multipass.rbins = [];
  
  param.multipass.baseline_master_idx = 2;
  param.multipass.master_idx = 2;
  
  param.multipass.pass_en_mask = [];
  param.multipass.output_fn_midfix = [];
  param.multipass.coregistration_time_shift = [1 0 0 0 -2 0];
  param.multipass.comp_mode = 2;
  param.multipass.time_gate = [2e-6 13e-6];
end

if strcmpi(example_str,'camp_century_2014_same_day')
  %% Camp Century: 2014 Greenland P3 same day
  param.multipass.fn = fullfile(gRadar.out_path,'rds','2014_Greenland_P3','CSARP_multipass','camp_century_2014_same_day');
  
  param.multipass.rbins = [];
  
  param.multipass.baseline_master_idx = 8;
  param.multipass.master_idx = 8;
  
  param.multipass.pass_en_mask = [];
  param.multipass.output_fn_midfix = [];
  param.multipass.coregistration_time_shift = [];
  param.multipass.comp_mode = 1:4;
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
end

if strcmpi(example_str,'camp_century_2011_2012_2013_2014')
  %% Camp Century: 2011, 2012, 2013, 2014 Greenland P3
  param.multipass.fn = fullfile(gRadar.out_path,'rds','2014_Greenland_P3','CSARP_multipass','camp_century_2011_2012_2013_2014');
  
  param.multipass.rbins = 220:420;
  
  param.multipass.baseline_master_idx = 8;
  param.multipass.master_idx = 8;
  
  param.multipass.pass_en_mask = [];
  param.multipass.output_fn_midfix = [];
  param.multipass.coregistration_time_shift = [];
  param.multipass.comp_mode = 1;
  param.multipass.time_gate = [];
end

if strcmpi(example_str,'summit_2012_2014_allwf')
  %% Summit Camp: 2012-2014
  param.multipass.fn = fullfile(gRadar.out_path,'rds','2014_Greenland_P3','CSARP_multipass','summit_2012_2014_allwf');
  
  param.multipass.rbins = [];
  
  if 1
    % 2014 only: useful for estimating cross-track slope, equalization
    param.multipass.comp_mode = 2;
    param.multipass.baseline_master_idx = 8;
    param.multipass.master_idx = 8;
    param.multipass.output_fn_midfix = '_2014';
    param.multipass.pass_en_mask = false(1,30);
    param.multipass.pass_en_mask(1:15) = true;
  elseif 0
    % All passes: useful for coregistration, equalization and differential interferometry
    param.multipass.comp_mode = 3;
    param.multipass.slope_correction_en = true;
    param.multipass.baseline_master_idx = 8;
    param.multipass.master_idx = 8;
    param.multipass.output_fn_midfix = '';
    param.multipass.pass_en_mask = false(1,30);
    param.multipass.pass_en_mask(1:30) = true;
  elseif 0
    % 2012 only using 2014 as master: useful for 2012 equalization
    param.multipass.comp_mode = 2;
    param.multipass.baseline_master_idx = 8;
    param.multipass.master_idx = 15+8;
    param.multipass.output_fn_midfix = '_2012';
    param.multipass.pass_en_mask = false(1,30);
    param.multipass.pass_en_mask(16:30) = true;
  else 0
    % 2012 only using 2012 as master: useful for 2012 equalization
    param.multipass.comp_mode = 2;
    param.multipass.baseline_master_idx = 15+8;
    param.multipass.master_idx = 15+8;
    param.multipass.output_fn_midfix = '_2012master';
    param.multipass.pass_en_mask = false(1,30);
    param.multipass.pass_en_mask(16:30) = true;
  end
  
  param.multipass.coregistration_time_shift = [];
  param.multipass.time_gate = [];
  
  param.multipass.equalization = 10.^(zeros(1,30)) ...
    .* exp(1i*([7.9 22.5 19.7 22.7 29.9 14.9 22.3 0.0 1.5 5.4 13.4 19.2 17.1 17.8 22.7 167.4 166.3 177.1 164.3 -177.6 165.9 171.6 155.3 154.6 157.5 164.6 179.0 176.7 174.8 -129.8]/180*pi));
  
  param.multipass.debug_plots = {'debug','coherent'};
end


if strcmpi(example_str,'egig_2018_allwf')
  %% EGIG line: 2011-2018
  param.multipass.fn = fullfile(gRadar.out_path,'rds','2018_Greenland_P3','CSARP_multipass','egig_2018_allwf');
  
  param.multipass.rbins = [];
  
  param.multipass.comp_mode = 2;
  param.multipass.baseline_master_idx = 8;
  param.multipass.master_idx = 8;
  param.multipass.output_fn_midfix = '_2018';
  param.multipass.pass_en_mask = true(1,30);
  
  param.multipass.coregistration_time_shift = [];
  param.multipass.time_gate = [];
  
  param.multipass.equalization = 10.^(zeros(1,30)) ...
    .* exp(1i*(zeros(1,30)/180*pi));
  
  param.multipass.debug_plots = {'debug','coherent'};
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
%multipass.multipass(param, param_override);
multipass.multipass
