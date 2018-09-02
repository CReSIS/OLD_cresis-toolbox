% Script run_collate_deconv.m
%
% Runs collate_deconv.m
%
% Author: Jilu Li, John Paden

%% USER SETTINGS
% =========================================================================

% params = read_param_xls(ct_filename_param('kuband_param_2009_Greenland_P3.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('kuband_param_2009_Antarctica_DC8.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('kuband_param_2010_Greenland_DC8.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('kuband_param_2010_Greenland_P3.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('kuband_param_2010_Antarctica_DC8.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('kuband_param_2011_Greenland_P3.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('kuband_param_2011_Antarctica_DC8.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('kuband_param_2012_Greenland_P3.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('kuband_param_2012_Antarctica_DC8.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('kuband_param_2014_Greenland_P3.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('kuband_param_2015_Greenland_C130.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('kuband_param_2016_Greenland_P3.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('kuband_param_2016_Antarctica_DC8.xls'),'',{'analysis_spec' 'analysis'});

% params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('rds_param_2016_Greenland_Polar6.xls'),'',{'analysis_spec' 'analysis'});

% params = read_param_xls(ct_filename_param('snow_param_2009_Greenland_P3.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2009_Antarctica_DC8.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2010_Greenland_DC8.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2010_Greenland_P3.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2010_Antarctica_DC8.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2011_Greenland_P3.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2011_Antarctica_DC8.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2012_Greenland_P3.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2012_Antarctica_DC8.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2014_Antarctica_DC8.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2014_Greenland_P3.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2015_Greenland_C130.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2016_Greenland_P3.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2016_Antarctica_DC8.xls'),'',{'analysis_spec' 'analysis'});
params = read_param_xls(ct_filename_param('snow_param_2017_Greenland_P3.xls'),'',{'analysis_spec' 'analysis'});
% params = ct_set_params(params,'cmd.generic',0);

% params = ct_set_params(params,'cmd.generic',1,'day_seg','20170314_01'); % Good

% params = ct_set_params(params,'cmd.generic',1,'day_seg','20170320');

% params = ct_set_params(params,'cmd.generic',1,'day_seg','20170310_02');
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20170314_02');
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20170324_02');

% 2-8 GHz: 20170410_01, 20170323_0[23]

% 2018: [-40 5 -30 -40 inf inf]

% 2-18 GHz Deconvolution Settings (3 sets)
params = ct_set_params(params,'analysis.cmd{1}.f0',2.85e9);
params = ct_set_params(params,'analysis.cmd{1}.f1',14e9);
params = ct_set_params(params,'analysis.cmd{1}.abs_metric',[58 4.5 -25 -35 inf inf]);
params = ct_set_params(params,'analysis.cmd{1}.SL_guard_bins',6);
param_override.collate_deconv.out_dir = 'analysis_uwb';
% params = ct_set_params(params,'analysis.cmd{1}.day_segs',{'20170324_01'});

% params = ct_set_params(params,'analysis.cmd{1}.f0',2.85e9);
% params = ct_set_params(params,'analysis.cmd{1}.f1',7.5e9);
% params = ct_set_params(params,'analysis.cmd{1}.abs_metric',[58 9.2 -25 -35 inf inf]);
% params = ct_set_params(params,'analysis.cmd{1}.SL_guard_bins',10);
% param_override.collate_deconv.out_dir = 'analysis';

% params = ct_set_params(params,'analysis.cmd{1}.f0',12e9);
% params = ct_set_params(params,'analysis.cmd{1}.f1',14e9);
% params = ct_set_params(params,'analysis.cmd{1}.abs_metric',[58 23 -25 -28 inf inf]);
% params = ct_set_params(params,'analysis.cmd{1}.SL_guard_bins',20);
% param_override.collate_deconv.out_dir = 'analysis_kuband';

% 2-8 GHz Deconvolution Settings
% params = ct_set_params(params,'analysis.cmd{1}.f0',2.85e9);
% params = ct_set_params(params,'analysis.cmd{1}.f1',7.5e9);
% params = ct_set_params(params,'analysis.cmd{1}.abs_metric',[65 4.5 -25 -35 inf inf]);
% params = ct_set_params(params,'analysis.cmd{1}.SL_guard_bins',6);
% param_override.collate_deconv.out_dir = 'analysis';
% params = ct_set_params(params,'analysis.cmd{1}.day_segs',{'20170324_01'});

param_override.collate_deconv.gps_time_penalty = 1/(8*3600);

param_override.collate_deconv.cmd_idx = 1;
% debug_level:
% 1. Set to 3 the very first time this is run to set Nt_shorten
% 2. Set to 4 the very first time this is run to set rbins
% 3. Set to 1 the first time for each segment to make sure good waveforms exist
% 4. Set to 0 for routine operation/re-running
param_override.collate_deconv.debug_level = 1;
param_override.collate_deconv.imgs = 1;
param_override.collate_deconv.wf_adcs = [];
param_override.collate_deconv.stage_one_en = false; % All deconv segments must go through stage one.
param_override.collate_deconv.stage_two_en = true; % Stage 2 only required for fmcw and not rds or accum
param.collate_deconv.Mt = 10; % Oversampling amount for peak measurements (Mt=10 recommended)

param_override.collate_deconv.preserve_old = false;
param_override.collate_deconv.file_lock = false;

%% Automated Section
% =====================================================================

% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

% Process each of the segments
ctrl_chain = {};
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  %collate_deconv(param,param_override);
  collate_deconv
  
end
