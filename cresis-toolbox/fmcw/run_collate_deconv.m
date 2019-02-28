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
%  params = read_param_xls(ct_filename_param('snow_param_2017_Greenland_P3.xls'),'',{'analysis_spec' 'analysis'});
params = read_param_xls(ct_filename_param('snow_param_2018_Alaska_SO.xls'),'',{'analysis_spec' 'analysis'});
physical_constants;

% params = ct_set_params(params,'cmd.generic',0);
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20180529_01');

stage_one_en = true; % All deconv segments must go through stage one.

CORR_METRIC_THRESHOLD = 0.996; % Found through experimentation
CORR_METRIC_TIME_CONSTANT = 60; % Found through experimentation
TWTT_GROUPS_PER_NZ = 5; % Number of two way travel time groups per Nyquist zone
Mt = 8; % Amount to over-sample when estimating peaks of lobes

% stage_two_en: Only enable stage 2 when stage one done for all segments and make
%   sure to enable the generic column for all segments. Stage 2 is not run for rds and accum.
stage_two_en = true;

spec_file_input_type = 'analysis'; % e.g. set to 'noise' to input from CSARP_noise folder
spec_file_output_type = 'noise'; % e.g. set to 'noise' to output to CSARP_noise folder

% debug_level:
% 1. Set to 3 the very first time this is run to set Nt_shorten
% 2. Set to 4 the very first time this is run to set rbins
% 3. Set to 1 the first time for each segment to make sure good waveforms exist
% 4. Set to 0 for routine operation/re-running
debug_level = 0;

preserve_old = false; % Set to true once you have final deconv files you do not want to overwrite

imgs = 1;
wf_adcs = 1;

%% AUTOMATED SECTION
% =========================================================================

for img = imgs
  for wf_adc = wf_adcs
    collate_deconv;
  end
end

return;
