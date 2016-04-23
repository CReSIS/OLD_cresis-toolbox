% Script run_collate_deconv.m
%
% Runs collate_deconv.m
%
% Author: Jilu Li, John Paden

%% USER SETTINGS
% =========================================================================

param_fn = ct_filename_param('rds_param_2016_Greenland_Polar6.xls');

% param_fn = ct_filename_param('snow_param_2009_Greenland_P3.xls');
% param_fn = ct_filename_param('snow_param_2009_Antarctica_DC8.xls');
% param_fn = ct_filename_param('snow_param_2010_Greenland_DC8.xls');
% param_fn = ct_filename_param('snow_param_2010_Greenland_P3.xls');
% param_fn = ct_filename_param('snow_param_2010_Antarctica_DC8.xls');
% param_fn = ct_filename_param('snow_param_2011_Greenland_P3.xls');
% param_fn = ct_filename_param('snow_param_2011_Antarctica_DC8.xls');
% param_fn = ct_filename_param('snow_param_2012_Greenland_P3.xls');
% param_fn = ct_filename_param('snow_param_2012_Antarctica_DC8.xls');
% param_fn = ct_filename_param('snow_param_2014_Greenland_P3.xls');
% param_fn = ct_filename_param('snow_param_2015_Greenland_C130.xls');

% param_fn = ct_filename_param('kuband_param_2009_Greenland_P3.xls');
% param_fn = ct_filename_param('kuband_param_2009_Antarctica_DC8.xls');
% param_fn = ct_filename_param('kuband_param_2010_Greenland_DC8.xls');
% param_fn = ct_filename_param('kuband_param_2010_Greenland_P3.xls');
% param_fn = ct_filename_param('kuband_param_2010_Antarctica_DC8.xls');
% param_fn = ct_filename_param('kuband_param_2011_Greenland_P3.xls');
% param_fn = ct_filename_param('kuband_param_2011_Antarctica_DC8.xls');
% param_fn = ct_filename_param('kuband_param_2012_Greenland_P3.xls');
% param_fn = ct_filename_param('kuband_param_2012_Antarctica_DC8.xls');
% param_fn = ct_filename_param('kuband_param_2014_Greenland_P3.xls');
% param_fn = ct_filename_param('kuband_param_2015_Greenland_C130.xls');

analysis_sheet = 'analysis_spec';

physical_constants;

stage_one_en = true;
CORR_METRIC_THRESHOLD = 0.996; % Found through experimentation
CORR_METRIC_TIME_CONSTANT = 60; % Found through experimentation
TWTT_GROUPS_PER_NZ = 5; % Number of two way travel time groups per Nyquist zone
Mt = 8; % Amount to over-sample when estimating peaks of lobes

stage_two_en = true; % It is important to enable all segments in the param sheet at once for this stage

spec_file_input_type = 'noise'; % e.g. set to 'noise' to input from CSARP_noise folder
spec_file_output_type = 'noise'; % e.g. set to 'noise' to output to CSARP_noise folder

debug_level = 1; % Set to zero to run with no plots/outputs/stops

preserve_old = true; % Set to true to not overwrite old deconv file

imgs = 3;
wf_adcs = 1;

%% AUTOMATED SECTION
% =========================================================================

for img = imgs
  for wf_adc = wf_adcs
    collate_deconv;
  end
end

return;
