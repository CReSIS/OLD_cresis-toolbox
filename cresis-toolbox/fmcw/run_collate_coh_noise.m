% script run_collate_coh_noise
%
% Runs collate_coh_noise
%
% Authors: John Paden

%% USER SETTINGS
% =========================================================================

% param_fn = ct_filename_param('snow_param_2009_Greenland_P3.xls');
% param_fn = ct_filename_param('snow_param_2009_Antarctica_DC8.xls');
% param_fn = ct_filename_param('snow_param_2010_Greenland_DC8.xls');
% param_fn = ct_filename_param('snow_param_2010_Greenland_P3.xls');
% param_fn = ct_filename_param('snow_param_2010_Antarctica_DC8.xls');
% param_fn = ct_filename_param('snow_param_2011_Greenland_P3.xls');
param_fn = ct_filename_param('snow_param_2011_Antarctica_DC8.xls');
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

coh_ave_file_input_type = 'noise';
coh_ave_file_output_type = 'noise';

analysis_sheet_name = {'analysis_coh_noise','analysis'};
day_seg = '';

debug_level = 0;

%% AUTOMATED SECTION

collate_coh_noise;

return;
