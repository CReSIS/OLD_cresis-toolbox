% script run_collate_coh_noise
%
% Runs collate_coh_noise
%
% Authors: John Paden

%% USER SETTINGS
% =========================================================================

% params = read_param_xls(ct_filename_param('snow_param_2009_Greenland_P3.xls'),'',{'analysis_coh_noise','analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2009_Antarctica_DC8.xls'),'',{'analysis_coh_noise','analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2010_Greenland_DC8.xls'),'',{'analysis_coh_noise','analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2010_Greenland_P3.xls'),'',{'analysis_coh_noise','analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2010_Antarctica_DC8.xls'),'',{'analysis_coh_noise','analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2011_Greenland_P3.xls'),'',{'analysis_coh_noise','analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2011_Antarctica_DC8.xls'),'',{'analysis_coh_noise','analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2012_Greenland_P3.xls'),'',{'analysis_coh_noise','analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2012_Antarctica_DC8.xls'),'',{'analysis_coh_noise','analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2014_Greenland_P3.xls'),'',{'analysis_coh_noise','analysis'});
% params = read_param_xls(ct_filename_param('snow_param_2015_Greenland_C130.xls'),'',{'analysis_coh_noise','analysis'});
params = read_param_xls(ct_filename_param('snow_param_2016_Antarctica_DC8.xls'),'20161022_05',{'analysis_coh_noise','analysis'});
params.cmd.generic = 1;

% params = read_param_xls(ct_filename_param('kuband_param_2009_Antarctica_DC8.xls'),'',{'analysis_coh_noise','analysis'});
% params = read_param_xls(ct_filename_param('kuband_param_2010_Greenland_DC8.xls'),'',{'analysis_coh_noise','analysis'});
% params = read_param_xls(ct_filename_param('kuband_param_2010_Greenland_P3.xls'),'',{'analysis_coh_noise','analysis'});
% params = read_param_xls(ct_filename_param('kuband_param_2010_Antarctica_DC8.xls'),'',{'analysis_coh_noise','analysis'});
% params = read_param_xls(ct_filename_param('kuband_param_2011_Greenland_P3.xls'),'',{'analysis_coh_noise','analysis'});
% params = read_param_xls(ct_filename_param('kuband_param_2011_Antarctica_DC8.xls'),'',{'analysis_coh_noise','analysis'});
% params = read_param_xls(ct_filename_param('kuband_param_2012_Greenland_P3.xls'),'',{'analysis_coh_noise','analysis'});
% params = read_param_xls(ct_filename_param('kuband_param_2012_Antarctica_DC8.xls'),'',{'analysis_coh_noise','analysis'});
% params = read_param_xls(ct_filename_param('kuband_param_2014_Greenland_P3.xls'),'',{'analysis_coh_noise','analysis'});
% params = read_param_xls(ct_filename_param('kuband_param_2015_Greenland_C130.xls'),'',{'analysis_coh_noise','analysis'});
params = read_param_xls(ct_filename_param('snow_param_2016_Antarctica_DC8.xls'),'20161022_05',{'analysis_coh_noise','analysis'});

coh_ave_file_input_type = 'noise';
coh_ave_file_output_type = 'noise';

debug_level = 0;
imgs = [1];
wf_adcs = [1];

%% AUTOMATED SECTION

for img = imgs
  for wf_adc = wf_adcs
    collate_coh_noise;
  end
end

return;
