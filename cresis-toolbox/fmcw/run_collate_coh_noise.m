% script run_collate_coh_noise
%
% Runs collate_coh_noise
%
% Authors: John Paden

%% USER SETTINGS
% =========================================================================

param_override = [];

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
% params = read_param_xls(ct_filename_param('snow_param_2016_Antarctica_DC8.xls'),'',{'analysis_coh_noise','analysis'});

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
% params = read_param_xls(ct_filename_param('kuband_param_2016_Antarctica_DC8.xls'),'',{'analysis_coh_noise','analysis'});

params = read_param_xls(ct_filename_param('rds_param_2016_Greenland_TOdtu.xls'),'','analysis');
params = ct_set_params(params,'get_heights.coh_noise_method',17);
params = ct_set_params(params,'get_heights.coh_noise_arg',{5,251,hanning(251),'analysis'});
params = ct_set_params(params,'cmd.generic',0,'day_seg','2016110[12]');
% params = ct_set_params(params,'cmd.generic',1,'day_seg','2016111[01]_01');

param_override.collate_coh_noise.in_dir = 'analysis';
param_override.collate_coh_noise.out_dir = 'analysis';
param_override.collate_coh_noise.cmd_idx = 1;
param_override.collate_coh_noise.debug_level = 0; % Set to 1 to view threshold and minimum samples result, set to 2 to also view Doppler plots

imgs = [1];
wf_adcs = [1];

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
  for img = imgs
    for wf_adc = wf_adcs
      collate_coh_noise(param,param_override);
    end
  end
  
end