% script run_update_surface_twtt_delta
%
% Runs update_surface_twtt_delta.m
%
% Author: John Paden

%% User Settings
% ----------------------------------------------------------------------
% params = read_param_xls(ct_filename_param('kuband_param_2014_Greenland_P3.xls'));
% params = read_param_xls(ct_filename_param('rds_param_2015_Greenland_Polar6.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2009_Greenland_P3.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2010_Greenland_DC8.xls'));
params = read_param_xls(ct_filename_param('snow_param_2017_Greenland_P3.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2011_Greenland_P3.xls'));
% params = read_param_xls(ct_filename_param('snow_param_2012_Greenland_P3.xls'));

params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20170510_02');
params = ct_set_params(params,'cmd.frms',[5]);
% params = ct_set_params(params,'cmd.generic',1);
% params = ct_set_params(params,'cmd.generic',0,'cmd.notes','Do not process');

% params = ct_set_params(params,'update_surface_twtt_delta.data_types',{'deconv'});
% params = ct_set_params(params,'update_surface_twtt_delta.data_types',{'CSARP_post/qlook','CSARP_post/deconv','CSARP_post/qlook_uwb','CSARP_post/qlook_kuband'});
% params = ct_set_params(params,'update_surface_twtt_delta.data_types',{'CSARP_post/qlook_uwb','CSARP_post/qlook_kuband'});
params = ct_set_params(params,'update_surface_twtt_delta.data_types',{'CSARP_post/qlook'});
params = ct_set_params(params,'update_surface_twtt_delta.imgs',[0]);
params = ct_set_params(params,'update_surface_twtt_delta.update_adc_gains_dB',1);

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
  update_surface_twtt_delta(param,param_override);
end
