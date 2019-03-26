% Script run_collate_deconv.m
%
% Runs collate_deconv.m
%
% Author: Jilu Li, John Paden

%% USER SETTINGS
% =========================================================================

params = read_param_xls(ct_filename_param('snow_param_2017_Greenland_P3.xls'),'',{'analysis_spec' 'analysis'});
% params = read_param_xls(ct_filename_param('rds_param_2018_Greenland_P3.xls'),'',{'analysis_spec' 'analysis'});


% params = ct_set_params(params,'cmd.generic',1);
% params = ct_set_params(params,'cmd.generic',0,'cmd.notes','Do not process');

params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20170410_01');

params = ct_set_params(params,'collate_deconv.f0',2.85e9);
params = ct_set_params(params,'collate_deconv.f1',7.5e9);
params = ct_set_params(params,'collate_deconv.abs_metric',[58 9.8 -25 -35 inf inf]);
params = ct_set_params(params,'collate_deconv.SL_guard_bins',10);
param_override.collate_deconv.out_dir = 'analysis';


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
