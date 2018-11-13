% script compress_echogram
%
% Author: John Paden

%% User Settings
% ====================================================================

params = read_param_xls(ct_filename_param('snow_param_2017_Greenland_P3.xls'));

param_override.compress_echogram.echogram_dir = 'qlook';
param_override.compress_echogram.out_dir = 'CSARP_post/qlook';

% param_override.compress_echogram.echogram_dir = 'deconv';
% param_override.compress_echogram.out_dir = 'CSARP_post/deconv';
% 
% param_override.compress_echogram.echogram_dir = 'qlook_kuband';
% param_override.compress_echogram.out_dir = 'CSARP_post/qlook_kuband';
% 
% param_override.compress_echogram.echogram_dir = 'qlook_uwb';
% param_override.compress_echogram.out_dir = 'CSARP_post/qlook_uwb';

% params = ct_set_params(params,'cmd.generic',0);
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20170424_02');
% params = ct_set_params(params,'cmd.generic',1,'cmd.mission_names','^sea.*');
% params = ct_set_params(params,'cmd.generic',1,'cmd.mission_names','(?(?!^sea.*)^.*)');
params = ct_set_params(params,'cmd.generic',1);
params = ct_set_params(params,'cmd.generic',0,'cmd.notes','Do not process');
params = ct_set_params(params,'cmd.frms',[]);

%% Automated Section
% =====================================================================
% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

%% Process each of the segments
% =====================================================================
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  compress_echogram(param,param_override);
end
