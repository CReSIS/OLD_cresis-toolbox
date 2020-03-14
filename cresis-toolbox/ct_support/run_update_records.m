% script run_records_update
%
% This script setsup the parameters and calls records_update.  Make
% a local copy of the file in your personal folder.
%
% Author: John Paden
%
% See also records_update.m

%% User Settings
% =========================================================================
param_override = [];

% Parameters spreadsheet to use for updating
%   1. Segment and frame list are taken from the parameter sheet
%   2. For GPS update, GPS time offsets are pulled from the parameter sheet
params = read_param_xls(ct_filename_param('snow_param_2018_Greenland_P3.xls'));
params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20180315');

%% Automated section
% =========================================================================

global gRadar;

% Input checking
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  records_update(param,param_override);
end
