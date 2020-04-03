% script run_make_layer_files
%
% Calls the make_layer_files function
%
% Author: John Paden
%
% See also: make_layer_files

%% User Settings
% =========================================================================
param_override = [];

% Parameters spreadsheet to use for updating
params = read_param_xls(ct_filename_param('rds_param_2019_Antarctica_Ground.xls'));
params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20191231');
params = ct_set_params(params,'cmd.frms',[]);

% .out_path: string containing file path where layer files will be stored
% to; string is passed to ct_filename_out to form the file path.
%  Default is 'layer'.
param_override.make_layer_files.out_path = 'layer';
% param_override.make_layer_files.out_path = 'CSARP_post/layer';

% Uncomment to delete all layer files and start over:
param_override.make_layer_files.update_mode = 1;

% Uncomment to even update files that already exist:
% param_override.make_layer_files.update_mode = 2;

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
  
  make_layer_files(param,param_override);
  %make_layer_files
end
