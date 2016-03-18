% script run_collate_equalization
%
% Script for running collate_equalization which takes outputs from analysis
% surf and finds the time delay, phase, and amplitude differences between
% different wf_adc pairs.
%
% Author: John Paden
%
% See also: collate_equalization

if 1
  %% 2015_Greenland_Polar6 Setup
  
  param_fn = ct_filename_param('rds_param_2015_Greenland_Polar6.xls');
  params = read_param_xls(param_fn,'20150913_02',{'analysis_surf' 'analysis'});
  params.cmd.generic = 1;
  
  input_fn_dir = 'noise';
  
  debug_level = 1;
end

for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  fprintf('=====================================================================\n');
  fprintf('%s: %s (%s)\n', 'collate_equalization', param.day_seg, datestr(now,'HH:MM:SS'));
  fprintf('=====================================================================\n');
  
  collate_equalization;

  fprintf('  Done (%s)\n', datestr(now));
end
