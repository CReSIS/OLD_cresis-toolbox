% script run_collate_equal
%
% Script for running collate_equal.m which takes outputs from analysis
% waveform command and finds the time delay, phase, and amplitude
% differences between different wf_adc pairs.
%
% See "Receiver equalization" wiki page for details.
%
% Author: John Paden
%
% See also: collate_equal
param_override = [];

params = read_param_xls(ct_filename_param('rds_param_2011_Greenland_P3.xls'),'20110413_01',{'analysis_equal' 'analysis'});

params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20110413_01');

% param_override.collate_equal.wf_adcs = {{[2:16]},{[2:16]},{[2:16]}};
% param_override.collate_equal.rlines = [8302:12688];
% param_override.collate_equal.ref = 3;

for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  % collate_equal(param,param_override); % Execute as a function
  collate_equal; % Execute as a script (easier debugging)
end
