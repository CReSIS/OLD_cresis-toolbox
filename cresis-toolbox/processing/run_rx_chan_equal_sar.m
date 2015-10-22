% script run_basic_rx_chan_equalization
%
% This script is for helping with setting the receiver coefficients
% from f-k processed data. It requires loading one waveform and N receive 
% channels and then use these data to compute the receiver coefficients.
%
% Author: Peng Seng Tan

fprintf('\n\n========================================================\n');
fprintf('Run receiver equalization using SAR images\n');
fprintf('========================================================\n');

%% User Settings

params = read_param_xls('/N/u/jpaden/Quarry/scripts/branch/params-cr1/rds_param_2013_Greenland_P3.xls',[],'equal');

clear('param_override');
param_override = [];
% param_override.sched.type = 'no scheduler';
param_override.sched.rerun_only = false;


%% Automated Section
% =====================================================================
% Create param structure array
% =====================================================================
tic;
global gRadar;

% Input checking
if ~exist('params','var')
  error('Use run_rx_chan_equal_sar: A struct array of parameters must be passed in\n');
end
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

% =====================================================================
% =====================================================================
% Process each of the days
% =====================================================================
% =====================================================================
for param_idx = 1:length(params)
  param = params(param_idx);
  cmd = param.cmd;
  if isfield(cmd,'generic') && cmd.generic
    rx_chan_equal_sar(param,param_override);
  end  
end

