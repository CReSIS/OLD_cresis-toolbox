% script run_collate_fast_time
%
% Runs collate_fast_time
%  
% Bug: The frequency axis shows the fft_length for that img-wf_adc
%
% Authors: John Paden, Hara Madhav Talasila

%% USER SETTINGS
% =========================================================================

param_override = [];
param_sheet_name = 'rds_param_2019_Greenland_P3.xls';
param_fn = ct_filename_param(param_sheet_name);
params = read_param_xls(param_fn,'',{'analysis'});

param_override.analysis.enable_visible_plot = 1; % Set 1/0; Default: 0
param_override.radar.noise_fig = 2*ones(1,8); % NF linear scale, not dB
%Default noise_fig = [1.7934 1.8666 1.8038 1.6262 1.9503 1.9399 1.9608 1.9713];

% param_override.analysis.imgs = {[ones([4 1]),(1:4).'],...
% %                                [2*ones([7 1]),(1:7).']};%,...
%                                [3*ones([3 1]),(1:3).']};

% param_override.analysis.imgs = {[3*ones([1 1]),(1).']}; % For 3-1
                         
% param_override.analysis.imgs = {[2*ones([3 1]),(3:5).']}; % For 3-[3:5]

%% Automated Section
% =========================================================================

% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

% Process each of the segments
for param_idx =1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  collate_fast_time;
end
