
% Script for running block_data
%
% Authors: fill in...
%
% See also: run_FUNCTION.m, run_unblock.m, 

%% User Setup
% =====================================================================
param_override = [];

params = read_param_xls(ct_filename_param('rds_param_2012_Greenland_P3.xls'));

% Syntax for running a specific segment and frame by overriding parameter spreadsheet values
%params = read_param_xls(ct_filename_param('rds_param_2009_Antarctica_TO.xls'),'20091228_01');
params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20120330_03'); % 20120330_03
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20120507_05');
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20120507_07');
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20120508_04');
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20120508_06');

params = ct_set_params(params,'cmd.frms',[]);

% Set override parameters using ct_set_params, by setting param_override, or by setting FUNCTION_params variable which eventually is copied to param_override.FUNCTION).


param_override.block_data = [] ;
param_override.block_data.block_size = 256; % Set the value for each block
param_override.block_data.block_overlap = 0.5; % Set the % of overlap between each block
param_override.block_data.top_gap = 70; % Number of rows of data before first layer
% 30 for rds with no surface flattening
% 150 for snow with surface flattening


% User parameters

param_override.block_data.filter_echo_en = 0; % Enable echogram filtering before detrending
param_override.block_data.detrend_en = 1; % Enable/Disable detrending of each block


param_override.block_data.bottom_pad = 100; % Number of rows of data used to pad after the last layer

param_override.block_data.surface_flat_en = 0; % Enable/disable surface flattening. Filtering is autmaticcally disabled when flattening is off
param_override.block_data.surface_rel_layers_flat_en = 0; % Enable if layers should be relative to the flattened surface.
param_override.block_data.surface_filter_len = 1; % Typically,some surface filtering must be done if surface_flat is enabled 

param_override.block_data.pre_detrend_filter_en = 0; % Filter entire frame before "blocking"
param_override.block_data.post_detrend_filter_en = 1; % filter each block 

param_override.block_data.early_trunc = 0; % Truncate data after flattening
param_override.block_data.late_trunc = 1; % Usually 1 for rds: truncates data after normalizing


param_override.block_data.echo_path = 'CSARP_post/standard'; % snow = 'qlook',rds='CSARP_post/standard'; % Echogram source e.g rds uses 'CSARP\standard' => ct_filename_out(param,'CSARP\standard')

param_override.block_data.out_fn ='new'; % Specify desired output path i.e fn passed into ct_filename_tmp 
%param_override.sched.type = 'no scheduler'; % Example to override default cluster settings

%Paramaters of "layer_params" argument of
% OpsLoadLayers(param,"layer_params")

param_override.block_data.layers_name = 'layers';
param_override.block_data.layers_source = 'layerdata'; % snow ='layerData_koenig', 'layerData_ver1';
param_override.block_data.layerdata_source = 'layer_MacGregor';% 'layerData_koenig'; %'layerData_ver1';
param_override.block_data.regexp = 'L'; % 'L' for Macgregor and 'snow' for snow radar

% Detrend,filtering and normalization parameters
param_override.block_data.norm_detrend_params.order = 6;
param_override.block_data.norm_detrend_params.method ='polynomial';

param_override.block_data.norm_detrend_params.filter_len = 11;
param_override.block_data.norm_detrend_params.scale_min = 0.15;
param_override.block_data.norm_detrend_params.scale_max = 0.85;
param_override.block_data.norm_detrend_params.offset = 10; % Offset to be applied to the detrend curve

param_override.block_data.norm_detrend_params.roll_comp_en = 1; % 1  to apply roll compensation
param_override.block_data.norm_detrend_params.norm_window = [];



% param_override.block_data.layers_prefix = 'layers';

param_override.block_data.debug_plot = 0 ; % Debug plots (True or False)
param_override.block_data.detrend_debug = 0;


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
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  block_data(param,param_override);
end

% Post process code (only include if necessary)