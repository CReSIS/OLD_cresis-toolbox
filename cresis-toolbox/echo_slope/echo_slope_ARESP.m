function [slope, slope_corr] = echo_slope_ARESP(param, param_override)
% echo_slope_radon(param,param_override)
%
% Function uses ARESP technique to produce layer slope data
%
% param: struct with processing parameters
%         -- OR --
%         function handle to script with processing parameters
% param_override: parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Example:
%  See run_FUNCTION.m for how to run this function directly.
%  This function may be called from the run_master.m script using the
%  param spreadsheet and the cmd.generic column.
%
% Authors: Kevin Moore
%
% See also: 


%% General Setup
% =====================================================================

param = merge_structs(param, param_override);

% fprintf('=====================================================================\n');
% fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
% fprintf('=====================================================================\n');

%% Input checks
% =====================================================================
if ~isfield(param,mfilename) || isempty(param.(mfilename))
  param.(mfilename) = [];
end

if ~isfield(param.echo_slope_ARESP,'out_path') || isempty(param.echo_slope_ARESP.out_path)
  param.echo_slope_radon.out_path = 'echo_slope_radon';
end



mdata = load('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_post/CSARP_standard/20140508_01/Data_20140508_01_057.mat');



echo_slope_ARESP_task(mdata, param);

end

