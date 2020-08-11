% script run_BW_window_gen
%
% Script for running BW_window_gen
%
% Authors: John Paden, Hara Madhav Talasila
%
% See also: run_BW_window_gen.m, BW_window_gen.m

%% User Setup
% =====================================================================

param_override = [];

% params = read_param_xls(ct_filename_param('snow_param_2012_Greenland_P3.xls'));
% params = ct_set_params(params,'radar.wfs.DDC_valid',1);

% params = read_param_xls(ct_filename_param('snow_param_2013_Greenland_P3.xls'));
% params = ct_set_params(params,'radar.wfs.DDC_valid',[1 4 8 16 32]);

params = read_param_xls(ct_filename_param('snow_param_2014_Greenland_P3.xls'));
params = ct_set_params(params,'radar.wfs.DDC_valid',[1 2 4 8 16]);

% params = read_param_xls(ct_filename_param('snow_param_2019_Arctic_Polar6.xls'));
% params = ct_set_params(params,'radar.wfs(1).DDC_valid',[1 2 4 8 16]);
% params = ct_set_params(params,'radar.wfs(2).DDC_valid',[1 2 4 8 16]);

params = ct_set_params(params,'cmd.generic',1);
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20120330_03');
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20120330_04');


warning off;

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
fprintf('%s\t%s\t%s\n', 'day_seg', 'BW_window_GHz', 'BW_window');
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  if isempty(regexpi(param.cmd.notes,'do not process'))
    try
      BW_window_gen(param,param_override);
    catch
      fprintf('%s\t\n', param.day_seg);
    end
  else
    fprintf('%s\t\n', param.day_seg);
  end
end
