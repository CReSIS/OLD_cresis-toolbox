% =========================================================================
% OPS BULK ALTER PATH RESOLUTION SCRIPT
%
% Changes the path resolution of segments in the database based season
% param sheets.
%
% To use this script manually edit the input fileds under USER INPUT.
%
% Data inputs required:
%   paramFn: param sheet. Generic column flags segments to alter. 
%   sysName: System ('rds','accum', 'kuband', 'snow')
%   resolution: the along-track path resolution (meters)
%

% Authors: Trey Stafford
%
% see also opsAlterPathResolution.m
%
% =========================================================================

%% USER INPUT
% ----------------------------------------------------------------
% paramFn: FILENAME.xls OF EXCEL CReSIS PARAMS SHEET
paramFn = 'rds_param_2012_Greenland_P3.xls';

% ----------------------------------------------------------------
% sysName: SYSTEM NAME ('rds','snow','accum','kuband')
sysName = 'rds';

% ----------------------------------------------------------------
% resolution: the along-track path resolution (meters). This sets the spacing of
% verticies along the flight path linestring.
resolution = 500;

%% AUTOMATED SECTION
alterParam.properties.resolution = resolution;
params = read_param_xls(ct_filename_param(paramFn));

fprintf('Altering segments ...\n');
for param_idx = 1:length(params)
  param = params(param_idx);
  if param.cmd.generic == 1
    if ~isempty(regexpi(param.cmd.notes,'do not process'))
      warning('You have enabled a segment with ''do not process'' in the cmd.notes, dbcont to continue');
      keyboard
    end
    %Get information from the param sheet. 
    alterParam.properties.segment = param.day_seg;
    alterParam.properties.season = param.season_name;
    %Alter the current segment.
    [status,message] = opsAlterPathResolution(sysName,alterParam);
    fprintf('%s\n',message)
  end
end
clear alterParam params 