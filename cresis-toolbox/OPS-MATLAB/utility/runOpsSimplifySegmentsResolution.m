% =========================================================================
% OPS BULK SIMPLIFY SEGMENTS RESOLUTION
%
% Changes the path resolution of segments in the database using ST_SimplifyPreserveTopology
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

clear alterParam params 
% sysName: SYSTEM NAME ('rds','snow','accum','kuband')
sysName = 'rds';

% ----------------------------------------------------------------
% resolution: the resolution in meters to which to simplify the segments using the Postgis
% command, ST_SimplifyPreserveTopology. This removes any points which can be removed while
% maintaining the original line within resolution meters.
resolution = 100;

alterParam.properties.resolution = resolution;

%% Update every segment in the db
if 1
  query = sprintf('SELECT id from %s_segments;', sysName);
  [~,data] = opsQuery(query);
  segments = [data{:}];
  alterParam.properties.segment_id = segments;
end

%% Update from param sheet
if 0

  % paramFn: FILENAME.xls OF EXCEL CReSIS PARAMS SHEET
  paramFn = 'rds_param_2015_Greenland_C130.xls';

  params = read_param_xls(ct_filename_param(paramFn));
  params = ct_set_params(params,'cmd.generic',0);
  params = ct_set_params(params,'cmd.generic',1,'day_seg','20150417_02');

  segments = {};

  for param_idx = 1:length(params)
    param = params(param_idx);
    if param.cmd.generic == 1
      if ~isempty(regexpi(param.cmd.notes,'do not process'))
        warning('You have enabled a segment with ''do not process'' in the cmd.notes, dbcont to continue');
        keyboard
      end
      segments{end + 1} = param.day_seg;
    end
  end
  alterParam.properties.segment = segments;
end

%% Call opsSimplifySegmentsResolution

%Alter the current segments
[status,message] = opsSimplifySegmentsResolution(sysName,alterParam);
if status ~= 1
  fprintf('\n');
  warning(message);
else
  fprintf(message);
  fprintf('\n');
end

clear alterParam params 
