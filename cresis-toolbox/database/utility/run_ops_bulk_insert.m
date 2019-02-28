% CReSIS DATABASE LOADER
%
% Loads all data from CReSIS layerData into the database format
% Loads line_paths and point_paths (ops_create_path.m)
% Loads layer_points etc... (ops_create_layer_points.m)

%% USER INPUT

% SET the process command type 
%   'path' : Inserts path into database (point_paths,segments,season info,etc...)
%   'layer': Inserts layerData into database (requires path)
%   'atm': Insert the ATM layer into database (requires path)
%   'both': Executes 'path' and 'layer'
%   'all': Executes 'path', 'layer', and 'atm'
settings.process_type = 'layer';

% SET the layer name filter
%   undefined or empty --> All layers
%   Should be a function which returns logical true when passed a character
%   array containing a layer that you want inserted.
%   inline('~isempty(regexp(x,''^lm.*''))') --> Only layers starting with "lm"
%   inline('~isempty(regexp(x,''(^surface$|^bottom$)''))') --> Only surface and bottom
settings.layer_filter = inline('~isempty(regexp(x,''(^surface$|^bottom$)''))');

% SET the absolute path with filename to the param spreadsheet
settings.param_fn = ct_filename_param('rds_param_2017_Antarctica_Basler.xls');

% SET the absolute path for logs
settings.log_base_path = gRadar.tmp_path;

% SET the system name (rds, snow, accum, kuband) 
settings.sys_name = 'rds';

% SET the decimation spacing for the gps data (100m is standard for CReSIS)
settings.path_decimation_spacing = 100;

% SET the location (arctic, antarctic)
settings.location = 'antarctic';

% SET the layer post directory:
%   'layerData' for normal location
%   'qlook' for normal location using radar echograms
%   fullfile('CSARP_post','layerData') for the posting directory
%   fullfile('CSARP_post','qlook') for the qlook directory in posting
settings.layer_post_directory = 'layerData';

%% AUTOMATED SECTION

% DO SOME BASIC ERROR CHECKING
if ~any(strcmp(settings.location,{'arctic','antarctic'}));
  tmp = questdlg({'Location must be ''arctic'' or ''antarctic''','',sprintf('You entered: %s',settings.location),''},'INVALID LOCATION','OK','OK');
  error('INVALID LOCATION');
elseif ~any(strcmp(settings.sys_name,{'rds','snow','accum','kuband'}));
  tmp = questdlg({'System must be ''rds'',''snow'',''accum'' or ''kuband''','',sprintf('You entered: %s',settings.sys_name),''},'INVALID SYSTEM','OK','OK');
  error('INVALID SYSTEM');
end

global gRadar;
if isvarname('gRadar')
  settings = mergestruct(settings,gRadar);
else
  error('gRadar IS NOT A GLOBAL VARIABLE. RUN STARTUP.M');
end

ops_bulk_insert(settings);

param = [];
% Analyze tables after insertion of data. 
if strcmpi(settings.process_type, 'both')
    param.properties.tables = {'crossovers', 'point_paths', 'segments', 'locations', 'seasons',...
        'radars','frames','echograms','layer_groups','layers','layer_points'};
elseif strcmpi(settings.process_type, 'path')
    param.properties.tables = {'crossovers', 'point_paths', 'segments', 'locations',...
        'seasons','radars','frames','echograms'};
else
    param.properties.tables = {'locations','layer_groups','layers','layer_points'};
end
[status, message] = ops_analyze_tables(settings.sys_name, param);
fprintf('%s\n', message);