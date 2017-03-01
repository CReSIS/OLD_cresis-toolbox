% =========================================================================
% Delete OPS BULK DATA
%
% Remove data in bulk from the OPS.
%
% To use this script manually edit the input fileds under USER INPUT.
%
% Data inputs required:
%   path: records and frames files
%   layerData: layerData files or echogram files
%   atm: atm L2 nadir interpolated files
%
% Multi-Layer data is supported, the layerData files must contain a name
% field with a unique name identifying each layer.
%
% Authors: Kyle W. Purdon, Trey Stafford, John Paden
%
% see also opsBulkDelete.m opsCreatePath.m opsCreaLayerPoints.m layerDataToOps.m
%
% =========================================================================

%% USER INPUT

% Initialize param structure
settings = [];

% ----------------------------------------------------------------
% runType: THE TYPE OF DATA TO DELETE
%   0: DELETE PATHS AND LAYERS
%   1: DELETE LAYERS ONLY
settings.properties.only_layer_points = 0;

% ----------------------------------------------------------------
% paramFn: FILENAME.xls OF EXCEL CReSIS PARAMS SHEET
settings.properties.season = 'rds_param_2016_Greenland_P3';
settings.properties.segment = {'20160509_01','20160509_02','20160509_05','20160509_06','20160510_01'};
% ----------------------------------------------------------------
% sys: SYSTEM NAME ('rds','snow','accum','kuband')
sys = 'rds';


%% AUTOMATED SECTION
% ERROR CHECK THE INPUT
if ~any(strcmp(sys,{'rds','snow','accum','kuband'}));
  tmp = questdlg({'System must be ''rds'',''snow'',''accum'' or ''kuband''','',sprintf('You entered: %s',settings.sysName),''},'INVALID SYSTEM','OK','OK');
  error('INVALID SYSTEM');
end

% GET THE CReSIS GLOBAL
global gRadar;
if isvarname('gRadar')
  settings = mergestruct(settings,gRadar);
else
  warning('gRadar IS NOT A GLOBAL VARIABLE.')
  fprintf('Type ''dbcont'' to run startup.m and continue\n or ''dbquit'' to fix this yourself ...\n');
  keyboard;
  startup
end

% RUN OPS BULK DELETE
[status,data] = opsDeleteBulk(sys,settings);