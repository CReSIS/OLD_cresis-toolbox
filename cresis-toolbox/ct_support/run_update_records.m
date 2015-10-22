% script run_update_records
%
% This script setsup the parameters and calls update_records.  Make
% a local copy of the file in your personal folder.
%
% Author: John Paden
%
% See also update_records.m

clear param;

% ======================================================================
% User Settings
% ======================================================================

% Parameters spreadsheet to use for updating
%   1. Segment and frame list are taken from the parameter sheet
%   2. For GPS update, GPS time offsets are pulled from the parameter sheet
param.param_fn = ct_filename_param('rds_param_2009_Greenland_TO.xls');

% param.skip_phrase: All segments will be skipped with this phrase in their
%   verification field. "do not process" is the standard. Leave this field
%   blank to do all segments.
param.skip_phrase = 'do not process';

% param.save_changes: Logical, For debugging purposes, you can turn the file save on/off
param.save_changes = false;

% param.debug_level: default is 1, anything higher is for debugging
param.debug_level = 1;

param.layers.update_en = false;
param.layers.source_type = 'ops'; % 'ops' or 'layerData'
param.layers.path = ''; % Required for 'layerData'
param.layers.layer_names = {'surface' 'bottom'}; % Required for 'ops'
param.layers.gaps_dist = [300 60]; % Required for 'ops'

param.gps_time.update_en = false;
param.gps_time.path = '';
param.gps_time.time_offset = -14;
param.gps_time.special_en = false;

param.gps_source.update_en = false;
param.gps_source.path = '';

update_records(param);

return;
