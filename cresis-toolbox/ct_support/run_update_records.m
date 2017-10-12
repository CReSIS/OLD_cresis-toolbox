% script run_update_records
%
% This script setsup the parameters and calls update_records.  Make
% a local copy of the file in your personal folder.
%
% Author: John Paden
%
% See also update_records.m

% ======================================================================
% User Settings
% ======================================================================

% Parameters spreadsheet to use for updating
%   1. Segment and frame list are taken from the parameter sheet
%   2. For GPS update, GPS time offsets are pulled from the parameter sheet
params = read_param_xls(ct_filename_param('kuband_param_2016_Antarctica_DC8.xls'),'20161112_03','post');
params.cmd.generic = 1;
% params = read_param_xls(ct_filename_param('rds_param_2009_Greenland_TO.xls'),'','post');

update_records(params);
