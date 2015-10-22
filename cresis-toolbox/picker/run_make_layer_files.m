% script run_make_layer_files
%
% Calls the make_layer_files function
%
% Author: John Paden
%
% See also: make_layer_files

% ======================================================================
% User Settings
% ======================================================================

% Parameters spreadsheet to use for updating
%   1. Segment and frame list are taken from the parameter sheet
%   2. For GPS update, GPS time offsets are pulled from the parameter sheet
param.param_fn = ct_filename_param('accum_param_2015_Greenland_Ground.xls');

% param.skip_phrase: All segments will be skipped with this phrase in their
%   verification field. "do not process" is the standard. Leave this field
%   blank to do all segments.
param.skip_phrase = 'do not process';

% param.save_changes: Logical, For debugging purposes, you can turn the file save on/off
param.save_changes = true;

% Flag which prevents overwriting layer files which already exist
param.do_not_overwrite_layer_files = false;

% Flag for updating GPS values (useful when you do not want to overwrite
% the files, but you do want the GPS info to be updated). Note that
% existing layers will be re-interpolated to the new GPS time so the
% layer data does change some.
param.update_gps = true;

% Flag for modifying the GPS time so that the min/max time match between the
% old and new layer files. This is useful when timing offsets are
% present, but the data themselves don't have an offset.
param.adjust_gps_time = true;

% .combine_chan_wf_input = ct_filename_out path argument for which 
%   radar echograms to use for grabbing the initial surface values and the
%   time axis from. Typical values are shown here.
param.echogram_input = 'qlook';
% param.echogram_input = 'standard';
% param.echogram_input = 'CSARP_post/qlook';
% param.echogram_input = 'CSARP_post/standard';
param.layer_output = 'layerData';
% param.layer_output = 'CSARP_post/layerData';

make_layer_files(param);

return;




