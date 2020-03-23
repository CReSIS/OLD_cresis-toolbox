% function layer_update(param,param_override)
% layer_update(param,param_override)
%
% This function updates the layer according to the param struct. The
% primary purpose is to update an old layer file to the new file format.
%
% param: parameter spreadsheet structure array
%
% Examples: See run_all_layer_update.m
%
% Author: John Paden
%
% See also: run_all_layer_update, layer_update

param = merge_structs(param,param_override);

dbstack_info = dbstack;
fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', dbstack_info(1).name, param.day_seg, datestr(now,'HH:MM:SS'));
fprintf('=====================================================================\n');

%% Input checks

if ~isfield(param.layer_update,'in_path') || isempty(param.layer_update.in_path)
  param.layer_update.in_path = 'layer';
end

if ~isfield(param.layer_update,'out_path') || isempty(param.layer_update.out_path)
  param.layer_update.out_path = param.layer_update.in_path;
end

if ~isfield(param.layer_update,'frames_records_en') || isempty(param.layer_update.frames_records_en)
  param.layer_update.frames_records_en = true;
end

%% Update layers

layers = layerdata(param,param.layer_update.in_path);
layers.check_all();
layers.save(param.layer_update.out_path);
