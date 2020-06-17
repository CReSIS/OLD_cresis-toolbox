function make_layer_files(param,param_override)
% make_layer_files(param,param_override)
%
% Makes layer files for the picker.m program. This should be run
% from run_make_layer_files.m (that script sets up all the control
% variables).
%
% Author: John Paden
%
% See also: run_make_layer_files

param = merge_structs(param,param_override);

dbstack_info = dbstack;
fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', dbstack_info(1).name, param.day_seg, datestr(now,'HH:MM:SS'));
fprintf('=====================================================================\n');

%% Input Checks

if ~isfield(param.make_layer_files,'update_mode') || isempty(param.make_layer_files.update_mode)
  param.make_layer_files.update_mode = 0;
  % In all modes, if file does not exist, then a new file is created. If a
  % file exists, how it is updated is controlled by the update_mode field.
  % The update_mode field is a nonnegative integer scalar indicating:
  %
  % 0: no updates are made, layer file is created only if no layer file
  % exists or if there were problems with an existing layer file. This is
  % the default setting.
  %
  % 1: layer file is overwritten with new information (blank layer file)
  % and old layer data is lost
  %
  % 2: Old GPS time field is adjusted according to changes in
  % param.records.gps.time_offset and then the position and layer
  % information is reinterpolated from the old GPS time field to the new
  % GPS time field.
end

if ~isfield(param.make_layer_files,'out_path') || isempty(param.make_layer_files.out_path)
  param.make_layer_files.out_path = 'layer';
end

%% Make or update layer files
layers = layerdata(param,param.make_layer_files.out_path);
if param.make_layer_files.update_mode == 1
  if isempty(param.cmd.frms)
    layers.delete_all();
  else
    layers.delete_layer_file(param.cmd.frms);
  end
end
if param.make_layer_files.update_mode == 2
  layers.update_gps_all();
end
layers.check_all();
if isempty(layers.layer_organizer.lyr_name)
  layer_organizer = [];
  layer_organizer.lyr_name = {'surface','bottom'};
  layer_organizer.lyr_group_name = {'standard','standard'};
  layers.insert_layers(layer_organizer);
  layers.update_layer(1:length(layers.frames.frame_idxs),'surface',layers.frames.gps_time([1 end]),[0 0]);
  layers.update_layer(1:length(layers.frames.frame_idxs),'bottom',layers.frames.gps_time([1 end]),[NaN NaN]);
end
layers.save();
