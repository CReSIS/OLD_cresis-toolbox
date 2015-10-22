function picker_update_layer(cur_layer,manual_derived,rlines,vals,begin_command)
% picker_update_layer(cur_layer,manual_derived,rlines,vals,begin_command)
%
% Function called from picker_pick_button when a tool is used.
%
% cur_layer = current layer (1 or 2)
% manual_derived = update manual points (1) or update derived points (2)
% rlines = Nx1 vector containung which range lines to update
% vals = Nx1 vector containing values for each range line to update
% begin_command = often a single user command initiates several calls
%   to picker_update_layer, only the first of these several calls should
%   have the begin_command set to true. This is necessary so that when
%   the users chooses to undo the command, all the subcommands get undone
%   rather than just some of them
%
% Author: John Paden

global gCtrl;
global hui;

% Create a new undo stack element
undo.rlines = rlines;
undo.cur_layer = cur_layer;
undo.manual_derived = manual_derived;
undo.begin_command = begin_command;

quality = get(hui.fig.ctrl_panel.qualityPM,'Value');

% Update layer handle
layer_data = get(hui.pickfig.layer_h(2*(cur_layer-1)+manual_derived),'YData');
undo.layer_data = layer_data(rlines);
layer_data(rlines) = vals;
set(hui.pickfig.layer_h(2*(cur_layer-1)+manual_derived),'YData', layer_data);

% Update quality handles
for idx = 1:3
  quality_data = get(hui.pickfig.quality_h(3*(cur_layer-1)+idx),'YData');
  undo.quality_data{idx} = quality_data(rlines);
  if idx == quality
    quality_data(rlines) = vals;
  else
    quality_data(rlines) = inf;
  end
  set(hui.pickfig.quality_h(3*(cur_layer-1)+idx),'YData', quality_data);
end

% Add modification to undo stack
if gCtrl.undo_stack_ptr == 0
  gCtrl.undo_stack_ptr = gCtrl.undo_stack_ptr + 1;
  gCtrl.undo_stack = undo;
elseif gCtrl.undo_stack_ptr < 100
  gCtrl.undo_stack_ptr = gCtrl.undo_stack_ptr + 1;
  gCtrl.undo_stack(gCtrl.undo_stack_ptr:end) = [];
  gCtrl.undo_stack(gCtrl.undo_stack_ptr) = undo;
else
  gCtrl.undo_stack(1:99) = gCtrl.undo_stack(2:100);
  gCtrl.undo_stack(100) = undo;
end

return;
