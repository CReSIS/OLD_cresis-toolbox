function picker_undo
% picker_undo
%
% Function called from picker_pick_key when undo command is called
% (Ctrl-Z)
%
% Author: John Paden

global gCtrl;
global hui;

if gCtrl.undo_stack_ptr == 0
  fprintf('No commands to undo\n');
  return;
end

% Undo last change
fprintf('Undoing last change\n');
while gCtrl.undo_stack_ptr > 0
  
  % Pop the last change off of the undo stack
  undo = gCtrl.undo_stack(gCtrl.undo_stack_ptr);
  rlines = undo.rlines;
  cur_layer = undo.cur_layer;
  manual_derived = undo.manual_derived;

  % Update layer handle
  layer_data = get(hui.pickfig.layer_h(2*(cur_layer-1)+manual_derived),'YData');
  layer_data(rlines) = undo.layer_data;
  set(hui.pickfig.layer_h(2*(cur_layer-1)+manual_derived),'YData', layer_data);

  % Update quality handles
  for idx = 1:3
    quality_data = get(hui.pickfig.quality_h(3*(cur_layer-1)+idx),'YData');
    quality_data(rlines) = undo.quality_data{idx};
    set(hui.pickfig.quality_h(3*(cur_layer-1)+idx),'YData', quality_data);
  end
  
  gCtrl.undo_stack_ptr = gCtrl.undo_stack_ptr - 1;
  if gCtrl.undo_stack(gCtrl.undo_stack_ptr+1).begin_command
    % Reached the beginning of the last command so we can stop undoing
    break;
  end
end
fprintf('  Done\n');


return;
