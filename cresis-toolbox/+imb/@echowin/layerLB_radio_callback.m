function layerLB_radio_callback(obj,source,event)

% layerLB_radio_callback(obj,source,event)
%
% Layer listbox callback for radio button (layer visibility controls).
% Called when a radio button is clicked. When key shortcuts are used to change
% layer selection, other functions are used to change the selection state
% variable (see @echowin/key_press.m).
%
% Only one layer is allowed to be selected at once, except for with surface
% and bottom which can be simultaneously selected.
%

if strcmp(get(obj.left_panel.layer_panel.slider,'Enable'),'on')
  slider_max = get(obj.left_panel.layer_panel.slider,'Max');
  top_idx = slider_max+1-obj.left_panel.layer_panel.last_pointer;
else
  top_idx = 1;
end
table_idx = find(source==[obj.left_panel.layer_panel.table.handles{1:obj.left_panel.layer_panel.MAX_ROW,2}]);
val = get(source,'Value');

% check if this echogram has surface and bottom layers
surf_bot_exist_check = ~isnan(obj.left_panel.layer_panel.surface) & ~isnan(obj.left_panel.layer_panel.bottom);
button_idx = top_idx + table_idx - 1;
% if it does, check if user is selecting 'both'
if surf_bot_exist_check...
    && (button_idx==obj.left_panel.layer_panel.surface&&obj.left_panel.layer_panel.selected_layers(obj.left_panel.layer_panel.bottom)==true)...
    || (button_idx==obj.left_panel.layer_panel.bottom&&obj.left_panel.layer_panel.selected_layers(obj.left_panel.layer_panel.surface)==true)
  % this case indicates that the user has either bottom or surface selected
  % currently and is selecting the other (picking 'both')
  obj.left_panel.layer_panel.selected_layers(button_idx)=val;
else
  % otherwise, since layers can't be selected at the same time, only select
  % this layer
  obj.left_panel.layer_panel.selected_layers = (1:length(obj.left_panel.layer_panel.selected_layers)==button_idx & val).';
  % unselect other layers in the table
  if val
    for idx=1:min(obj.left_panel.layer_panel.MAX_ROW,length(obj.left_panel.layer_panel.selected_layers))
      set(obj.left_panel.layer_panel.table.handles{idx,2},'Value',obj.left_panel.layer_panel.selected_layers(top_idx+idx-1));
    end
  end
end

% draw updated colors
obj.set_visibility();

% next 3 lines are a workaround due to limitations of matlab's gui
% after a button has been pressed, the button keeps focus unless the user
% clicks elsewhere
% as a result, the next time the user presses a key shortcut, it will
% toggle the button because it is still selected
% this is notably annoying when spacebar is pressed to toggle layer
% visibility, because this callback gets called in addition to the echogram
% key press function and the result is always different from what the user
% would have anticipated
% if this code generates an error, it can be removed, but the user needs to
% click in the echogram figure after every time a button is pressed
try
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
javaFrame = get(obj.h_fig,'JavaFrame');
javaFrame.getAxisComponent.requestFocus;
catch
  obj.status_text_set(sprintf('Focus error, click inside echogram window before using key shortcuts'),'replace');
end

end

