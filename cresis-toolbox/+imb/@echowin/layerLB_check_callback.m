function layerLB_check_callback(obj,source,event)

% layerLB_check_callback(obj,source,event)
%
% Layer listbox callback for check boxes (layer visibility controls).
% Called when a checkbox is clicked. When key shortcuts are used to change
% layer visibility, other functions are used to change the visibility state
% variable (see @echowin/key_press.m)
%


if strcmp(get(obj.left_panel.layer_panel.slider,'Enable'),'on')
  slider_max = get(obj.left_panel.layer_panel.slider,'Max');
  top_idx = slider_max+1-obj.left_panel.layer_panel.last_pointer;
else
  top_idx = 1;
end
table_idx = find(source==[obj.left_panel.layer_panel.table.handles{1:obj.left_panel.layer_panel.MAX_ROW,1}]);
val = get(source,'Value');

box_idx = top_idx + table_idx - 1;

obj.left_panel.layer_panel.visible_layers(box_idx) = val;

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

