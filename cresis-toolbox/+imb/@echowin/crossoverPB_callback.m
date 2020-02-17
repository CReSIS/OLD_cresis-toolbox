function crossoverPB_callback(obj,hObj,event)
% crossoverPB_callback(obj,hObj,event)
%
% Shows the cross over window

% Ensure focus stays on figure to prevent hotkeys registering with this
% uicontrol.
uicontrol(obj.right_panel.status_panel.statusText);

if ~obj.crossovers.en
  obj.crossovers.en = true;
  obj.load_crossovers();
  obj.crossovers.gui.visibility_toggle(true);
end

obj.crossovers.gui.figure_visibility_toggle(true);
