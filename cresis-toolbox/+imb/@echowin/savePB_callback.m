function savePB_callback(obj,hObj,event)
% savePB_callback(obj,hObj,event)

% Ensure focus stays on figure to prevent hotkeys registering with this
% uicontrol.
uicontrol(obj.right_panel.status_panel.statusText);

obj.status_text_set(sprintf('(%s) Saving...', datestr(now,'HH:MM:SS')),'replace');
drawnow;

obj.busy_mode = true;
set(obj.h_fig,'Pointer','watch');
drawnow;

imb.save_undo_stack(obj.undo_stack);
obj.eg.layers.saved.lyr_age = obj.eg.layers.lyr_age;
obj.eg.layers.saved.lyr_desc = obj.eg.layers.lyr_desc;
obj.eg.layers.saved.lyr_group_name = obj.eg.layers.lyr_group_name;
obj.eg.layers.saved.lyr_id = obj.eg.layers.lyr_id;
obj.eg.layers.saved.lyr_name = obj.eg.layers.lyr_name;
obj.eg.layers.saved.lyr_order = obj.eg.layers.lyr_order;

set(obj.h_fig,'Pointer','Arrow');
obj.busy_mode = false;

obj.status_text_set(sprintf('Done %s', datestr(now,'HH:MM:SS')),'append');
