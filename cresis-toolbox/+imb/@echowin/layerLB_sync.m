function layerLB_sync(obj,type,index)
% layerLB_sync(obj,type,index)
%
% Synchronizes the listbox with the echogram when visibility/selections
% need to be changed indirectly through ctrl clicking and key presses.
% type:
%  string identifying what was changed (either 'vis' or 'sel')
% index:
%  number identifying what index was modified (just used in selections)
%

if strcmp(type,'vis')
  % just update the table
  % find the index of the top element in the listbox
  if strcmp(get(obj.left_panel.layer_panel.slider,'Enable'),'on')
    slider_max = get(obj.left_panel.layer_panel.slider,'Max');
    top_idx = slider_max+1-obj.left_panel.layer_panel.last_pointer;
  else
    top_idx = 1;
  end
  for idx=1:min(obj.left_panel.layer_panel.MAX_ROW,length(obj.left_panel.layer_panel.visible_layers))
    set(obj.left_panel.layer_panel.table.handles{idx,1},'Value',obj.left_panel.layer_panel.visible_layers(top_idx+idx-1));
  end
elseif strcmp(type,'sel')
  % find the index of the top element in the listbox
  if strcmp(get(obj.left_panel.layer_panel.slider,'Enable'),'on')
    slider_max = get(obj.left_panel.layer_panel.slider,'Max');
    top_idx = slider_max+1-obj.left_panel.layer_panel.last_pointer;
  else
    top_idx = 1;
  end
  
  % change the state variables
  val = 1; % sync always turns the layer on (never toggles)
  old_idx = find(obj.left_panel.layer_panel.selected_layers);
  if isempty(old_idx)
    old_idx = 1;
  end
  if length(index) > 1
    % if length 2+ is passed in, select both (will be for layer groups)
    lyr_idxs = (1:length(obj.left_panel.layer_panel.selected_layers)).';
    sb_mask = zeros(length(lyr_idxs),1);
    for idx = 1:length(index)
      sb_mask = sb_mask | lyr_idxs == index(idx);
    end
    obj.left_panel.layer_panel.selected_layers(sb_mask) = val;
    obj.left_panel.layer_panel.selected_layers(~sb_mask) = 0;
  else
    button_idx = index;
    
    % select this layer and unselect others
    obj.left_panel.layer_panel.selected_layers = (1:length(obj.left_panel.layer_panel.selected_layers)==button_idx & val).';
  end
  % determine if a radio button needs to be updated
  selection_visible = top_idx <= min(index) & max(index) <= top_idx+4;
  old_visible = top_idx <= min(old_idx) & max(old_idx) <= top_idx+4;
  if selection_visible || old_visible
    for idx=1:min(obj.left_panel.layer_panel.MAX_ROW,length(obj.left_panel.layer_panel.selected_layers))
      set(obj.left_panel.layer_panel.table.handles{idx,2},'Value',obj.left_panel.layer_panel.selected_layers(top_idx+idx-1));
    end
  end
end

drawnow;

return

