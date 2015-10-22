function loadPB_callback(obj,hObj,event)
% mapwin.loadPB_callback(obj,hObj,event)
%
% "load" pushbutton callback which loads new echo windows, called when user
% presses load button or double clicks a frame. Loads the "obj.cur_sel"
% frame.

%% Setup

% Check to make sure a frame has been selected before we load
if isempty(obj.cur_sel.frame_name)
  uiwait(msgbox('No frame selected, select frames with ctrl+left-click','Error loading','modal'));
  return;
end

% Change the pointer to a watch
set(obj.h_fig,'Pointer','watch');
drawnow;

%% Get the echowin that the current frame will be loaded into
idx = get(obj.top_panel.picker_windowPM,'Value');
if idx > 1
  %% Loading into a pre-existing window
  echo_idx = idx-1;
  exists_flag = true;
  cancel_operation = obj.echowin_list(echo_idx).undo_stack_modified_check();
  if cancel_operation
    set(obj.h_fig,'Pointer','Arrow');
    return;
  end
  obj.echowin_list(echo_idx).cmds_set_undo_stack([]);
else
  %% Loading into a new window, Create a new echowin class
  % Get the index for this new echowin
  echo_idx = length(obj.echowin_list) + 1;
  exists_flag = false;
  new_echowin = imb.echowin([],obj.default_params.echowin);
  
  % Add the new class instance to the echowin_list
  obj.echowin_list(echo_idx) = new_echowin;
  
  % Add a new entry in the picker window popup menu and set the active
  % entry to this new entry
  menu_string = get(obj.top_panel.picker_windowPM,'String');
  if strcmpi(class(obj.echowin_list(echo_idx).h_fig),'double')
    menu_string{end+1} = sprintf('%d: Echo',obj.echowin_list(echo_idx).h_fig);
  else
    menu_string{end+1} = sprintf('%d: Echo',obj.echowin_list(echo_idx).h_fig.Number);
  end
  set(obj.top_panel.picker_windowPM,'String',menu_string);
  set(obj.top_panel.picker_windowPM,'Value',echo_idx+1);
  
  % Set up the listeners
  addlistener(obj.echowin_list(echo_idx),'close_window',@obj.close_echowin);
  addlistener(obj.echowin_list(echo_idx),'update_echowin_flightline',@obj.update_echowin_flightlines);
  addlistener(obj.echowin_list(echo_idx),'update_cursors',@obj.update_echowin_cursors);
  addlistener(obj.echowin_list(echo_idx),'update_map_selection',@obj.update_map_selection_echowin);
  addlistener(obj.echowin_list(echo_idx),'open_crossover_event',@obj.open_crossover_echowin);
  
  % Create a selection plot that identifies the echowin on the map
  obj.echowin_maps(echo_idx).h_cursor = plot(obj.map_panel.h_axes, [NaN],[NaN],'kx','LineWidth',2,'MarkerSize',10);
  obj.echowin_maps(echo_idx).h_line = plot(obj.map_panel.h_axes, [NaN],[NaN],'g.');
  obj.echowin_maps(echo_idx).h_text = text(0, 0, '', 'parent', obj.map_panel.h_axes);
end

% Draw the echo class in the selected echowin
param.sources = obj.cur_map_pref_settings.sources;
param.layers = obj.cur_map_pref_settings.layers;
param.cur_sel = obj.cur_sel;
param.cur_sel.location = obj.cur_map_pref_settings.mapzone;
param.cur_sel.radar_name = obj.cur_map_pref_settings.system;  % hack for ct_filename_out to work
param.system = obj.cur_map_pref_settings.system;
try
  obj.echowin_list(echo_idx).draw(param);
catch ME
  % Draw failed... close echo window and report error
  obj.close_echowin(obj.echowin_list(echo_idx));
  set(obj.h_fig,'Pointer','Arrow');
  rethrow(ME);
end

%-------------------------------------------------------------------------
%% Create link between the echowin and undo_stack list

% Look through the unique identifiers in the undo stack document list to
% see if an undo stack already exists for this echowin's system-segment
% combination
match_idx = [];
for stack_idx = 1:length(obj.undo_stack_list)
  if strcmpi(obj.undo_stack_list(stack_idx).unique_id{1},obj.cur_map_pref_settings.system) ...
      && obj.undo_stack_list(stack_idx).unique_id{2} == obj.cur_sel.segment_id
    % An undo stack already exists for this system-segment pair
    match_idx = stack_idx;
    break;
  end
end

if isempty(match_idx)
  % An undo stack does not exist for this system-segment pair, so create a
  % new undo stack
  param.id = {obj.cur_map_pref_settings.system obj.cur_sel.segment_id};
  obj.undo_stack_list(end+1) = imb.undo_stack(param);
  match_idx = length(obj.undo_stack_list);
end

% Attach echowin to the undo stack
obj.echowin_list(echo_idx).cmds_set_undo_stack(obj.undo_stack_list(match_idx));

%% Cleanup
set(obj.h_fig,'Pointer','Arrow');

return

