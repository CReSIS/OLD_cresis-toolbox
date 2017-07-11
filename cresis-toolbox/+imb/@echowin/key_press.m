function key_press(obj,src,event)
% key_press(obj,src,event)
%
% Support function for mapwin class.
%
% KEY FUNCTIONS:
% -----------------------------------------------------------------------
% Space bar:
%   No mod:  cycles between all layers on or all layers off
%   Shift:   cycles between quality and regular coloring
%   Control: cycles between all layers on or only master layer on
% -----------------------------------------------------------------------
% ` (not-tilde):
%   No mod:  toggles point display and line display for all layers
%   Shift:   none
%   Control: none
% -----------------------------------------------------------------------

modifiers = get(event.Source,'CurrentModifier');
obj.shift_pressed = ismember('shift',   modifiers);  % true/false
obj.control_pressed  = ismember('control', modifiers);  % true/false
obj.alt_pressed   = ismember('alt',     modifiers);  % true/false

% Check to make sure that a key was pressed and not
% just a modifier (e.g. shift, ctrl, alt)
if ~isempty(event.Key) && ~strcmpi(event.Key,'shift') && ~strcmpi(event.Key,'alt') && ~strcmpi(event.Key,'control')
  
  % Switch layer number function (allows large numbers to be entered
  % as long as each number is entered within 0.5 seconds of the last
  % number).
  if length(event.Key) == 1 && event.Key >= '0' && event.Key <= '9'
    cur_time = now;
    done = false;
    if cur_time - obj.switch_layers.old_time < 1/86400
      % Accumulate event characters
      obj.switch_layers.accumulated_event_characters ...
        = [obj.switch_layers.accumulated_event_characters event.Key];
      desired_layer_num = str2double(obj.switch_layers.accumulated_event_characters);
      if desired_layer_num > 0 && desired_layer_num <= length(obj.left_panel.layer_panel.layer_names)
        done = true;
      end
    end
    if ~done
      % If the accumulated number does not exist OR no number was entered
      % in the last 0.5 seconds, then we start over.
      obj.switch_layers.accumulated_event_characters = event.Key;
      desired_layer_num = str2double(obj.switch_layers.accumulated_event_characters);
      if desired_layer_num > 0 && desired_layer_num <= length(obj.left_panel.layer_panel.selected_layers)
        done = true;
      end
    end
    if done
      % If the number was valid, we update the selection
      obj.layerLB_sync('sel',desired_layer_num);
      obj.set_visibility();
    end
    obj.switch_layers.old_time = cur_time;
  end
  
  
  % see event.Modifier for modifiers
  switch event.Key
    case 'f1'
      % Print out help for this window
      obj.status_text_set(sprintf('Printing help information to MATLAB console...'),'replace');
      fprintf('--------------------Echowin help--------------------\n\n');
      
      fprintf('Mouse controls:\n\n');
      fprintf('Always active regardless of mode:\n');
      fprintf('Mouse wheel scroll: Zooms in or out, centered at cursor position\n');
      fprintf('Ctrl + any click: Select closest layer and crossover\n');
      fprintf('Shift + any click: Set cursor\n');
      
      fprintf('Zoom mode active (press z to toggle zoom mode):\n');
      fprintf('Left click: Zoom in at point\n');
      fprintf('Left click and drag: select a region to zoom into\n');
      fprintf('Right click: Zoom out at point\n');
      fprintf('Double click: Full zoom out\n\n');
      
      fprintf('Tool mode active (zoom turned off):\n');
      fprintf('Left click: Appy tool to a point\n');
      fprintf('Right click and drag: Delete selected region\n');
      fprintf('Alt + Left click and drag: Apply tool to a region\n\n');
      tool_list = get(obj.left_panel.toolPM,'String');
      tool_idx = get(obj.left_panel.toolPM,'Value');
      fprintf('Currently selected tool (change in top left pulldown menu) is: %s\n',tool_list{tool_idx});
      fprintf('For this tool, the functions are...\n');
      switch tool_idx
        case 1 % interpolate
          fprintf('Left click: Enter point. Open param window (p) to set search range, enabling max point functionality (disabled by default)\n');
          fprintf('Alt + Left click and drag: Interpolate selected region. Change interpolate mode in param window (p)\n');
          fprintf('n.b. Leading edge tool (mentioned in param window) is not implemented. Leading edge and max track interpolation is not implemented\n\n');
        case 2 % quality
          fprintf('Left click: No function.\n');
          fprintf('Alt + Left click and drag: Change quality of layer points in selected region. Quality changes to the selection made in the quality pulldown menu in the left panel\n\n');
        case 3 % snake
          fprintf('Left click: Enter point. Open param window (p) to set search range, enabling max point functionality (enabled by default)\n');
          fprintf('Alt + Left click and drag: Snake selected region. Be sure ''basic'' mode is selected in the param window (p) pulldown menu\n');
          fprintf('n.b. Crandall tool (selectable in param window pulldown menu) is undocumented. Panton tool is not implemented\n\n');
        case 4 % browse
          fprintf('Left click: Open ascope window.\n');
          fprintf('Alt + Left click and drag: No function\n\n');
        case 5 % convert layers
          fprintf('Left click: No function\n');
          fprintf('Left click and drag: Convert layers within selected region to the layers specified in the param window (p)\n\n');
         case 6 % HMM detection
          fprintf('Left click: Enter point.\n');
          fprintf('Left click and drag: Perform HMM detection on selected region.\n\n');
      end
      
      fprintf('Keyboard controls:\n\n');
      fprintf('F1: This help display\n');
      fprintf('0-9: Type a sequence of numbers to switch to that layer. See layers listbox in the left panel for numbering\n');
      fprintf('a: Select both surface and bottom layers\n');
      fprintf('b: Switch to [b]rowse tool\n');
      fprintf('c: Switch to [c]onvert layers tool\n');
      fprintf('ctrl-c: Copy status bar to clipboard\n');
      fprintf('shift-c: Toggle crossover visibility (load cross overs for the first time)\n');
      fprintf('e: Switch to interpolate tool\n');
      fprintf('f: [F]lip the horizontal axis direction\n');
      fprintf('i: Switch to [i]nterpolate tool\n');
      fprintf('p: Open [p]aram window for current tool. Press p again with the window selected to close it\n');
      fprintf('q: Switch to [q]uality tool\n');
      fprintf('Shift+q: Change currently selected [q]uality (cycles between good, moderate, derived\n');
      fprintf('r: [R]edo last operation\n');
      fprintf('s: Switch to [s]nake tool\n');
      fprintf('Ctrl+s: Save [s]creenshot of current echowin (BETA)\n');
      fprintf('Shift+s: [S]ave layers\n');      
      fprintf('u: [U]ndo last operation\n');      
      fprintf('z: Enable [z]oom mode\n');
      fprintf('Ctrl+z: [Z]oom all the way out\n');
      fprintf('`: Toggle between points between using dots and dots/lines\n');
      fprintf('shift-` (~): Toggle between manual points being displayed or hidden\n');
      fprintf('Spacebar: Toggle between all layers visible and no layers visible\n');
      fprintf('Ctrl+Spacebar: Toggle between all layers visible and only currently selected layer visible\n');
      fprintf('Shift+Spacebar: Enable or disable quality display mode\n');
      fprintf('Up: Pan view field up\n');
      fprintf('Left: Pan view field left. If at left edge of frame, advance to previous frame\n');
      fprintf('Down: Pan view field down\n');
      fprintf('Right: Pan view field right. If at right edge of frame, advance to next frame\n\n\n');
      
      fprintf('Layers listbox info:\n\n');
      fprintf('The layers listbox is located in the left panel below the picker push button controls.\n');
      fprintf('The first column in the layers listbox controls visibility for each layer listed.\n');
      fprintf('The second column controls which layer is currently selected for edits.\n');
      fprintf('Normally only one layer can be selected at a time, but both surface and bottom can be selected at once, too.\n');
      fprintf('These button controls synchronize with other key and mouse controls that change layer selection and visibility.\n\n');
      
      fprintf('https://wiki.cresis.ku.edu/cresis/Data_Picking_Tutorial\n');
      
      
    case 8 % Backspace
    case 127 % Delete
      
    case 'backquote' % Toggle display of manual points
      if any(strcmp('shift',event.Modifier))
        obj.show_manual_pts = ~obj.show_manual_pts;
        obj.status_text_set(sprintf('Toggling manual point visibility'),'replace');
        obj.set_visibility();
      else
        obj.show_dots_only = ~obj.show_dots_only;
        obj.status_text_set(sprintf('Toggling dots/lines and dots'),'replace');
        obj.set_visibility();
      end
      
    case 'a'
      if ~isnan(obj.left_panel.layer_panel.surface) && ~isnan(obj.left_panel.layer_panel.bottom)
        obj.layerLB_sync('sel',[obj.left_panel.layer_panel.surface obj.left_panel.layer_panel.bottom]);
        obj.set_visibility();
      else
        obj.status_text_set(sprintf('Surface and bottom layers do not both exist'),'replace');
      end
      
    case 'b'
      set(obj.left_panel.toolPM,'Value',4);
      tmp = obj.tool_list{4}; obj.left_click = @tmp.left_click;
      tmp = obj.tool_list{4}; obj.left_click_and_drag = @tmp.left_click_and_drag;
      tmp = obj.tool_list{4}; obj.right_click_and_drag = @tmp.right_click_and_drag;
      obj.toolPM_callback();      
      
    case 'c'
      if any(strcmp('control',event.Modifier))
        obj.status_text_copy_callback();
      elseif any(strcmp('shift',event.Modifier))
        if ~obj.crossovers_en
          obj.crossovers_en = true;
          obj.load_crossovers();
        end
        % enable or disable crossovers
        obj.eg.crossovers.gui.visibility_toggle();
      else
        % convert layer tool hotkey
        set(obj.left_panel.toolPM,'Value',5);
        tmp = obj.tool_list{5}; obj.left_click = @tmp.left_click;
        tmp = obj.tool_list{5}; obj.left_click_and_drag = @tmp.left_click_and_drag;
        tmp = obj.tool_list{5}; obj.right_click_and_drag = @tmp.right_click_and_drag;
        obj.toolPM_callback();      
      end
      
    case 'e'
      set(obj.left_panel.toolPM,'Value',1);
      tmp = obj.tool_list{1}; obj.left_click = @tmp.left_click;
      tmp = obj.tool_list{1}; obj.left_click_and_drag = @tmp.left_click_and_drag;
      tmp = obj.tool_list{1}; obj.right_click_and_drag = @tmp.right_click_and_drag;
      obj.toolPM_callback();      
      
    case 'f'
      %% Flip the horizontal axis
      if strcmpi(get(obj.right_panel.axes.handle,'XDir'),'normal')
        set(obj.right_panel.axes.handle,'XDir','reverse');
      else
        set(obj.right_panel.axes.handle,'XDir','normal');
      end
      
    case 'i'
      set(obj.left_panel.toolPM,'Value',1);
      tmp = obj.tool_list{1}; obj.left_click = @tmp.left_click;
      tmp = obj.tool_list{1}; obj.left_click_and_drag = @tmp.left_click_and_drag;
      tmp = obj.tool_list{1}; obj.right_click_and_drag = @tmp.right_click_and_drag;
      obj.toolPM_callback();
      
    case 'p'
      %% Toggle param window
      obj.paramPB_callback([],[]);
      
    case 'q'
      if isempty(event.Modifier)
        set(obj.left_panel.toolPM,'Value',2);
        tmp = obj.tool_list{2}; obj.left_click = @tmp.left_click;
        tmp = obj.tool_list{2}; obj.left_click_and_drag = @tmp.left_click_and_drag;
        tmp = obj.tool_list{2}; obj.right_click_and_drag = @tmp.right_click_and_drag;
        obj.toolPM_callback();        
      elseif obj.shift_pressed && ~obj.alt_pressed && ~obj.control_pressed
        new_quality = 1+mod(get(obj.left_panel.qualityPM,'Value'), ...
          length(get(obj.left_panel.qualityPM,'String')));
        set(obj.left_panel.qualityPM,'Value',new_quality);
      end
      
    case 'r'
      %% Redo last tool operation
      obj.undo_stack.redo();
      
    case 's'
      if isempty(event.Modifier)
        set(obj.left_panel.toolPM,'Value',3);
        tmp = obj.tool_list{3}; obj.left_click = @tmp.left_click;
        tmp = obj.tool_list{3}; obj.left_click_and_drag = @tmp.left_click_and_drag;
        tmp = obj.tool_list{3}; obj.right_click_and_drag = @tmp.right_click_and_drag;
        obj.toolPM_callback();        
      elseif ~obj.shift_pressed && ~obj.alt_pressed && obj.control_pressed
        %% Save screenshot of current echowin
        if ~exist('gRadar','var')
          global gRadar;
        end        
        if ~exist(sprintf('%spicker_screenshots/',gRadar.path_override),'dir')
          mkdir(sprintf('%spicker_screenshots/',gRadar.path_override));
        end
        fn = sprintf('%spicker_screenshots/%s_%s_%s_%s_%s.png',...
          gRadar.path_override,datestr(now,'yyyymmdd_HHMMSS'),...
          obj.eg.cur_sel.radar_name,...
          obj.eg.cur_sel.season_name,obj.eg.cur_sel.day_seg,...
          num2str([obj.eg.frame_idxs(1) obj.eg.frame_idxs(end)],'%d-%d'));
        obj.status_text_set(sprintf('Saving screenshot %s...',fn),'replace');
        set(obj.h_fig,'paperpositionmode','auto');
        print('-dpng',sprintf('-f%d',obj.h_fig),fn);
      elseif obj.shift_pressed && ~obj.alt_pressed && ~obj.control_pressed
        %% Save current layers
        obj.savePB_callback();
      end
      
    case 'd'
      if isempty(event.Modifier)
        set(obj.left_panel.toolPM,'Value',6);
        tmp = obj.tool_list{3}; obj.left_click = @tmp.left_click;
        tmp = obj.tool_list{3}; obj.left_click_and_drag = @tmp.left_click_and_drag;
        tmp = obj.tool_list{3}; obj.right_click_and_drag = @tmp.right_click_and_drag;
        obj.toolPM_callback();
      end
      
    case 'u'
      %% Undo last tool operation
      obj.undo_stack.pop();
      
    case 'z'
      if obj.control_pressed
        %% Zoom reset
        obj.redraw(-inf,inf,-inf,inf,struct('clipped',true));
      else
        %% Toggle Zoom mode
        obj.status_text_set(sprintf('Toggling zoom mode'),'replace');
        obj.zoom_mode = ~obj.zoom_mode;
        if obj.zoom_mode
          set(obj.h_fig,'Pointer','custom');
        else
          set(obj.h_fig,'Pointer','Arrow');
        end
      end
      
    case 'space'
      % Space bar cycles between "currently visible layers on" and "all
      %   layers off"
      % Shift-space bar cycles between quality and regular coloring
      %   Toggles obj.tool.quality_en
      % Ctrl-Space bar cycles between "currently visible layers on" and "all
      %   layers on"
      if obj.shift_pressed
        % Toggle visibility of quality layers
        obj.tool.quality_en = ~obj.tool.quality_en;
      elseif obj.control_pressed
        % Toggle visibility between "currently visible layers on" and "all
        % layers on"
        if ~all(obj.left_panel.layer_panel.visible_layers)
          % Save this visibility state for later
          obj.tool.old_visibility = obj.left_panel.layer_panel.visible_layers;
          % Since not all layers are on, turn them all on
          obj.left_panel.layer_panel.visible_layers(:) = 1;
        else
          if isfield(obj.tool,'old_visibility') ...
              && length(obj.tool.old_visibility) == length(obj.left_panel.layer_panel.visible_layers)
            % If an old state exists, then just turn these layers on
            obj.left_panel.layer_panel.visible_layers = obj.tool.old_visibility;
          else
            % Otherwise... just turn them all off
            obj.left_panel.layer_panel.visible_layers(:) = 0;
          end
        end
        obj.layerLB_sync('vis',[]);
      else
        % Toggle visibility between "currently visible layers on" and "all
        % layers off"
        if any(obj.left_panel.layer_panel.visible_layers)
          % Save this visibility state for later
          obj.tool.old_visibility = obj.left_panel.layer_panel.visible_layers;
          % Since some layers are on, turn them all off
          obj.left_panel.layer_panel.visible_layers(:) = 0;
        else
          if isfield(obj.tool,'old_visibility') ...
              && length(obj.tool.old_visibility) == length(obj.left_panel.layer_panel.visible_layers)
            % If an old state exists, then just turn these layers on
            obj.left_panel.layer_panel.visible_layers = obj.tool.old_visibility;
          else
            % Otherwise... just turn them all back on
            obj.left_panel.layer_panel.visible_layers(:) = 1;
          end
        end
        obj.layerLB_sync('vis',[]);
      end
      obj.set_visibility();
      
    case 'downarrow' % Down-arrow: Move Echogram Down
      % check if echogram is selected
      cur_axis = [get(obj.right_panel.axes.handle,'Xlim') ...
        get(obj.right_panel.axes.handle,'YLim')];
      % Only redraw if not already at limit
      if (strcmp(obj.eg.y_order,'reverse') && cur_axis(4) < max(obj.eg.image_yaxis([1 end]))) ...
          || (strcmp(obj.eg.y_order,'normal') && cur_axis(3) > min(obj.eg.image_yaxis([1 end])))
        y_extent = cur_axis(4) - cur_axis(3);
        if strcmp(obj.eg.y_order,'reverse')
          cur_axis(3:4) = cur_axis(3:4) + y_extent*0.25;
        elseif strcmp(obj.eg.y_order,'normal')
          cur_axis(3:4) = cur_axis(3:4) - y_extent*0.25;
        end
        
        % Convert x_min, x_max to GPS time
        xlims = interp1(obj.eg.image_xaxis,obj.eg.image_gps_time,cur_axis(1:2),'linear','extrap');
        
        % Draw data with new axis, but do not allow new data to be loaded (i.e.
        % clip new axis to limits of loaded data
        obj.redraw(xlims(1),xlims(2),cur_axis(3),cur_axis(4),struct('clipped',true));
      end
      
    case 'uparrow' % Up-arrow: Move Echogram Up
      % check if echogram is selected
      cur_axis = [get(obj.right_panel.axes.handle,'Xlim') ...
        get(obj.right_panel.axes.handle,'YLim')];
      % Only redraw if not already at limit
      if (strcmp(obj.eg.y_order,'reverse') && cur_axis(3) > min(obj.eg.image_yaxis([1 end]))) ...
          || (strcmp(obj.eg.y_order,'normal') && cur_axis(4) < max(obj.eg.image_yaxis([1 end])))
        y_extent = cur_axis(4) - cur_axis(3);
        if strcmp(obj.eg.y_order,'reverse')
          cur_axis(3:4) = cur_axis(3:4) - y_extent*0.25;
        elseif strcmp(obj.eg.y_order,'normal')
          cur_axis(3:4) = cur_axis(3:4) + y_extent*0.25;
        end
        
        % Convert x_min, x_max to GPS time
        xlims = interp1(obj.eg.image_xaxis,obj.eg.image_gps_time,cur_axis(1:2),'linear','extrap');
        
        % Draw data with new axis, but do not allow new data to be loaded (i.e.
        % clip new axis to limits of loaded data
        obj.redraw(xlims(1),xlims(2),cur_axis(3),cur_axis(4),struct('clipped',true));
      end
      
    case 'rightarrow' % Right arrow
      cur_axis = [get(obj.right_panel.axes.handle,'Xlim') ...
        get(obj.right_panel.axes.handle,'YLim')];
      
      x_extent = cur_axis(2) - cur_axis(1);
      cur_axis(1:2) = cur_axis(1:2) + x_extent*0.25;
      
      % Convert x_min, x_max to GPS time
      xlims = interp1(obj.eg.image_xaxis,obj.eg.image_gps_time,cur_axis(1:2),'linear','extrap');
      
      % Draw data with new axis
      obj.redraw(xlims(1),xlims(2),cur_axis(3),cur_axis(4),struct('clipped',false));
      
    case 'leftarrow' % Left arrow
      cur_axis = [get(obj.right_panel.axes.handle,'Xlim') ...
        get(obj.right_panel.axes.handle,'YLim')];
      
      x_extent = cur_axis(2) - cur_axis(1);
      cur_axis(1:2) = cur_axis(1:2) - x_extent*0.25;
      
      % Convert x_min, x_max to GPS time
      xlims = interp1(obj.eg.image_xaxis,obj.eg.image_gps_time,cur_axis(1:2),'linear','extrap');
      
      % Draw data with new axis
      obj.redraw(xlims(1),xlims(2),cur_axis(3),cur_axis(4),struct('clipped',false));
      
  end
  obj.shift_pressed = false;
  obj.control_pressed = false;
  obj.alt_pressed = false;
end