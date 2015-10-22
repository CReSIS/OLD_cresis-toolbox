function picker_pick_key(src,event)
% picker_pick_key(src,event)
%
% Support function for picker_pick (the pick figure window)
%
% Author: John Paden, Aric Beaver

% Check to make sure that a key was pressed and not
% just a modifier (e.g. shift, ctrl, alt)
if ~(strcmpi(event.Key,'alt') || strcmpi(event.Key,'control') || strcmpi(event.Key,'shift') || strcmpi(event.Key,'return'))
  
  global hui; % hui: user interface handles
  global gCtrl; % Geobase: contains all the geographic info
  
  [x,y,but] = get_mouse_info(hui.pickfig.handle,hui.pickfig.axes.handle);
  
  if 0 % This code for debugging
    if ischar(event.Key)
      fprintf('Pick: x = %.2f, y = %.2f us\n  key = %s %d\n', ...
        x, y, event.Key, event.Character);
    else
      fprintf('Pick: x = %.2f, y = %.2f us\n  key = %d %d\n', ...
        x, y, event.Key, event.Character);
    end
    if ~isempty(event.Modifier)
      fprintf('  Modifiers ');
      for ind = 1:length(event.Modifier)
        fprintf('%s ', event.Modifier{ind});
      end
      fprintf('\n');
    end
  end
  
  if isempty(event.Character)
    increment_size = 3;
    if strcmpi(event.Key,'home')
      cur_caxis = caxis;
      cur_caxis(2) = cur_caxis(2) + increment_size;
      caxis(cur_caxis);
    elseif strcmpi(event.Key,'end')
      cur_caxis = caxis;
      if cur_caxis(2) > cur_caxis(1) + increment_size
        cur_caxis(2) = cur_caxis(2) - increment_size;
      end
      caxis(cur_caxis);
    elseif strcmpi(event.Key,'pageup')
      cur_caxis = caxis;
      if cur_caxis(2) > cur_caxis(1) + increment_size
        cur_caxis(1) = cur_caxis(1) + increment_size;
      end
      caxis(cur_caxis);
    elseif strcmpi(event.Key,'pagedown')
      cur_caxis = caxis;
      cur_caxis(1) = cur_caxis(1) - increment_size;
      caxis(cur_caxis);
    end
  end
  
  % A few of the keys act like tools and this is the setup
  % for tool operation (ctrl-up, ctrl-down, ctrl-left, ctrl-right)
  begin_command = true;
  cur_layers = get(hui.fig.ctrl_panel.layerPM,'Value');
  if cur_layers == 3
    cur_layers = [1 2];
  end
  
  
  % see event.Modifier for modifiers
  switch event.Character
    % ==========================================================
    % Tool Shortcuts
    % menuString{1} = '(b)rowse';
    % menuString{2} = '(e)nter pnt';
    % menuString{3} = '(m)ax pnt';
    % menuString{4} = '(d)elete';
    % menuString{5} = '(q)uality';
    % menuString{6} = '(i)nterp (w)';
    % menuString{7} = '(s)nake';
    % menuString{8} = '(M)ax';
    % menuString{9} = '(l)eading edge';
    % menuString{10} = 'landmar(k)';
    % menuString{11} = '(c)onvert layers';
    case 'b'
      set(hui.fig.ctrl_panel.toolPM,'Value',1);
      gCtrl.tool.selection.x = NaN;
    case 'e'
      set(hui.fig.ctrl_panel.toolPM,'Value',2);
      gCtrl.tool.selection.x = NaN;
    case 'm'
      set(hui.fig.ctrl_panel.toolPM,'Value',3);
      gCtrl.tool.selection.x = NaN;
    case 'd'
      set(hui.fig.ctrl_panel.toolPM,'Value',4);
      gCtrl.tool.selection.x = NaN;
    case 'q'
      set(hui.fig.ctrl_panel.toolPM,'Value',5);
      gCtrl.tool.selection.x = NaN;
    case 'i'
      set(hui.fig.ctrl_panel.toolPM,'Value',6);
      gCtrl.tool.selection.x = NaN;
    case 'w'
      % More user friendly, don't have to move hand from (d) or (e)
      set(hui.fig.ctrl_panel.toolPM,'Value',6);
      gCtrl.tool.selection.x = NaN;
    case 's'
      set(hui.fig.ctrl_panel.toolPM,'Value',7);
      gCtrl.tool.selection.x = NaN;
    case 'M'
      set(hui.fig.ctrl_panel.toolPM,'Value',8);
      gCtrl.tool.selection.x = NaN;
    case 'l'
      set(hui.fig.ctrl_panel.toolPM,'Value',9);
      gCtrl.tool.selection.x = NaN;
    case 'k'
      set(hui.fig.ctrl_panel.toolPM,'Value',10);
      gCtrl.tool.selection.x = NaN;
      set(hui.landfig.handle,'Visible','on');
      picker_pick_landmarks; % turn on correct landmarks (frame switching)
    case 'c'
      % Convert layer
      set(hui.fig.ctrl_panel.toolPM,'Value',11);
      gCtrl.tool.layer_switch = false;
      gCtrl.tool.selection.x = NaN;
    case 'C'
      % Replace layer
      set(hui.fig.ctrl_panel.toolPM,'Value',11);
      gCtrl.tool.layer_switch = true;
      gCtrl.tool.selection.x = NaN;
      
      
    case 'A'
      % Open A-scope of data where cursor is at
      if ishandle(5)
        figure(5);
        hold on;
        title_str = get(get(gca,'Title'),'String');
        title_str(:,end+1) = sprintf('\n');
        title_str = textscan(title_str.','%s','delimiter','\n')
        pick_title = title_str{1}{1};
        if length(title_str{1}) > 1
          view_title = title_str{1}{2};
        else
          view_title = '';
        end
      else
        figure(5); clf;
        pick_title = '';
        view_title = '';
      end
      A = get(hui.pickfig.image.h,'CData');
      Time = get(hui.pickfig.image.h,'YData');
      plot(Time,A(:,gCtrl.pick.cur_idx),'b');
      grid on;
      pick_title = get(get(hui.pickfig.axes.handle,'Title'),'String');
      title(sprintf('%s\n%s', pick_title, view_title),'Interpreter','none');
      xlabel('Time (us)');
      ylabel('Relative power (dB)');

    case 'a'
      % Open A-scope of data where cursor is at
      if ishandle(6)
        figure(6);
        hold on;
        title_str = get(get(gca,'Title'),'String');
        title_str(:,end+1) = sprintf('\n');
        title_str = textscan(title_str.','%s','delimiter','\n')
        pick_title = title_str{1}{1};
        if length(title_str{1}) > 1
          view_title = title_str{1}{2};
        else
          view_title = '';
        end
      else
        figure(6); clf;
        pick_title = '';
        view_title = '';
      end
      A = get(hui.pickfig.image.h,'CData');
      Time = get(hui.pickfig.image.h,'YData');
      plot(Time,A(:,gCtrl.pick.cur_idx),'b');
      grid on;
      pick_title = get(get(hui.pickfig.axes.handle,'Title'),'String');
      title(sprintf('%s\n%s', pick_title, view_title),'Interpreter','none');
      xlabel('Time (us)');
      ylabel('Relative power (dB)');
      
    case 'f' % flip the horizontal axis
      if strcmpi(get(hui.pickfig.axes.handle,'XDir'),'normal')
        set(hui.pickfig.axes.handle,'XDir','reverse');
      else
        set(hui.pickfig.axes.handle,'XDir','normal');        
      end
    case 'S' %Save current layers
      picker_save_layers;    
    case '1'
      set(hui.fig.ctrl_panel.layerPM,'Value',1);
      gCtrl.tool.selection.x = NaN;
    case '2'
      set(hui.fig.ctrl_panel.layerPM,'Value',2);
      gCtrl.tool.selection.x = NaN;
    case '3'
      set(hui.fig.ctrl_panel.layerPM,'Value',3);
      gCtrl.tool.selection.x = NaN;
      
    case 'u' % Undo last tool operation
      picker_undo;
      
    case '!' % Multiple selection and layer 1 (surface) off/on
      if strcmp(event.Modifier{1},'shift') == 1 && strcmp(event.Modifier{end},'control') == 1
        gCtrl.tool.layer_multiple = 1;  
      elseif strcmp(event.Modifier{1},'shift') == 1
        gCtrl.tool.layer_visible_pick(3) = ~gCtrl.tool.layer_visible_pick(3);        
        picker_pick_layer_visible;
      end
        
    case '@' % Multiple selection and layer 2 (bottom) off/on
      if strcmp(event.Modifier{1},'shift') == 1 && strcmp(event.Modifier{end},'control') == 1
        gCtrl.tool.layer_multiple = 2;
      elseif strcmp(event.Modifier{1},'shift') == 1  
        gCtrl.tool.layer_visible_pick(4) = ~gCtrl.tool.layer_visible_pick(4);        
        picker_pick_layer_visible;
      end
      
    case '#' % Multiple selection and (not implemented yet view only layer)
      if strcmp(event.Modifier{1},'shift') == 1 && strcmp(event.Modifier{end},'control') == 1
        gCtrl.tool.layer_multiple = 3;
      elseif strcmp(event.Modifier{1},'shift') == 1  
        gCtrl.tool.layer_visible_pick(5) = ~gCtrl.tool.layer_visible_pick(5);        
        picker_pick_layer_visible;
      end
      
    case '$' % Multiple selection
      gCtrl.tool.layer_multiple = 4;
      
    case ' '
      if ~isempty(event.Modifier)
        gCtrl.tool.layer_visible_pick(2) = 1;
        if gCtrl.tool.layer_visible_pick(1) == 1
          gCtrl.tool.layer_visible_pick(1) = 2
        else
          gCtrl.tool.layer_visible_pick(1) = 1
        end
      else
        % Layer on/off
        gCtrl.tool.layer_visible_pick(2) = ~gCtrl.tool.layer_visible_pick(2);
      end
      picker_pick_layer_visible;
      
    case 'Q'
      new_quality = 1+mod(get(hui.fig.ctrl_panel.qualityPM,'Value'), ...
        length(get(hui.fig.ctrl_panel.qualityPM,'String')));
      set(hui.fig.ctrl_panel.qualityPM,'Value',new_quality);
      
      % ==========================================================
      % Other Shortcuts
    case 8 % Backspace
    case 127 % Delete
    case 27 % Escape: Reset tool rubberband/box select
      gCtrl.tool.selection.x = NaN;
    case 28 % Left-arrow
      if ~isempty(event.Modifier) && strcmpi(event.Modifier,'control')
        dT = median(diff(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time));
        dX = max(1,round(0.5/dT));
        for cur_layer = cur_layers
          % Update manual layer
          layer_data = get(hui.pickfig.layer_h(2*cur_layer-1),'YData');
          layer_data(1:end-dX) = layer_data(1+dX:end);
          picker_update_layer(cur_layer,1,1:length(layer_data),layer_data,begin_command);
          begin_command = false;
          % Update automated layer
          layer_data = get(hui.pickfig.layer_h(2*cur_layer-1+1),'YData');
          layer_data(1:end-dX) = layer_data(1+dX:end);
          picker_update_layer(cur_layer,2,1:length(layer_data),layer_data,begin_command);
        end
      elseif ~isempty(event.Modifier) && strcmpi(event.Modifier,'shift')
        % -----------------------------------------------------
        % Move to previous frame, but keep current axes size
        cur_axis = axis(hui.pickfig.axes.handle);
        gCtrl.source.cur_sel = gCtrl.source.cur_pick;
        if gCtrl.source.cur_sel > 1
          set(hui.fig.ctrl_panel.framesLB,'Value',gCtrl.source.cur_sel - 1);
          gCtrl.source.cur_sel = get(hui.fig.ctrl_panel.framesLB,'Value');
          picker_pick(1);
        end
        x_extent = cur_axis(2) - cur_axis(1);
        new_axis = cur_axis;
        new_axis(2) = length(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time);
        new_axis(1) = new_axis(2) - x_extent;
        axis(new_axis);
      else
        % -----------------------------------------------------
        % Shift axis to the left by 50%
        cur_axis = axis(hui.pickfig.axes.handle);
        x_extent = cur_axis(2) - cur_axis(1);
        new_axis = cur_axis;length(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time);
        new_axis(1) = new_axis(1) - 0.5*x_extent;
        new_axis(2) = new_axis(2) - 0.5*x_extent;
        if new_axis(1) < 1
          new_axis(1) = 1;
          new_axis(2) = new_axis(1) + x_extent;
        end
        axis(new_axis);
      end
    case 29 % Right-arrow
      if ~isempty(event.Modifier) && strcmpi(event.Modifier,'control')
        dT = median(diff(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time));
        dX = max(1,round(0.5/dT));
        for cur_layer = cur_layers
          % Update manual layer
          layer_data = get(hui.pickfig.layer_h(2*cur_layer-1),'YData');
          layer_data(1+dX:end) = layer_data(1:end-dX);
          picker_update_layer(cur_layer,1,1:length(layer_data),layer_data,begin_command);
          begin_command = false;
          % Update automated layer
          layer_data = get(hui.pickfig.layer_h(2*cur_layer-1+1),'YData');
          layer_data(1+dX:end) = layer_data(1:end-dX);
          picker_update_layer(cur_layer,2,1:length(layer_data),layer_data,begin_command);
        end
      elseif ~isempty(event.Modifier) && strcmpi(event.Modifier,'shift')
        % -----------------------------------------------------
        % Move to next frame, but keep current axes size
        cur_axis = axis(hui.pickfig.axes.handle);
        gCtrl.source.cur_sel = gCtrl.source.cur_pick;
        if gCtrl.source.cur_sel < length(get(hui.fig.ctrl_panel.framesLB,'String'))
          set(hui.fig.ctrl_panel.framesLB,'Value',gCtrl.source.cur_sel + 1);
          gCtrl.source.cur_sel = get(hui.fig.ctrl_panel.framesLB,'Value');
          picker_pick(1);
        end
        x_extent = cur_axis(2) - cur_axis(1);
        new_axis = cur_axis;
        new_axis(1) = 1;
        new_axis(2) = new_axis(1) + x_extent;
        axis(new_axis);
      else
        % -----------------------------------------------------
        % Shift axis to the right by 50%
        cur_axis = axis(hui.pickfig.axes.handle);
        x_extent = cur_axis(2) - cur_axis(1);
        new_axis = cur_axis;
        new_axis(1) = new_axis(1) + 0.5*x_extent;
        new_axis(2) = new_axis(2) + 0.5*x_extent;
        if new_axis(2) > length(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time);
          new_axis(2) = length(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time);
          new_axis(1) = new_axis(2) - x_extent;
        end
        axis(new_axis);
      end
    case 30 % Up-arrow
      % -----------------------------------------------------
      % Shift layer data up
      if ~isempty(event.Modifier) && strcmpi(event.Modifier,'control')
        yAxis = get(hui.pickfig.image.h,'YData');
        dY = yAxis(2)-yAxis(1);
        for cur_layer = cur_layers
          % Update manual layer
          layer_data = get(hui.pickfig.layer_h(2*cur_layer-1),'YData');
          layer_data = layer_data - dY/2;
          picker_update_layer(cur_layer,1,1:length(layer_data),layer_data,begin_command);
          begin_command = false;
          % Update automated layer
          layer_data = get(hui.pickfig.layer_h(2*cur_layer-1+1),'YData');
          layer_data = layer_data - dY/2;
          picker_update_layer(cur_layer,2,1:length(layer_data),layer_data,begin_command);
        end
      elseif ~isempty(event.Modifier) && strcmpi(event.Modifier,'shift')
        % -----------------------------------------------------
        % Shift between source types
        if gCtrl.source.cur_src ~= 1 
          gCtrl.source.cur_src = gCtrl.source.cur_src-1;
        else
          gCtrl.source.cur_src = length(get(hui.fig.ctrl_panel.sourceLB,'String'));
        end
        set(hui.fig.ctrl_panel.sourceLB,'Value',gCtrl.source.cur_src)
        picker_pick(1);
      else
        % -----------------------------------------------------
        % Shift axis up by 25%
        cur_axis = axis(hui.pickfig.axes.handle);
        y_extent = cur_axis(4) - cur_axis(3);
        new_axis = cur_axis;
        new_axis(3) = new_axis(3) - 0.25*y_extent;
        new_axis(4) = new_axis(4) - 0.25*y_extent;
        if new_axis(3) < gCtrl.pick.time(1)*1e6
          new_axis(3) = gCtrl.pick.time(1)*1e6;
          new_axis(4) = new_axis(3) + y_extent;
        end
        axis(new_axis);
      end
    case 31 % Down-arrow
      if ~isempty(event.Modifier) && strcmpi(event.Modifier,'control')
        yAxis = get(hui.pickfig.image.h,'YData');
        dY = yAxis(2)-yAxis(1);
        for cur_layer = cur_layers
          % Update manual layer
          layer_data = get(hui.pickfig.layer_h(2*cur_layer-1),'YData');
          layer_data = layer_data + dY/2;
          picker_update_layer(cur_layer,1,1:length(layer_data),layer_data,begin_command);
          begin_command = false;
          % Update automated layer
          layer_data = get(hui.pickfig.layer_h(2*cur_layer-1+1),'YData');
          layer_data = layer_data + dY/2;
          picker_update_layer(cur_layer,2,1:length(layer_data),layer_data,begin_command);
        end
      elseif ~isempty(event.Modifier) && strcmpi(event.Modifier,'shift')
        % -----------------------------------------------------
        % Shift between source types
        if gCtrl.source.cur_src ~= length(get(hui.fig.ctrl_panel.sourceLB,'String'))
          gCtrl.source.cur_src = gCtrl.source.cur_src+1;
        else
          gCtrl.source.cur_src = 1;
        end
        set(hui.fig.ctrl_panel.sourceLB,'Value',gCtrl.source.cur_src)
        picker_pick(1);  
      else
        % -----------------------------------------------------
        % Shift axis down by 25%
        cur_axis = axis(hui.pickfig.axes.handle);
        y_extent = cur_axis(4) - cur_axis(3);
        new_axis = cur_axis;
        new_axis(3) = new_axis(3) + 0.25*y_extent;
        new_axis(4) = new_axis(4) + 0.25*y_extent;
        if new_axis(4) > gCtrl.pick.time(end)*1e6
          new_axis(4) = gCtrl.pick.time(end)*1e6;
          new_axis(3) = new_axis(4) - y_extent;
        end
        axis(new_axis);
      end
    case 'n' % Move pick window to next frame
      gCtrl.source.cur_sel = gCtrl.source.cur_pick;
      if gCtrl.source.cur_sel < length(get(hui.fig.ctrl_panel.framesLB,'String'))
        set(hui.fig.ctrl_panel.framesLB,'Value',gCtrl.source.cur_sel + 1);
        gCtrl.source.cur_sel = get(hui.fig.ctrl_panel.framesLB,'Value');
        menuString = gCtrl.source.src_disp{gCtrl.source.cur_sel};
        if gCtrl.source.cur_src > length(menuString)
          gCtrl.source.cur_src = 1;
          set(hui.fig.ctrl_panel.sourceLB,'Value',gCtrl.source.cur_src);
        end
        picker_pick(1);
      end
    case 'p' % Move pick window to previous frame
      gCtrl.source.cur_sel = gCtrl.source.cur_pick;
      if gCtrl.source.cur_sel > 1
        set(hui.fig.ctrl_panel.framesLB,'Value',gCtrl.source.cur_sel - 1);
        gCtrl.source.cur_sel = get(hui.fig.ctrl_panel.framesLB,'Value');
        menuString = gCtrl.source.src_disp{gCtrl.source.cur_sel};
        if gCtrl.source.cur_src > length(menuString)
          gCtrl.source.cur_src = 1;
          set(hui.fig.ctrl_panel.sourceLB,'Value',gCtrl.source.cur_src);
        end
        picker_pick(1);
      end
    case {'z','Z'} % zoom mode
      new_ylim = picker_ylimits(get(hui.fig.ctrl_panel.ylim_TE,'String'), ...
        gCtrl.pick.time([1 end])*1e6);
      axis([1 length(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time) new_ylim]);
  end
end

return;

