function picker_view_key(src,event)
% pick_view_key(src,event)
%
% Support function for picker_view (the view figure window)
%
% Author: John Paden, Aric Beaver

% Check to make sure that a key was pressed and not
% just a modifier (e.g. shift, ctrl, alt)
if ~(strcmpi(event.Key,'alt') || strcmpi(event.Key,'control') || strcmpi(event.Key,'shift'))

  global hui; % hui: user interface handles
  global gCtrl; % Geobase: contains all the geographic info

  [x,y,but] = get_mouse_info(hui.viewfig.handle,hui.viewfig.axes.handle);

  if 0 % This code for debugging
    if ischar(event.Key)
      fprintf('View: x = %.2f, y = %.2f us\n  key = %s %d\n', ...
        x, y, event.Key, event.Character);
    else
      fprintf('View: x = %.2f, y = %.2f us\n  key = %d %d\n', ...
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
    if strcmpi(event.Key,'home')
      cur_caxis = caxis;
      cur_caxis(2) = cur_caxis(2) + 1;
      caxis(cur_caxis);
    elseif strcmpi(event.Key,'end')
      cur_caxis = caxis;
      if cur_caxis(2) > cur_caxis(1) + 1
        cur_caxis(2) = cur_caxis(2) - 1;
      end
      caxis(cur_caxis);
    elseif strcmpi(event.Key,'pageup')
      cur_caxis = caxis;
      if cur_caxis(2) > cur_caxis(1) + 1
        cur_caxis(1) = cur_caxis(1) + 1;
      end
      caxis(cur_caxis);
    elseif strcmpi(event.Key,'pagedown')
      cur_caxis = caxis;
      cur_caxis(1) = cur_caxis(1) - 1;
      caxis(cur_caxis);
    end
  end
  
  % see event.Modifier for modifiers
  switch event.Character
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
      A = get(hui.viewfig.image.h,'CData');
      Time = get(hui.viewfig.image.h,'YData');
      plot(Time,A(:,gCtrl.view.cur_idx),'r');
      grid on;
      view_title = get(get(hui.viewfig.axes.handle,'Title'),'String');
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
      A = get(hui.viewfig.image.h,'CData');
      Time = get(hui.viewfig.image.h,'YData');
      plot(Time,A(:,gCtrl.view.cur_idx),'r');
      grid on;
      view_title = get(get(hui.viewfig.axes.handle,'Title'),'String');
      title(sprintf('%s\n%s', pick_title, view_title),'Interpreter','none');
      xlabel('Time (us)');
      ylabel('Relative power (dB)');
      
    case 'f'
      if strcmpi(get(hui.viewfig.axes.handle,'XDir'),'normal')
        set(hui.viewfig.axes.handle,'XDir','reverse');
      else
        set(hui.viewfig.axes.handle,'XDir','normal');        
      end

    case '!' % Multiple selection and layer 1 (surface) off/on
      if strcmp(event.Modifier{1},'shift') == 1
        gCtrl.tool.layer_visible_pick(3) = ~gCtrl.tool.layer_visible_pick(3);        
        picker_view_layer_visible;
      end
        
    case '@' % Multiple selection and layer 2 (bottom) off/on
      if strcmp(event.Modifier{1},'shift') == 1  
        gCtrl.tool.layer_visible_pick(4) = ~gCtrl.tool.layer_visible_pick(4);        
        picker_view_layer_visible;
      end
      
    case '#' % Multiple selection and (not implemented yet view only layer)
      if strcmp(event.Modifier{1},'shift') == 1  
        gCtrl.tool.layer_visible_pick(5) = ~gCtrl.tool.layer_visible_pick(5);        
        picker_view_layer_visible;
      end
      
    case ' '
      if ~isempty(event.Modifier)
        gCtrl.tool.layer_visible_view(2) = 1;
        if gCtrl.tool.layer_visible_view(1) == 1
          gCtrl.tool.layer_visible_view(1) = 2;
        else
          gCtrl.tool.layer_visible_view(1) = 1;
        end
      else
        % Layer on/off
        gCtrl.tool.layer_visible_view(2) = ~gCtrl.tool.layer_visible_view(2);
      end
      picker_view_layer_visible;
      
    case 8 % Backspace
    case 127 % Delete
    case 28 % Left-arrow
      if ~isempty(event.Modifier)
        % -----------------------------------------------------
        % Move to previous frame, but keep current axes size
        cur_axis = axis(hui.viewfig.axes.handle);
        gCtrl.source.cur_sel = gCtrl.source.cur_view;
        if gCtrl.source.cur_sel > 1
          set(hui.fig.ctrl_panel.framesLB,'Value',gCtrl.source.cur_sel - 1);
          gCtrl.source.cur_sel = get(hui.fig.ctrl_panel.framesLB,'Value');
          picker_view(1);
        end
        x_extent = cur_axis(2) - cur_axis(1);
        new_axis = cur_axis;
        new_axis(2) = length(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time);
        new_axis(1) = new_axis(2) - x_extent;
        axis(new_axis);
      else
        % -----------------------------------------------------
        % Shift axis to the left by 50%
        cur_axis = axis(hui.viewfig.axes.handle);
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
      if ~isempty(event.Modifier)
        % -----------------------------------------------------
        % Move to next frame, but keep current axes size
        cur_axis = axis(hui.viewfig.axes.handle);
      gCtrl.source.cur_sel = gCtrl.source.cur_view;
      if gCtrl.source.cur_sel < length(get(hui.fig.ctrl_panel.framesLB,'String'))
        set(hui.fig.ctrl_panel.framesLB,'Value',gCtrl.source.cur_sel + 1);
        gCtrl.source.cur_sel = get(hui.fig.ctrl_panel.framesLB,'Value');
        picker_view(1);
      end
        x_extent = cur_axis(2) - cur_axis(1);
        new_axis = cur_axis;
        new_axis(1) = 1;
        new_axis(2) = new_axis(1) + x_extent;
        axis(new_axis);
      else
        % -----------------------------------------------------
        % Shift axis to the right by 50%
        cur_axis = axis(hui.viewfig.axes.handle);
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
      if ~isempty(event.Modifier) && strcmpi(event.Modifier,'shift')
        % -----------------------------------------------------
        % Shift between source types
        if gCtrl.source.cur_src ~= 1 
          gCtrl.source.cur_src = gCtrl.source.cur_src-1;
        else
          gCtrl.source.cur_src = length(get(hui.fig.ctrl_panel.sourceLB,'String'));
        end
        set(hui.fig.ctrl_panel.sourceLB,'Value',gCtrl.source.cur_src)
        picker_view(1);
      else
        % -----------------------------------------------------
        % Shift axis up by 25%
        cur_axis = axis(hui.viewfig.axes.handle);
        y_extent = cur_axis(4) - cur_axis(3);
        new_axis = cur_axis;
        new_axis(3) = new_axis(3) - 0.25*y_extent;
        new_axis(4) = new_axis(4) - 0.25*y_extent;
        if new_axis(3) < gCtrl.view.time(1)*1e6
          new_axis(3) = gCtrl.view.time(1)*1e6;
          new_axis(4) = new_axis(3) + y_extent;
        end
        axis(new_axis);
      end
    case 31 % Down-arrow
      if ~isempty(event.Modifier) && strcmpi(event.Modifier,'shift')
        % -----------------------------------------------------
        % Shift between source types
        if gCtrl.source.cur_src ~= length(get(hui.fig.ctrl_panel.sourceLB,'String'))
          gCtrl.source.cur_src = gCtrl.source.cur_src+1;
        else
          gCtrl.source.cur_src = 1;
        end
        set(hui.fig.ctrl_panel.sourceLB,'Value',gCtrl.source.cur_src)
        picker_view(1);  
      else
        % -----------------------------------------------------
        % Shift axis down by 25%
        cur_axis = axis(hui.viewfig.axes.handle);
        y_extent = cur_axis(4) - cur_axis(3);
        new_axis = cur_axis;
        new_axis(3) = new_axis(3) + 0.25*y_extent;
        new_axis(4) = new_axis(4) + 0.25*y_extent;
        if new_axis(4) > gCtrl.view.time(end)*1e6
          new_axis(4) = gCtrl.view.time(end)*1e6;
          new_axis(3) = new_axis(4) - y_extent;
        end
        axis(new_axis);
      end
    case 'n' % Move view window to next frame
      gCtrl.source.cur_sel = gCtrl.source.cur_view;
      if gCtrl.source.cur_sel < length(get(hui.fig.ctrl_panel.framesLB,'String'))
        set(hui.fig.ctrl_panel.framesLB,'Value',gCtrl.source.cur_sel + 1);
        gCtrl.source.cur_sel = get(hui.fig.ctrl_panel.framesLB,'Value');
        picker_view(1);
      end
    case 'p' % Move view window to previous frame
      gCtrl.source.cur_sel = gCtrl.source.cur_view;
      if gCtrl.source.cur_sel > 1
        set(hui.fig.ctrl_panel.framesLB,'Value',gCtrl.source.cur_sel - 1);
        gCtrl.source.cur_sel = get(hui.fig.ctrl_panel.framesLB,'Value');
        picker_view(1);
      end
    case {'z','Z'} % zoom mode
      new_ylim = picker_ylimits(get(hui.fig.ctrl_panel.ylim_TE,'String'), ...
        gCtrl.view.time([1 end])*1e6);
      axis([1 length(gCtrl.source.geo{gCtrl.source.cur_view}.GPS_time) new_ylim]);
  end
end
