function key_press(obj,src,event)
% see event.Modifier for modifiers

current_object = gco;
xlims = obj.xlims; xlims = sort(xlims([1 end]));
ylims = obj.ylims; ylims = sort(ylims([1 end]));

switch event.Key
  case 'f1'
    % Print out help for this window
    
  case 'z'
    %% toggle zoom mode
    obj.zoom_mode = ~obj.zoom_mode;
    if obj.zoom_mode
      set(obj.h_fig,'pointer','custom');
    else
      set(obj.h_fig,'pointer','arrow');
    end
    
  case 'downarrow' % Down-arrow: Pan down
    if any(current_object == obj.left_panel.ascopeLB) || any(current_object == obj.left_panel.ascopeLB)
      return
    end
    zoom_arrow(event,struct('h_axes',obj.h_axes, ...
      'xlims',xlims,'ylims',ylims));
    
  case 'uparrow' % Up-arrow: Pan up
    if any(current_object == obj.left_panel.ascopeLB) || any(current_object == obj.left_panel.ascopeLB)
      return
    end
    zoom_arrow(event,struct('h_axes',obj.h_axes, ...
      'xlims',xlims,'ylims',ylims));
    
  case 'rightarrow' % Right arrow: Pan right
    zoom_arrow(event,struct('h_axes',obj.h_axes, ...
      'xlims',xlims,'ylims',ylims));
    
  case 'leftarrow' % Left arrow: Pan left
    zoom_arrow(event,struct('h_axes',obj.h_axes, ...
      'xlims',xlims,'ylims',ylims));
end
end