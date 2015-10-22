% Class column_browser
%
% Class which creates a GUI for browsing the columns of a matrix that has
% been loaded into an image handle.
%
% obj = column_browser(img_title,h_img,cur_column,hide_only)

classdef (HandleCompatible = true) column_browser < handle
  properties
    h_fig
    h_gui
    h_img
    h_line
    cur_column
    hide_only
    get_data_fh
    normalize_en
    zoom_mode
    xdir_reverse
    x
    y
    xlims
    ylims
  end
  
  methods
    function obj = column_browser(img_title,h_img,cur_column,hide_only,get_data_fh)
      if ~exist('get_data_fh','var')
        get_data_fh = [];
      end
      
      obj.xdir_reverse = true;
      obj.normalize_en = false;
      obj.zoom_mode = true;
      obj.h_img = h_img;
      obj.cur_column = cur_column;
      obj.hide_only = hide_only;
      obj.get_data_fh = get_data_fh;
      
      obj.create_ui(img_title);
      
    end
    
    function delete(obj)
      % Delete the map figure handle
      try
        delete(obj.h_fig);
      end
    end
    
    function close_win(obj,h_obj,event)
      delete(obj);
    end
    
    function button_down(obj,h_obj,event)
      [obj.x,obj.y,but] = get_mouse_info(obj.h_fig,obj.h_gui.h_axes);
      rbbox;
    end
    
    function button_up(obj,h_obj,event)
      [x,y,but] = get_mouse_info(obj.h_fig,obj.h_gui.h_axes);
      
      if obj.zoom_mode
        zoom_button_up(x,y,but,struct('x',obj.x,'y',obj.y,'h_axes',obj.h_gui.h_axes,'xlims',obj.xlims,'ylims',obj.ylims));
      end
      
    end
    
    function button_scroll(obj,h_obj,event)
      zoom_button_scroll(event,struct('h_fig',obj.h_fig,'h_axes',obj.h_gui.h_axes));
    end
    
    function key_press(obj,src,event)
      
      if any(strcmp('shift',event.Modifier))
        shift_pressed = true;
      else
        shift_pressed = false;
      end
      
      if any(strcmp('ctrl',event.Modifier))
        ctrl_pressed = true;
      else
        ctrl_pressed = false;
      end
      
      % Check to make sure that a key was pressed and not
      % just a modifier (e.g. shift, ctrl, alt)
      if ~isempty(event.Key)
        
        % see event.Modifier for modifiers
        switch event.Key
          
          case 'z'
            if ctrl_pressed
              %% zoom reset
              axis tight;
            else
              %% toggle zoom mode
              obj.zoom_mode = ~obj.zoom_mode;
              if obj.zoom_mode
                set(obj.h_fig,'pointer','custom');
              else
                set(obj.h_fig,'pointer','arrow');
              end
            end
          
          case 'n'
            obj.cur_column = obj.cur_column + 1;
            obj.update_plot();
            
          case 'p'
            obj.cur_column = obj.cur_column - 1;
            obj.update_plot();

          case 'downarrow' % Down-arrow: Move Echogram Down
            zoom_arrow(event,struct('h_axes',obj.h_gui.h_axes));
            
          case 'uparrow' % Up-arrow: Move Echogram Up
            zoom_arrow(event,struct('h_axes',obj.h_gui.h_axes));
            
          case 'rightarrow' % Right arrow
            zoom_arrow(event,struct('h_axes',obj.h_gui.h_axes));
            
          case 'leftarrow' % Left arrow
            zoom_arrow(event,struct('h_axes',obj.h_gui.h_axes));
            
        end
        
      end
    end
    
    function update_plot(obj)
      obj.cur_column
      if length(obj.h_gui.h_plot) ~= length(obj.h_img)
        delete(obj.h_gui.h_plot);
        update_plots = true;
      else
        update_plots = false;
      end
      colors = {'k','r','g','c','b'};
      obj.xlims = [inf -inf];
      obj.ylims = [inf -inf];
      for idx = 1:length(obj.h_img)
        colors_idx = mod(idx-1,length(colors))+1;
        if ~isempty(obj.get_data_fh)
          [XData,YData] = obj.get_data_fh(obj.h_img(idx),cur_column);
        else
          YData = get(obj.h_img(idx),'CData');
          YData = YData(:,obj.cur_column);
          XData = get(obj.h_img(idx),'YData');
        end
        if min(XData(isfinite(XData))) < obj.xlims(1)
          obj.xlims(1) = min(XData(isfinite(XData)));
        end
        if max(XData(isfinite(XData))) > obj.xlims(2)
          obj.xlims(2) = max(XData(isfinite(XData)));
        end
        if min(YData(isfinite(YData))) < obj.ylims(1)
          obj.ylims(1) = min(YData(isfinite(YData)));
        end
        if max(YData(isfinite(YData))) > obj.ylims(2)
          obj.ylims(2) = max(YData(isfinite(YData)));
        end
        if update_plots
          obj.h_gui.h_plot(idx) = plot(XData, YData, colors{colors_idx}, 'Parent', obj.h_gui.h_axes);
        else
          set(obj.h_gui.h_plot(idx),{'XData','YData'},{XData,YData});
        end
      end

      if obj.xdir_reverse
        set(obj.h_gui.h_axes,'XDir','reverse');
      end

      if obj.normalize_en
        obj.normalize_plot();
      else
        obj.add_markers();
        obj.ylims = obj.ylims + [-0.05 0.05] * diff(obj.ylims);
      end
      
    end
    
    function normalize_plot(obj)
      max_val = -inf;
      for idx = 1:length(obj.h_img)
        new_max_val = max(get(obj.h_gui.h_plot(idx),'YData'));
        if new_max_val > max_val
          max_val = new_max_val;
        end
      end
      obj.ylims = [inf -inf];
      for idx = 1:length(obj.h_img)
        YData = get(obj.h_gui.h_plot(idx),'YData');
        YData = YData - max(YData) + new_max_val;
        set(obj.h_gui.h_plot(idx),'YData',YData);
        if min(YData(isfinite(YData))) < obj.ylims(1)
          obj.ylims(1) = min(YData(isfinite(YData)));
        end
        if max(YData(isfinite(YData))) > obj.ylims(2)
          obj.ylims(2) = max(YData(isfinite(YData)));
        end
      end
      obj.ylims = obj.ylims + [-0.05 0.05] * diff(obj.ylims);
      
      obj.add_markers();
    end
        
    function add_markers(obj)
      idx = 1;
      XData = get(obj.h_gui.h_plot(idx),'XData');
      YData = get(obj.h_gui.h_plot(idx),'YData');
      [max_val,max_idx] = max(YData);
      
      max_idx = max_idx - 6 + find(YData(max_idx-5:max_idx) > max_val-3,1);
      peak_offset = find(diff(YData(max_idx+(0:20))) < 0,1);
      if isempty(peak_offset)
        peak_offset = 20;
      end
      max_idx = max_idx + peak_offset - 1;
      
      mean_val = nanmean(YData(1:50));
      thresh_idx = find(YData > mean_val+5,1);
      peak_offset = find(diff(YData(thresh_idx+(0:20))) < 0,1);
      if isempty(peak_offset)
        peak_offset = 20;
      end
      thresh_idx = thresh_idx + peak_offset - 1;
      
      for idx = 1:length(obj.h_gui.h_marker)
        try
          delete(obj.h_gui.h_marker);
        end
      end
      obj.h_gui.h_marker(1) = plot(XData([thresh_idx thresh_idx]),obj.ylims,'b');
      obj.h_gui.h_marker(2) = plot(XData([max_idx max_idx]),obj.ylims,'b');
    end

    create_ui(obj,img_title);
    
  end
end




