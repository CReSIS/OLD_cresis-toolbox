% Class plot_cross_overs2
%
% Class which creates a GUI for browsing ops_parse_crossovers.
%
% obj = plot_cross_overs2('C:\tmp\crossovers\crossovers.mat','C:\GIS_data\greenland\Landsat-7\mzl7geo_90m_lzw.tif');
%
% Authors: John Paden
%
% See also: ops_parse_crossovers.m

classdef (HandleCompatible = true) plot_cross_overs2 < handle
  properties
    h_fig % Map figure handle
    h_axes % Axes handle
    h_image % Image handle
    h_scatter % Scatter handle
    proj % Projection information
    h_gui % h_ctrl GUI handles struct
    h_ctrl % Control figure handle
    h_crossover % Plot handle to selected crossover
    
    crossover; % Cross over database
    season_names;
  end
  
  methods
    function obj = plot_cross_overs2(crossover_fn,geotiff_fn)
      if ~exist('geotiff_fn','var')
        geotiff_fn = '';
      end
      
      %% Load and cull cross overs
      load(crossover_fn);
      obj.crossover = data.out;
      obj.crossover.error = obj.crossover.bottom_error;
      
      good_crossover_mask = ~isnan(obj.crossover.error) ...
        & ~strcmpi('2008_Greenland_TO_wise',obj.crossover.season_1) ...
        & ~strcmpi('2008_Greenland_TO_wise',obj.crossover.season_2) ...
        & ~strcmpi('2010_Greenland_DC8',obj.crossover.season_1) ...
        & ~strcmpi('2010_Greenland_DC8',obj.crossover.season_2);
      
      obj.crossover.error = obj.crossover.error(good_crossover_mask);
      obj.crossover.X = obj.crossover.X(good_crossover_mask);
      obj.crossover.Y = obj.crossover.Y(good_crossover_mask);
      obj.crossover.frame_1 = obj.crossover.frame_1(good_crossover_mask);
      obj.crossover.frame_2 = obj.crossover.frame_2(good_crossover_mask);
      obj.crossover.gps_time_1 = obj.crossover.gps_time_1(good_crossover_mask);
      obj.crossover.gps_time_2 = obj.crossover.gps_time_2(good_crossover_mask);
      obj.crossover.season_1 = obj.crossover.season_1(good_crossover_mask);
      obj.crossover.season_2 = obj.crossover.season_2(good_crossover_mask);
      
      %% Print season statistics
      obj.season_names = unique(obj.crossover.season_1);
      for season_idx = 1:length(obj.season_names)
        season_name = obj.season_names{season_idx};
        fprintf('%d: %s\n', season_idx, season_name);
        
        season_mask = strcmpi(season_name,obj.crossover.season_1) ...
          & strcmpi(season_name,obj.crossover.season_2);
        
        fprintf('  Season Mean: %.0f m\n', mean(obj.crossover.error(season_mask)));
        fprintf('  Season Std: %.0f m\n', std(obj.crossover.error(season_mask)));
        
        season_mask = strcmpi(season_name,obj.crossover.season_1) ...
          | strcmpi(season_name,obj.crossover.season_2);
        
        fprintf('  Mean: %.0f m\n', mean(obj.crossover.error(season_mask)));
        fprintf('  Std: %.0f m\n', std(obj.crossover.error(season_mask)));
        
        [tmp,obj.crossover.idxs{season_idx}] = sort(abs(obj.crossover.error(season_mask)),'descend');
        season_mask_idxs = find(season_mask);
        obj.crossover.idxs{season_idx} = season_mask_idxs(obj.crossover.idxs{season_idx});
        for biggest_idxs_idx = 1:length(obj.crossover.idxs{season_idx})
          biggest_idx = obj.crossover.idxs{season_idx}(biggest_idxs_idx);
          if strcmpi(obj.crossover.season_1{biggest_idx}, season_name)
            other_season = obj.crossover.season_2{biggest_idx};
            other_frame = obj.crossover.frame_2{biggest_idx};
            this_frame = obj.crossover.frame_1{biggest_idx};
          else
            other_season = obj.crossover.season_1{biggest_idx};
            other_frame = obj.crossover.frame_1{biggest_idx};
            this_frame = obj.crossover.frame_2{biggest_idx};
          end
          obj.crossover.names{season_idx}{biggest_idxs_idx} ...
            = sprintf('%6.0f  %s  %s  %s', obj.crossover.error(biggest_idx), ...
            this_frame, other_frame, other_season);
        end
      end

      %% Create figure
      obj.h_fig = figure;
      set(obj.h_fig,'Numbertitle','off','Name',sprintf('%d: Crossover',obj.h_fig))
      obj.h_axes = axes('Parent',obj.h_fig);
      hold(obj.h_axes,'on');
      
      if isempty(geotiff_fn)
        obj.proj = [];
      else
        obj.proj = geotiffinfo(geotiff_fn);
        
        % Read the image
        [RGB, R, tmp] = geotiffread(geotiff_fn);
        if size(RGB,3) == 3 && strcmp(class(RGB),'uint16') && max(RGB(:)) <= 255
          RGB = uint8(RGB);
        end
        R = R/1e3;
        
        if strcmpi(class(RGB),'int16')
          RGB = double(RGB);
          RGB(RGB == 32767) = NaN;
          RGB = (RGB - min(RGB(:))) / (max(RGB(:)) - min(RGB(:)));
        end
        obj.h_image = mapshow(RGB,R,'Parent',obj.h_axes);
        xlabel('X (km)');
        ylabel('Y (km)');
      end
      
      obj.h_scatter = scatter(obj.crossover.X/1e3, ...
        obj.crossover.Y/1e3,[],obj.crossover.error, ...
        'Parent',obj.h_axes,'Marker','+');
      
      obj.h_crossover = plot(obj.h_axes,0,0,'kd','MarkerSize',7,'LineWidth',2,'MarkerFaceColor','white');
      set(obj.h_crossover,'XData',[],'YData',[]);
      
      %% Set up general handles
      set(obj.h_fig,'WindowButtonDownFcn',@obj.button_down);
      set(obj.h_fig,'WindowButtonUpFcn',@obj.button_up);
      set(obj.h_fig,'WindowScrollWheelFcn',@obj.button_scroll);
      set(obj.h_fig,'WindowKeyPressFcn',@obj.key_press);
      set(obj.h_fig,'CloseRequestFcn',@obj.close_win);
      
      %% Create control figure
      obj.h_ctrl = figure;
      set(obj.h_ctrl,'Numbertitle','off','Name',sprintf('%d: Crossover Ctrl',obj.h_ctrl))
      h_ctrl_pos = get(obj.h_ctrl,'Position');
      set(obj.h_ctrl,'Position',[h_ctrl_pos(1:2) 400 400]);
      
      season_idx = 1; % Load the first season
      
      obj.h_gui.seasonLB = uicontrol('Parent',obj.h_ctrl);
      set(obj.h_gui.seasonLB,'Style','listbox');
      set(obj.h_gui.seasonLB,'HorizontalAlignment','Center');
      set(obj.h_gui.seasonLB,'FontName','fixed');
      set(obj.h_gui.seasonLB,'String',obj.season_names);
      set(obj.h_gui.seasonLB,'Value',season_idx);
      set(obj.h_gui.seasonLB,'Callback',@obj.seasonsLB_callback);
      set(obj.h_gui.seasonLB,'Max',1);
       
      obj.h_gui.crossLB = uicontrol('Parent',obj.h_ctrl);
      set(obj.h_gui.crossLB,'Style','listbox');
      set(obj.h_gui.crossLB,'HorizontalAlignment','Center');
      set(obj.h_gui.crossLB,'FontName','fixed');
      set(obj.h_gui.crossLB,'Callback',@obj.crossLB_callback);
      set(obj.h_gui.crossLB,'Max',1);
      
      seasonsLB_callback(obj);

      %% Setup main table
      obj.h_gui.h_table.ui = obj.h_ctrl;
      obj.h_gui.h_table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
      obj.h_gui.h_table.height_margin = NaN*zeros(30,30);
      obj.h_gui.h_table.false_height = NaN*zeros(30,30);
      row = 1; col = 1;
      obj.h_gui.h_table.handles{row,col}   = obj.h_gui.seasonLB;
      obj.h_gui.h_table.width(row,col)     = inf;
      obj.h_gui.h_table.height(row,col)    = inf;
      obj.h_gui.h_table.false_height(row,col) = 0;
      row = row + 1;
      obj.h_gui.h_table.handles{row,col}   = obj.h_gui.crossLB;
      obj.h_gui.h_table.width(row,col)     = inf;
      obj.h_gui.h_table.height(row,col)    = inf;
      obj.h_gui.h_table.false_height(row,col) = 0;
      obj.h_gui.h_table.width_margin(row,col) = 0;
      obj.h_gui.h_table.height_margin(row,col) = 0;
      obj.h_gui.h_table.width_margin ...
        = obj.h_gui.h_table.width_margin(1:row,1:col);
      obj.h_gui.h_table.height_margin ...
        = obj.h_gui.h_table.height_margin(1:row,1:col);
      obj.h_gui.h_table.false_height ...
        = obj.h_gui.h_table.false_height(1:row,1:col);
      clear row col
      table_draw(obj.h_gui.h_table);
      
      %% Set up general handles
      set(obj.h_ctrl,'MenuBar','none');
      set(obj.h_ctrl,'CloseRequestFcn',@obj.close_win);
    end
    
    function delete(obj)
      % Delete the ctrl figure handle
      try
        delete(obj.h_fig);
      end
      try
        delete(obj.h_ctrl);
      end
    end
    
    function close_win(obj,h_obj,event)
      delete(obj);
    end
    
    function button_down(obj,h_obj,event)
    end
    
    function button_up(obj,h_obj,event)
      [x,y,but] = get_mouse_info(obj.h_fig,obj.h_axes);
      fprintf('%.0f %.0f\n', x, y);
      if but == 1
        % Find the closest point and print frame information
        dist = (x - obj.crossover.X/1e3).^2 + (y - obj.crossover.Y/1e3).^2;
        [min_dist min_idx] = min(dist);
        fprintf('%s:%s %s:%s %.0f\n', obj.crossover.season_1{min_idx}, ...
          obj.crossover.frame_1{min_idx}, obj.crossover.season_2{min_idx}, ...
          obj.crossover.frame_2{min_idx}, obj.crossover.error(min_idx));
        
        season_idx = get(obj.h_gui.seasonLB,'Value');
        if strcmpi(obj.crossover.season_1{min_idx},obj.season_names{season_idx}) ...
            || strcmpi(obj.crossover.season_2{min_idx},obj.season_names{season_idx})
          crossLB_idx = find(min_idx == obj.crossover.idxs{season_idx});
          set(obj.h_gui.crossLB,'Value',crossLB_idx);
          crossLB_callback(obj);
        else
          season_idx = find(strcmpi(obj.crossover.season_1{min_idx},obj.season_names));
          set(obj.h_gui.seasonLB,'Value',season_idx);
          seasonsLB_callback(obj);
          crossLB_idx = find(min_idx == obj.crossover.idxs{season_idx});
          set(obj.h_gui.crossLB,'Value',crossLB_idx);
          crossLB_callback(obj);
        end
      end
      
    end
    
    function button_scroll(obj,h_obj,event)
      [x,y,but] = get_mouse_info(obj.h_fig,obj.h_axes);
      
      % Make sure that click is on the right side panel
      mouse_pos = get(obj.h_fig,'CurrentPoint');
      
      zooms = -1 + (event.VerticalScrollCount/2);
      
      cur_axis = [get(obj.h_axes,'Xlim') ...
        get(obj.h_axes,'YLim')];
      y_extent = cur_axis(4) - cur_axis(3);
      x_extent = cur_axis(2) - cur_axis(1);
      
      % Zoom so that the mouse pointer's position in the echogram does not change
      x_percent = (x-cur_axis(1))/x_extent;
      y_percent = (y-cur_axis(3))/y_extent;
      xlims = [x - x_extent*2^(zooms+1)*x_percent, x + x_extent*2^(zooms+1)*(1-x_percent)];
      ylims = [y - y_extent*2^(zooms+1)*y_percent, y + y_extent*2^(zooms+1)*(1-y_percent)];
      
      xlim(xlims);
      ylim(ylims);
      
    end
    
    function key_press(obj,src,event)
      
      % Check to make sure that a key was pressed and not
      % just a modifier (e.g. shift, ctrl, alt)
      if ~isempty(event.Key)
        
        % see event.Modifier for modifiers
        switch event.Key
          
          case 'downarrow' % Down-arrow: Move Echogram Down
            % check if echogram is selected
            cur_axis = [get(obj.h_axes,'Xlim') ...
              get(obj.h_axes,'YLim')];
            
            y_extent = cur_axis(4) - cur_axis(3);
            if strcmp(get(obj.h_axes,'YDir'),'reverse')
              cur_axis(3:4) = cur_axis(3:4) + y_extent*0.25;
            else
              cur_axis(3:4) = cur_axis(3:4) - y_extent*0.25;
            end
            
            axis( [cur_axis(1),cur_axis(2),cur_axis(3),cur_axis(4)] );
            
          case 'uparrow' % Up-arrow: Move Echogram Up
            % check if echogram is selected
            cur_axis = [get(obj.h_axes,'Xlim') ...
              get(obj.h_axes,'YLim')];
            
            y_extent = cur_axis(4) - cur_axis(3);
            if strcmp(get(obj.h_axes,'YDir'),'reverse')
              cur_axis(3:4) = cur_axis(3:4) - y_extent*0.25;
            else
              cur_axis(3:4) = cur_axis(3:4) + y_extent*0.25;
            end
            
            axis( [cur_axis(1),cur_axis(2),cur_axis(3),cur_axis(4)] );
            
          case 'rightarrow' % Right arrow
            % check if echogram is selected
            cur_axis = [get(obj.h_axes,'Xlim') ...
              get(obj.h_axes,'YLim')];
            
            x_extent = cur_axis(2) - cur_axis(1);
            cur_axis(1:2) = cur_axis(1:2) + x_extent*0.25;
            
            axis( [cur_axis(1),cur_axis(2),cur_axis(3),cur_axis(4)] );
            
          case 'leftarrow' % Left arrow
            % check if echogram is selected
            cur_axis = [get(obj.h_axes,'Xlim') ...
              get(obj.h_axes,'YLim')];
            
            x_extent = cur_axis(2) - cur_axis(1);
            cur_axis(1:2) = cur_axis(1:2) - x_extent*0.25;
            
            axis( [cur_axis(1),cur_axis(2),cur_axis(3),cur_axis(4)] );
            
          case 'z'
            if any(strcmp('control',event.Modifier))
              axis tight;
            end
            
        end
      end
    end
    
    function seasonsLB_callback(obj,h_obj,event)
      season_idx = get(obj.h_gui.seasonLB,'Value');
      set(obj.h_gui.crossLB,'String',obj.crossover.names{season_idx});
      set(obj.h_gui.crossLB,'Value',1);
      
      crossLB_callback(obj)
    end
    
    function crossLB_callback(obj,h_obj,event)
      season_idx = get(obj.h_gui.seasonLB,'Value');
      crossLB_idx = get(obj.h_gui.crossLB,'Value');
      cross_idx = obj.crossover.idxs{season_idx}(crossLB_idx);
      x = obj.crossover.X(cross_idx)/1e3;
      y = obj.crossover.Y(cross_idx)/1e3;
      set(obj.h_crossover,'XData',x,'YData',y);
      
      % Refocus plot around current pick if the pick is not in the view
      xlims = xlim(obj.h_axes);
      ylims = ylim(obj.h_axes);
      if x < xlims(1) || x > xlims(2) || y < ylims(1) || y > ylims(2)
        xlim(obj.h_axes,xlims + x - sum(xlims)/2);
        ylim(obj.h_axes,ylims + y - sum(ylims)/2);
      end
    end
  end
end
