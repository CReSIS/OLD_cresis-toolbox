classdef ice_mask_edit < handle
  
  properties
    
    %JORDAN
    slice
    ice_mask_fn
    mdata
    proj
    
    h_flight_dem_plot
    h_flight_mask_plot
    h_true_dem_plot
    h_false_dem_plot
    h_true_mask_plot
    h_false_mask_plot
    
    flight_line
    mdata_loaded
    intersections
    %JORDAN
    
    dem
    gray
    mask
    
    intensity_thresh
    intensity_thresh_default
    
    h_dem_fig
    h_dem_axes
    h_dem_plot
    
    h_mask_fig
    h_mask_axes
    h_mask_plot
    
    x
    y
    zoom_mode
    
    dem_x_mesh
    dem_y_mesh
    dem_x
    dem_y
    ice_x_mesh
    ice_y_mesh
    ice_x
    ice_y
    
    actions
    tool_list
    active_tool_idx
    
    gui
    local_undo_stack
    local_undo_flag
    
    cmd
    
    shift_pressed
    ctrl_pressed
    
  end
  
  events
    IceChange
    SliceChange
    Undo
    Redo
  end
  
  methods
%     function obj = ice_mask_edit(dem_x_mesh,dem_y_mesh,dem,ice_x_mesh,ice_y_mesh,ice_mask)
      function obj = ice_mask_edit(param)
      
      fprintf('In ice_mask_edit\n');  
      
      if isfield(param,'DEM')
        obj.dem = param.DEM;
      end
      
      if isfield(param,'R')
        obj.dem_x = param.R(3,1) + param.R(2,1)*(0:size(param.DEM,2)-1);
        obj.dem_y = param.R(3,2) + param.R(1,2)*(0:size(param.DEM,1)-1);
      end
      
      if isfield(param,'ice_mask')
        obj.mask = param.ice_mask.mask;
        obj.ice_x = param.ice_mask.X;
        obj.ice_y = param.ice_mask.Y;
      end
            
      if isfield(param,'ice_mask_fn');
        obj.ice_mask_fn = param.ice_mask_fn;
      end
      
      if isfield(param,'mdata')
%         obj.mdata.twtt = param.mdata.twtt;
%         obj.mdata.theta = param.mdata.theta;
%         obj.mdata.ice_mask = param.mdata.ice_mask;
%         obj.mdata.param_combine = param.mdata.param_combine;
%         rmfield(param,'mdata');
        obj.mdata = param.mdata;
      end
      
      if isfield(param,'proj')
        obj.proj = param.proj;
      end
      
      undo_param.id = [];
      obj.local_undo_stack = imb.undo_stack(undo_param);
      addlistener(obj.local_undo_stack,'synchronize_event',@obj.undo_sync);
      obj.local_undo_flag = 1;
      
      obj.dem = obj.dem(obj.dem_y>=min(obj.ice_y) & obj.dem_y<=max(obj.ice_y) , ...
        obj.dem_x>=min(obj.ice_x) & obj.dem_x<=max(obj.ice_x),:);
      obj.dem_x = obj.dem_x(obj.dem_x>=min(obj.ice_x) & obj.dem_x<=max(obj.ice_x));
      obj.dem_y = obj.dem_y(obj.dem_y>=min(obj.ice_y) & obj.dem_y<=max(obj.ice_y));
      
      obj.dem_x_mesh = repmat(obj.dem_x,length(obj.dem_y),1);
      obj.dem_y_mesh = repmat(obj.dem_y',1,length(obj.dem_x));
      obj.ice_x_mesh = repmat(obj.ice_x,length(obj.ice_y),1);
      obj.ice_y_mesh = repmat(obj.ice_y',1,length(obj.ice_x));
        
      obj.gray = rgb2gray(obj.dem);
      
      ice_mean = 242;
      ice_std = 14;
      rock_mean = 105;
      rock_std = 18;
      water_mean = 169;
      water_std = 15;
      obj.intensity_thresh_default = ice_mean-ice_std;
      obj.intensity_thresh = obj.intensity_thresh_default;
      
      obj.active_tool_idx = 0;
      obj.tool_list = [];
      tool_idx = 0;

      tool_idx = tool_idx + 1;
      obj.tool_list(tool_idx).event_key = 'm';
      obj.tool_list(tool_idx).shift_pressed = NaN;
      obj.tool_list(tool_idx).ctrl_pressed = NaN;
      obj.tool_list(tool_idx).type = 'mask';
      obj.tool_list(tool_idx).help_str = 'm: get ice mask';
      obj.tool_list(tool_idx).fh_callback = @obj.update_mask;
      
      obj.h_dem_fig = figure;
      set(obj.h_dem_fig,'DockControls','off');
      set(obj.h_dem_fig,'NumberTitle','off');
      set(obj.h_dem_fig,'ToolBar','none');
      set(obj.h_dem_fig,'MenuBar','none');
      set(obj.h_dem_fig,'Name','Satellite Image');
      
      obj.gui.left_panel = uipanel('parent',obj.h_dem_fig);
      obj.gui.right_panel = uipanel('parent',obj.h_dem_fig);
      obj.h_dem_axes = axes('Parent',obj.gui.right_panel,'YDir','normal');
      obj.h_dem_plot = imagesc(obj.dem_x,obj.dem_y,obj.dem,'Parent',obj.h_dem_axes);
      hold(obj.h_dem_axes,'on')
      
      obj.h_mask_fig = figure;
      set(obj.h_mask_fig,'DockControls','off');
      set(obj.h_mask_fig,'NumberTitle','off');
      set(obj.h_mask_fig,'ToolBar','none');
      set(obj.h_mask_fig,'MenuBar','none');
      set(obj.h_mask_fig,'Name','Ice Mask');
      
      obj.h_mask_axes = axes('Parent',obj.h_mask_fig,'YDir','normal');
      obj.h_mask_plot = imagesc(obj.ice_x,obj.ice_y,obj.mask);
      hold(obj.h_mask_axes,'on')
      
      set(obj.h_dem_fig,'WindowButtonUpFcn',@obj.button_up);
      set(obj.h_dem_fig,'WindowButtonDownFcn',@obj.button_down);
      set(obj.h_dem_fig,'WindowKeyPressFcn',@obj.key_press);
      set(obj.h_dem_fig,'WindowKeyReleaseFcn',@obj.key_release);
      
      % Set up zoom
      zoom_setup(obj.h_dem_fig);
      obj.zoom_mode = true;
      set(obj.h_dem_fig,'pointer','custom');
      
      zoom_setup(obj.h_mask_fig);
            
      obj.gui.table.ui = obj.h_dem_fig;
      obj.gui.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
      obj.gui.table.height_margin = NaN*zeros(30,30);
      obj.gui.table.false_width = NaN*zeros(30,30);
      obj.gui.table.false_height = NaN*zeros(30,30);
      obj.gui.table.offset = [0 0];
      row = 1;
      col = 1;
      obj.gui.table.handles{row,col}   = obj.gui.left_panel;
      obj.gui.table.width(row,col)     = 70;
      obj.gui.table.height(row,col)    = inf;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      
      row = 1;
      col = 2;
      obj.gui.table.handles{row,col}   = obj.gui.right_panel;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = inf;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      
      clear row col
      table_draw(obj.gui.table);
      
      obj.gui.threshSlider = uicontrol('parent',obj.gui.left_panel);
      set(obj.gui.threshSlider,'style','slider')
      set(obj.gui.threshSlider,'string','Theshold')
      set(obj.gui.threshSlider,'Min',0,'Max',255)
      set(obj.gui.threshSlider,'Value',obj.intensity_thresh_default);
      set(obj.gui.threshSlider,'Callback',@obj.update_threshold)
      set(obj.gui.threshSlider,'TooltipString','Changes intensity threshold for ice.');
      
      obj.gui.threshDefault = uicontrol('parent',obj.gui.left_panel);
      set(obj.gui.threshDefault,'style','pushbutton')
      set(obj.gui.threshDefault,'string','Default')
      set(obj.gui.threshDefault,'Callback',@obj.PB_default)
      
      obj.gui.threshDisp = uicontrol('parent',obj.gui.left_panel);
      set(obj.gui.threshDisp,'style','text')
      set(obj.gui.threshDisp,'String',sprintf('%0.0f',obj.intensity_thresh_default));
      
      
      obj.gui.left_table.ui =  obj.gui.left_panel;
      obj.gui.left_table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
      obj.gui.left_table.height_margin = NaN*zeros(30,30);
      obj.gui.left_table.false_width = NaN*zeros(30,30);
      obj.gui.left_table.false_height = NaN*zeros(30,30);
      obj.gui.left_table.offset = [0 0];
      
      row = 0;
      col = 0;
      
      row = row + 1;
      col = col + 1;
      obj.gui.left_table.handles{row,col}   = obj.gui.threshSlider;
      obj.gui.left_table.width(row,col)     = 20;
      obj.gui.left_table.height(row,col)    = inf;
      obj.gui.left_table.width_margin(row,col) = 1;
      obj.gui.left_table.height_margin(row,col) = 1;
      
      col = col + 1;
      obj.gui.left_table.handles{row,col} = obj.gui.threshDisp;
      obj.gui.left_table.width(row,col) = inf;
      obj.gui.left_table.height(row,col) = 20;
      obj.gui.left_table.width_margin(row,col) = 1;
      obj.gui.left_table.height_margin(row,col) = 1;
      
      row = row + 1;
      obj.gui.left_table.handles{row,col} = obj.gui.threshDefault;
      obj.gui.left_table.width(row,col) = inf;
      obj.gui.left_table.height(row,col) = 20;
      obj.gui.left_table.width_margin(row,col) = 1;
      obj.gui.left_table.height_margin(row,col) = 1;
      
      clear row col
      table_draw(obj.gui.left_table);
      
      if isfield(param,'mdata')
        obj.change_slice(1);
      end
      
    end
    
    
    function button_down(obj,h_obj,event)
      [obj.x,obj.y,but] = get_mouse_info(obj.h_dem_fig,obj.h_dem_axes);
      rbbox;
    end
    
    
    function button_up(obj,h_obj,event)
%       h_axes = get(h_obj,'Children');
      h_axes = obj.h_dem_axes;
      [x,y,but] = get_mouse_info(h_obj,h_axes);
      
      if obj.zoom_mode
        zoom_button_up(x,y,but,struct('x',obj.x,'y',obj.y, ...
        'h_axes',h_axes,'xlims',[min(obj.dem_x),max(obj.dem_x)],'ylims',[min(obj.dem_y),max(obj.dem_y)]));
        if strcmp(get(h_obj,'Name'),'Satellite Image')
          set(obj.h_mask_axes,'xlim',get(obj.h_dem_axes,'xlim'));
          set(obj.h_mask_axes,'ylim',get(obj.h_dem_axes,'ylim'));
        else
          set(obj.h_dem_axes,'xlim',get(obj.h_mask_axes,'xlim'));
          set(obj.h_dem_axes,'ylim',get(obj.h_mask_axes,'ylim'));
        end
        
      elseif but == 1 && x~=obj.x && y~=obj.y && obj.active_tool_idx == 1
        xlims = xlim(h_axes);
        ylims = ylim(h_axes);
        if x >= xlims(1) && x <= xlims(2) && y >= ylims(1) && y <= ylims(2)
          x_poly = [sort([x,obj.x]),sort([x,obj.x],'descend')];
          y_poly = [max(y,obj.y),max(y,obj.y),min(y,obj.y),min(y,obj.y)];
          
          update_mask(obj,x_poly,y_poly);
        end
      end
    end
    
    
    %% key_press
    function key_press(obj,src,event)
      
      if any(strcmp('shift',event.Modifier))
        obj.shift_pressed = true;
      else
        obj.shift_pressed = false;
      end
      
      if any(strcmp('control',event.Modifier))
        obj.ctrl_pressed = true;
      else
        obj.ctrl_pressed = false;
      end
      
      % Check to make sure that a key was pressed and not
      % just a modifier (e.g. shift, ctrl, alt)
      if ~isempty(event.Key)
        
        switch event.Key
            
          case 'z'
            % toggle zoom mode
            obj.zoom_mode = ~obj.zoom_mode;
            if obj.zoom_mode
              set(src,'pointer','custom');
            else
              set(src,'pointer','arrow');
            end
            
          case 'downarrow' % Down-arrow: Pan down
            zoom_arrow(event,struct('h_axes',obj.h_dem_axes, ...
              'xlims',[min(obj.dem_x),max(obj.dem_x)],'ylims',[min(obj.dem_y),max(obj.dem_y)]));
            set(obj.h_mask_axes,'xlim',get(obj.h_dem_axes,'xlim'));
            set(obj.h_mask_axes,'ylim',get(obj.h_dem_axes,'ylim'));
            
          case 'uparrow' % Up-arrow: Pan up
            zoom_arrow(event,struct('h_axes',obj.h_dem_axes, ...
              'xlims',[min(obj.dem_x),max(obj.dem_x)],'ylims',[min(obj.dem_y),max(obj.dem_y)]));
            set(obj.h_mask_axes,'xlim',get(obj.h_dem_axes,'xlim'));
            set(obj.h_mask_axes,'ylim',get(obj.h_dem_axes,'ylim'));
            
          case 'rightarrow' % Right arrow: Pan right
            zoom_arrow(event,struct('h_axes',obj.h_dem_axes, ...
              'xlims',[min(obj.dem_x),max(obj.dem_x)],'ylims',[min(obj.dem_y),max(obj.dem_y)]));
            set(obj.h_mask_axes,'xlim',get(obj.h_dem_axes,'xlim'));
            set(obj.h_mask_axes,'ylim',get(obj.h_dem_axes,'ylim'));
            
          case 'leftarrow' % Left arrow: Pan left
            zoom_arrow(event,struct('h_axes',obj.h_dem_axes, ...
              'xlims',[min(obj.dem_x),max(obj.dem_x)],'ylims',[min(obj.dem_y),max(obj.dem_y)]));
            set(obj.h_mask_axes,'xlim',get(obj.h_dem_axes,'xlim'));
            set(obj.h_mask_axes,'ylim',get(obj.h_dem_axes,'ylim'));
          
          case 'period'
            if ~obj.shift_pressed
              obj.slice = obj.slice + 1;
              obj.change_slice(obj.slice);
            else
              obj.slice = obj.slice + 10;
              obj.change_slice(obj.slice);
            end
            notify(obj,'SliceChange')
          
          case 'comma'
            if ~obj.shift_pressed
              obj.slice = obj.slice - 1;
              obj.change_slice(obj.slice);
            else
              obj.slice = obj.slice - 10;
              obj.change_slice(obj.slice);
            end
            notify(obj,'SliceChange')
            
          case 'p'
            if obj.active_tool_idx ~= 0
              obj.actions.getting_polygon = 1;
              poly_handle = impoly;
              position = wait(poly_handle);
              obj.actions.getting_polygon = 2;
              [polyPts] = getPosition(poly_handle);
              xPoly = polyPts(:,1);
              yPoly = polyPts(:,2);
              delete(poly_handle);
              obj.tool_list(obj.active_tool_idx).fh_callback(xPoly,yPoly);
            end
          case 'm'
%             if obj.active_tool_idx == 1;
%               set(src,'pointer','custom');
%               obj.active_tool_idx = 0;
%               obj.zoom_mode = 1;
%             else
              set(src,'pointer','arrow');
              obj.active_tool_idx = 1;
              obj.zoom_mode = 0;
%             end
            
          case 's'
            obj.change_slice(2000);

          case 'r'
            obj.local_undo_stack.redo();
            notify(obj,'Redo');
            
          case 'u'
            obj.local_undo_stack.pop();
            notify(obj,'Undo')
            
          otherwise
            
        end
        
      end
    end
    
    
    function key_release(obj,src,event)
      
      if any(strcmp('shift',event.Modifier))
        obj.shift_pressed = true;
      else
        obj.shift_pressed = false;
      end
      
      if any(strcmp('control',event.Modifier))
        obj.ctrl_pressed = true;
      else
        obj.ctrl_pressed = false;
      end
    end
    
    
    function update_mask(obj,xv,yv)
      
      x_p_lim = [min(xv),max(xv)];
      y_p_lim = [min(yv),max(yv)];
      
      gray_idx = find(obj.dem_y_mesh>=y_p_lim(1) & obj.dem_y_mesh<=y_p_lim(2) &...
        obj.dem_x_mesh>=x_p_lim(1) & obj.dem_x_mesh<=x_p_lim(2));
      
      gray_tmp = obj.gray(gray_idx);
%       gray_tmp = rgb2gray(obj.dem(gray_idx,:));
      x_tmp = obj.dem_x_mesh(gray_idx);
      y_tmp = obj.dem_y_mesh(gray_idx);
      
%       diff_x = abs(mean(diff(obj.dem_x)));
%       diff_y = abs(mean(diff(obj.dem_y)));
      
%       ice_block_idx = find(obj.ice_y_mesh>=y_p_lim(1) & obj.ice_y_mesh<=y_p_lim(2) &...
%         obj.ice_x_mesh>=x_p_lim(1) & obj.ice_x_mesh<=x_p_lim(2));
      ice_block_idx = find(obj.ice_y_mesh>=min(y_tmp) & obj.ice_y_mesh<=max(y_tmp) &...
        obj.ice_x_mesh>=min(x_tmp) & obj.ice_x_mesh<=max(x_tmp));
      
      if isempty(ice_block_idx)
        return
      end
%       ice_block_idx = find(obj.ice_y_mesh>=min(y_tmp)-diff_y/2 & obj.ice_y_mesh<=max(y_tmp)+diff_y/2 &...
%         obj.ice_x_mesh>=min(x_tmp)-diff_x/2 & obj.ice_x_mesh<=max(x_tmp)+diff_x/2);


      ice_x_tmp = obj.ice_x_mesh(ice_block_idx);
      ice_y_tmp = obj.ice_y_mesh(ice_block_idx);
      [ice_x_mesh,ice_y_mesh] = meshgrid(sort(unique(ice_x_tmp)),sort(unique(ice_y_tmp),'descend'));
%       ice_x_mesh = repmat(sort(unique(ice_x_tmp')),length(unique(ice_y_tmp)),1);
%       ice_y_mesh = repmat(sort(unique(ice_y_tmp),'descend'),1,length(unique(ice_x_tmp)));
      
%       inter_idx = find(obj.intersections(2,:,:)>=y_p_lim(1) & obj.intersections(2,:,:)<=y_p_lim(2) &...
%         obj.intersections(1,:,:)>=x_p_lim(1) & obj.intersections(1,:,:)<=x_p_lim(2));
      diff_x = abs(mean(diff(obj.ice_x)));
      diff_y = abs(mean(diff(obj.ice_y)));
            
      if length(ice_x_tmp)~=numel(ice_x_mesh)
        keyboard;
      end

      gray_tmp = griddata(x_tmp,y_tmp,double(gray_tmp),ice_x_mesh,ice_y_mesh);
      [ice_in] = inpolygon(ice_x_mesh,ice_y_mesh,xv,yv);
      gray_tmp_in = gray_tmp(ice_in);
      mask_tmp = gray_tmp_in>=obj.intensity_thresh;
      
      if all(obj.mask(ice_block_idx(ice_in))==mask_tmp)
        return
      end
      
      inter_idx = [];
      cmd = [];
      cmd{1}.redo.data_mask = [];
      cmd{1}.undo.data_mask = [];
      cmd{1}.undo.data_mask_idx = [];
      cmd{1}.redo.data_mask_idx = [];
      if ~isempty(obj.intersections)
        inter_idx = find(obj.intersections(2,:,:)>=min(ice_y_tmp)-diff_y/2 & obj.intersections(2,:,:)<=max(ice_y_tmp)+diff_y/2 &...
          obj.intersections(1,:,:)>=min(ice_x_tmp)-diff_x/2 & obj.intersections(1,:,:)<=max(ice_x_tmp)+diff_x/2);
        inter_idx = inter_idx(inpolygon(obj.intersections(1,inter_idx),obj.intersections(2,inter_idx),xv,yv));
      end
      if ~isempty(inter_idx)
%         ice_in_buff = logical(conv2(double(ice_in),ones(3),'same'));
%         inter_ice_val = griddata(ice_x_tmp(ice_in_buff),ice_y_tmp(ice_in_buff),double(mask_tmp),obj.intersections(1,inter_idx),obj.intersections(2,inter_idx),'nearest');
        inter_ice_val = griddata(ice_x_tmp(ice_in),ice_y_tmp(ice_in),double(mask_tmp),obj.intersections(1,inter_idx),obj.intersections(2,inter_idx),'nearest');

        inter_tmp_idx = ~isnan(inter_ice_val);
        inter_diff_idx = inter_idx(inter_tmp_idx);
        inter_val = inter_ice_val(inter_tmp_idx);

        data_mask_tmp = obj.mdata.ice_mask;
        data_mask_tmp(inter_diff_idx) = inter_val;
        cmd{1}.redo.data_mask = data_mask_tmp(inter_diff_idx);
        cmd{1}.undo.data_mask = obj.mdata.ice_mask(inter_diff_idx);
        cmd{1}.undo.data_mask_idx = inter_diff_idx;
        cmd{1}.redo.data_mask_idx = inter_diff_idx;
      end
      
      cmd{1}.undo.mask = obj.mask(ice_block_idx(ice_in));
      cmd{1}.redo.mask = mask_tmp;
      cmd{1}.undo.mask_idx = ice_block_idx(ice_in);
      cmd{1}.redo.mask_idx = ice_block_idx(ice_in);

      cmd{1}.type = 'ice_mask';
      
      if obj.local_undo_flag
        obj.local_undo_stack.push(cmd);
      else
        obj.cmd = cmd;
        notify(obj,'IceChange');
      end
      obj.cmd = [];
      
%       obj.mask(ice_block_idx(ice_in)) = mask_tmp;
      
    end
    
    
    function change_slice(obj,slice)
      if isempty(obj.mdata_loaded) || ~obj.mdata_loaded
        %% load data
        fprintf('Aligning TWTT data\n')
        % convert from FCS to proj
        origin_ecef = obj.mdata.param_combine.array_param.fcs{1}{1}.origin(:,:);
        physical_constants;
        [origin_lat,origin_lon] = ecef2geodetic(origin_ecef(1,:),origin_ecef(2,:),origin_ecef(3,:),WGS84.ellipsoid);
        origin_lat = origin_lat*180/pi;
        origin_lon = origin_lon*180/pi;
        % Convert from geodetic to projection
        obj.flight_line = zeros(2,size(obj.mdata.twtt,2));
        [obj.flight_line(1,:),obj.flight_line(2,:)] = projfwd(obj.proj,origin_lat,origin_lon);
        
        obj.h_flight_dem_plot = plot(obj.flight_line(1,:),obj.flight_line(2,:),'b','Parent',obj.h_dem_axes);
        obj.h_flight_mask_plot = plot(obj.flight_line(1,:),obj.flight_line(2,:),'b','Parent',obj.h_mask_axes);
        obj.h_true_dem_plot = plot(NaN,NaN,'b.','Parent',obj.h_dem_axes);
        obj.h_false_dem_plot = plot(NaN,NaN,'r.','Parent',obj.h_dem_axes);
        obj.h_true_mask_plot = plot(NaN,NaN,'b.','Parent',obj.h_mask_axes);
        obj.h_false_mask_plot = plot(NaN,NaN,'r.','Parent',obj.h_mask_axes);
        
%         intersection = zeros(2,length(obj.mdata.theta),size(obj.mdata.twtt,2));
        Nr = size(obj.mdata.twtt,2);
        Nt = length(obj.mdata.theta);
        
        intersect_ecef_all = zeros([3,Nt,Nr]);
        for rline = 1:Nr
          
          Tfcs_ecef = [obj.mdata.param_combine.array_param.fcs{1}{1}.x(:,rline), ...
            obj.mdata.param_combine.array_param.fcs{1}{1}.y(:,rline), ...
            obj.mdata.param_combine.array_param.fcs{1}{1}.z(:,rline)];
          
          dir = [zeros(1,Nt) ; sin(obj.mdata.theta) ; -cos(obj.mdata.theta)];
          for row = 1:Nt
            dir(:,row) = dir(:,row)/norm(dir(:,row));
          end
          
          intersect =  dir .* repmat(obj.mdata.twtt(:,rline)',3,1) * 3e8/2;
          
          intersect_ecef = Tfcs_ecef*intersect;
          intersect_ecef_all(1,:,rline) = intersect_ecef(1,:) + origin_ecef(1,rline);
          intersect_ecef_all(2,:,rline) = intersect_ecef(2,:) + origin_ecef(2,rline);
          intersect_ecef_all(3,:,rline) = intersect_ecef(3,:) + origin_ecef(3,rline);
                    
        end
        
        [intersect_lat,intersect_lon] = ecef2geodetic(intersect_ecef_all(1,:,:),intersect_ecef_all(2,:,:),intersect_ecef_all(3,:,:),WGS84.ellipsoid);
        intersect_lat = intersect_lat*180/pi;
        intersect_lon = intersect_lon*180/pi;
        [intersect_x,intersect_y] = projfwd(obj.proj,intersect_lat,intersect_lon);

        obj.intersections = [intersect_x;intersect_y];
        
        obj.mdata_loaded = 1;
              
%         obj.slice = slice;
%         intersection = obj.intersections(:,:,slice);
        
        set(obj.h_dem_axes,'xlim',[min(min(obj.intersections(1,:,:))) - 5000, ...
          max(max(obj.intersections(1,:,:))) + 5000]);
        set(obj.h_dem_axes,'ylim',[min(min(obj.intersections(2,:,:))) - 5000, ...
          max(max(obj.intersections(2,:,:))) + 5000]);
        set(obj.h_mask_axes,'xlim',[min(min(obj.intersections(1,:,:))) - 5000, ...
          max(max(obj.intersections(1,:,:))) + 5000]);
        set(obj.h_mask_axes,'ylim',[min(min(obj.intersections(2,:,:))) - 5000, ...
          max(max(obj.intersections(2,:,:))) + 5000]);
        
        fprintf('Finished Aligning Data\n');
      end
      
      if slice < 1
        slice = 1;
      elseif slice > size(obj.intersections,3)
        slice = size(obj.intersections,3);
      end
      
      obj.slice = slice;
      intersection = obj.intersections(:,:,slice);
      mask_tmp = logical(obj.mdata.ice_mask(:,slice));
%       obj.h_slice_plot = plot(intersection(1,:),intersection(2,:),'m.','Parent',obj.h_dem_axes);
      set(obj.h_true_dem_plot,'XData',intersection(1,mask_tmp),'YData',intersection(2,mask_tmp));
      set(obj.h_false_dem_plot,'XData',intersection(1,~mask_tmp),'YData',intersection(2,~mask_tmp));
      set(obj.h_true_mask_plot,'XData',intersection(1,mask_tmp),'YData',intersection(2,mask_tmp));
      set(obj.h_false_mask_plot,'XData',intersection(1,~mask_tmp),'YData',intersection(2,~mask_tmp));
    end
    
    function getEventData(obj,src,~)
      obj.change_slice(src.slice);
    end
    
    function undo_sync(obj,src,~)
      [cmds_list,cmds_direction] =  src.get_synchronize_cmds();
      if strcmp(cmds_direction,'redo')
        for cmd_idx = 1:length(cmds_list)
          for subcmd_idx = 1:length(cmds_list{cmd_idx})
            if strcmp(cmds_list{cmd_idx}{subcmd_idx}.type,'ice_mask')
              if isfield(cmds_list{cmd_idx}{subcmd_idx}.redo,'mask')
                obj.mask(cmds_list{cmd_idx}{subcmd_idx}.redo.mask_idx) = ...
                  cmds_list{cmd_idx}{subcmd_idx}.redo.mask;
              end
              if isfield(cmds_list{cmd_idx}{subcmd_idx}.redo,'data_mask')
                obj.mdata.ice_mask(cmds_list{cmd_idx}{subcmd_idx}.redo.data_mask_idx) ...
                  = cmds_list{cmd_idx}{subcmd_idx}.redo.data_mask;
              end
            end
          end
        end
      else
        for cmd_idx = 1:length(cmds_list)
          for subcmd_idx = 1:length(cmds_list{cmd_idx})
            if strcmp(cmds_list{cmd_idx}{subcmd_idx}.type,'ice_mask')
              if isfield(cmds_list{cmd_idx}{subcmd_idx}.undo,'mask')
                obj.mask(cmds_list{cmd_idx}{subcmd_idx}.undo.mask_idx) = ...
                  cmds_list{cmd_idx}{subcmd_idx}.undo.mask;
              end
              if isfield(cmds_list{cmd_idx}{subcmd_idx}.undo,'data_mask')
                obj.mdata.ice_mask(cmds_list{cmd_idx}{subcmd_idx}.undo.data_mask_idx) ...
                  = cmds_list{cmd_idx}{subcmd_idx}.undo.data_mask;
              end
            end
          end
        end
      end
      
      set(obj.h_mask_plot,'CData',obj.mask);
      intersection = obj.intersections(:,:,obj.slice);
      if ~isempty(obj.mdata_loaded) && obj.mdata_loaded
        twtt_mask_tmp = logical(obj.mdata.ice_mask(:,obj.slice));
        set(obj.h_true_dem_plot,'XData',intersection(1,twtt_mask_tmp),'YData',intersection(2,twtt_mask_tmp));
        set(obj.h_false_dem_plot,'XData',intersection(1,~twtt_mask_tmp),'YData',intersection(2,~twtt_mask_tmp));
        set(obj.h_true_mask_plot,'XData',intersection(1,twtt_mask_tmp),'YData',intersection(2,twtt_mask_tmp));
        set(obj.h_false_mask_plot,'XData',intersection(1,~twtt_mask_tmp),'YData',intersection(2,~twtt_mask_tmp));
      end
    end    
    
    function update_threshold(obj,source,~)
      val = source.Value;
      obj.intensity_thresh = val;
      set(obj.gui.threshDisp,'String',sprintf('%0.0f',val));
    end
    
    function PB_default(obj,source,~)
      val = obj.intensity_thresh_default;
      set(obj.gui.threshSlider,'Value',val);
      set(obj.gui.threshDisp,'String',sprintf('%0.0f',val));
      obj.intensity_thresh = val;
    end
    
    function cmd = edit_twtt(obj,theta,slice,val,val_old)
      x = obj.intersections(1,theta,slice);
      y = obj.intersections(2,theta,slice);

      ice_x = interp1(obj.ice_x,1:length(obj.ice_x),x,'nearest');
      ice_y = interp1(obj.ice_y,1:length(obj.ice_y),y,'nearest');
      ice_idx = sub2ind(size(obj.ice_x_mesh),ice_y,ice_x);

      ice_x_idx = [ice_x-1,ice_x,ice_x+1];
      ice_x_idx = ice_x_idx(ice_x_idx>0 & ice_x_idx<=length(obj.ice_x));
      ice_y_idx = [ice_y-1,ice_y,ice_y+1];
      ice_y_idx = ice_y_idx(ice_y_idx>0 & ice_y_idx<=length(obj.ice_x));
      [ice_x_idx,ice_y_idx] = meshgrid(ice_x_idx,ice_y_idx);
      ice_idx_exp = sub2ind(size(obj.ice_x_mesh),ice_y_idx,ice_x_idx);

      ice_x_mesh = obj.ice_x_mesh(ice_idx_exp);
      ice_y_mesh = obj.ice_y_mesh(ice_idx_exp);
      idx_mesh = reshape(1:numel(ice_x_mesh),size(ice_x_mesh));
      idx_match = griddata(ice_x_mesh,ice_y_mesh,idx_mesh,x,y,'nearest');

      inter_near_idx = find(obj.intersections(1,:,:) >= min(min(ice_x_mesh)) & ...
        obj.intersections(1,:,:) <= max(max(ice_x_mesh)) & ...
        obj.intersections(2,:,:) >= min(min(ice_y_mesh)) & ...
        obj.intersections(2,:,:) <= max(max(ice_y_mesh)));

      inter_tmp = obj.intersections(:,inter_near_idx);

      g = griddata(ice_x_mesh,ice_y_mesh,idx_mesh,inter_tmp(1,:),inter_tmp(2,:),'nearest');
%           g = interp2(ice_x_mesh,ice_y_mesh,idx_mesh,inter_tmp(1,:),inter_tmp(2,:),'nearest');

      inter_eff_idx = inter_near_idx(ismember(g,idx_match));
      
      cmd{1}.redo.data_mask = val*ones(1,length(inter_eff_idx));
      cmd{1}.undo.data_mask = val_old*ones(1,length(inter_eff_idx));
      cmd{1}.undo.data_mask_idx = inter_eff_idx;
      cmd{1}.redo.data_mask_idx = inter_eff_idx;

      cmd{1}.undo.mask = val_old;
      cmd{1}.redo.mask = val;
      cmd{1}.undo.mask_idx = ice_idx;
      cmd{1}.redo.mask_idx = ice_idx;
      cmd{1}.type = 'ice_mask';
      
    end
    
  end
  
end