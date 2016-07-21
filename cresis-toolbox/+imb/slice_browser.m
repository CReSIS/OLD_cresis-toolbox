classdef slice_browser < handle
  %slice_browser: Browse slices of a 3D image
  %
  % Example:
  % mdata = load('C:\tmp\Canada_Tomography\Data_img_01_20140401_03_044.mat');
  % h_control_fig = figure;
  % h_control_axes = axes('Parent',h_control_fig);
  % h_control_image = imagesc(lp(squeeze(mdata.Topography.img(:,33,:))),'Parent',h_control_axes);
  % try; delete(obj); end;
  % obj = slice_browser(mdata.Topography.img,h_control_image);
  
  properties
    select_mask
    data % N-dimensional matrix with last dimension equal to Nx
    slice % Integer from 1 to Nx
    layer % Layer structures
    layer_fn
    
    % GUI handles
    h_control_fig
    h_control_axes
    h_control_image
    h_control_plot

    h_fig
    h_axes
    h_image
    
    gui
    
    % Function handle hooks for customizing clip_matrix
    fh_button_up
    fh_key_press
    fh_button_motion
    
    % zoom_mode: boolean, x,y: used by zoom mode
    zoom_mode
    x
    y
    shift_pressed
    ctrl_pressed
  end
  
  methods
    %% constructor/slice_browser:
    function obj = slice_browser(data,h_control_image,param)
      if ~exist('param','var')
        param = [];
      end
      if ~isfield(param,'fh_button_up')
        param.fh_button_up = [];
      end
      if ~isfield(param,'fh_key_press')
        param.fh_key_press = [];
      end
      if ~isfield(param,'fh_button_motion')
        param.fh_button_motion = [];
      end
      
      obj.data = data;
      obj.slice = 1;
      
      % Load layer data
      if isfield(param,'layer_fn') && ~isempty(param.layer_fn)
        tmp = load(param.layer_fn);
        obj.layer_fn = param.layer_fn;
        obj.layer = tmp.layer;
      else
        obj.layer = [];
        obj.layer.x = [];
        obj.layer.y  = [];
        obj.layer.name = '';
        obj.layer.plot_name_values = [];
      end
      
      obj.h_control_image = h_control_image;
      obj.h_control_axes = get(obj.h_control_image,'Parent');
      obj.h_control_fig = get(obj.h_control_axes,'Parent');
      obj.fh_button_up = param.fh_button_up;
      obj.fh_key_press = param.fh_key_press;
      obj.fh_button_motion = param.fh_button_motion;
      
      obj.h_fig = figure;
      obj.gui.panel = uipanel('parent',obj.h_fig);
      obj.h_axes = axes('Parent',obj.gui.panel,'YDir','reverse');
      hold(obj.h_axes,'on');
      
      obj.h_image = imagesc(obj.data(:,:,obj.slice),'parent',obj.h_axes);
      for layer_idx = 1:numel(obj.layer)
        obj.layer(layer_idx).h_plot ...
          = plot(obj.layer(layer_idx).x(:,obj.slice), ...
          obj.layer(layer_idx).y(:,obj.slice), ...
          'parent',obj.h_axes,'color','black', ...
          obj.layer(layer_idx).plot_name_values{:});
      end
      obj.gui.h_select_plot = plot(NaN,NaN,'m.');
      
      set(obj.h_control_fig, 'WindowButtonUpFcn', @obj.control_button_up);
      
      hold(obj.h_control_axes,'on');
      obj.h_control_plot = plot(NaN,NaN,'parent',obj.h_control_axes,'Marker','x','Color','black','LineWidth',2,'MarkerSize',10);
      
      % Set figure call back functions
      set(obj.h_fig,'WindowButtonUpFcn',@obj.button_up);
      set(obj.h_fig,'WindowButtonDownFcn',@obj.button_down);
      set(obj.h_fig,'WindowButtonMotionFcn',@obj.button_motion);
      set(obj.h_fig,'WindowScrollWheelFcn',@obj.button_scroll);
      set(obj.h_fig,'WindowKeyPressFcn',@obj.key_press);
      set(obj.h_fig,'WindowKeyReleaseFcn',@obj.key_release);
      set(obj.h_fig,'CloseRequestFcn',@obj.close_win);
      
      % Set up zoom
      zoom_setup(obj.h_fig);
      obj.zoom_mode = true;
      set(obj.h_fig,'pointer','custom');
      
      obj.select_mask = logical(zeros(size(obj.data,2),1));
      
      % Set limits to size of data
      xlim([1 size(obj.data(:,:,obj.slice),2)]);
      ylim([1 size(obj.data(:,:,obj.slice),1)]);
      
      obj.gui.nextPB = uicontrol('parent',obj.h_fig);
      set(obj.gui.nextPB,'style','pushbutton')
      set(obj.gui.nextPB,'string','>(p)')
      set(obj.gui.nextPB,'Callback',@obj.Next_button_callback)
      
      obj.gui.prevPB = uicontrol('parent',obj.h_fig);
      set(obj.gui.prevPB,'style','pushbutton')
      set(obj.gui.prevPB,'string','<(n)')
      set(obj.gui.prevPB,'Callback',@obj.Prev_button_callback)
      
      obj.gui.savePB = uicontrol('parent',obj.h_fig);
      set(obj.gui.savePB,'style','pushbutton')
      set(obj.gui.savePB,'string','save')
      set(obj.gui.savePB,'Callback',@obj.save_button_callback)
      
      obj.gui.helpPB = uicontrol('parent',obj.h_fig);
      set(obj.gui.helpPB,'style','pushbutton')
      set(obj.gui.helpPB,'string','?')
      set(obj.gui.helpPB,'Callback',@obj.Help_button_callback)
      
      obj.gui.layerPM = uicontrol('parent',obj.h_fig);
      set(obj.gui.layerPM,'style','popup')
      set(obj.gui.layerPM,'string',{obj.layer.name})
      %set(obj.gui.layerPM,'Callback',@obj.layerPM_callback)
      
      obj.gui.plotPM = uicontrol('parent',obj.h_fig);
      set(obj.gui.plotPM,'style','popup')
      set(obj.gui.plotPM,'string',{'image','plot','scatter'})
      %set(obj.gui.plotPM,'Callback',@plotPM_callback)
      
      obj.gui.variablePM = uicontrol('parent',obj.h_fig);
      set(obj.gui.variablePM,'style','popup')
      set(obj.gui.variablePM,'string',{'data'})
      %set(obj.gui.variablePM,'Callback',@variablePM_callback)
      
      obj.gui.layerTXT = uicontrol('Style','text','string','layer');
      obj.gui.plotTXT = uicontrol('Style','text','string','plot');
      obj.gui.variableTXT = uicontrol('Style','text','string','variable');
      
      
      %% Create GUI Table
      obj.gui.table.ui = obj.h_fig;
      obj.gui.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
      obj.gui.table.height_margin = NaN*zeros(30,30);
      obj.gui.table.false_width = NaN*zeros(30,30);
      obj.gui.table.false_height = NaN*zeros(30,30);
      obj.gui.table.offset = [0 0];
      
      row = 0;
      col = 0;
      
      row = row + 1;
      col = col + 1;
      obj.gui.table.handles{row,col}   = obj.gui.nextPB;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      col = col + 1;
      obj.gui.table.handles{row,col}   = obj.gui.prevPB;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      col = col + 1;
      obj.gui.table.handles{row,col}   = obj.gui.layerTXT;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      col = col + 1;
      obj.gui.table.handles{row,col}   = obj.gui.plotTXT;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      
      col = col + 1;
      obj.gui.table.handles{row,col}   = obj.gui.variableTXT;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      
      col = 0
      row = row + 1;
      col = col + 1;
      obj.gui.table.handles{row,col}   = obj.gui.savePB;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      
      col = col + 1;
      obj.gui.table.handles{row,col}   = obj.gui.helpPB;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      
      col = col + 1;
      obj.gui.table.handles{row,col}   = obj.gui.layerPM;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      
      col = col + 1;
      obj.gui.table.handles{row,col}   = obj.gui.plotPM;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      
      col = col + 1;
      obj.gui.table.handles{row,col}   = obj.gui.variablePM;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      
      row = row + 1;
      col = 1;
      obj.gui.table.handles{row,col}   = obj.gui.panel;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = inf;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 5;
      
      clear row col
      table_draw(obj.gui.table);
      
      % Connect listeners to events
      
    end
    %% destructor/delete
    function delete(obj)
      delete(obj.h_fig)
    end
    %% close_win
    function close_win(obj,h_obj,event)
      try
        delete(obj);
      end
    end
    %% Buttonstuff
    function Next_button_callback(obj,source,callbackdata)
      obj.slice = obj.slice + 1;
      obj.update_slice();
    end
    function Prev_button_callback(obj,source,callbackdata)
      obj.slice = obj.slice -1;
      obj.update_slice();
    end
    function Help_button_callback(obj,source,callbackdata)
      obj.help_menu()
    end  
    function save_button_callback(obj,source,callbackdata)
      layer = obj.layer;
      save(param.layer_fn,'layer')  
    end
    %% control_button_up
    function control_button_up(obj,h_obj,event)
      [x,y,but] = get_mouse_info(obj.h_control_fig,obj.h_control_axes);
      
      obj.slice = ceil(x);
      obj.update_slice();
      
    end
    %% button_down
    function button_down(obj,h_obj,event)
      [obj.x,obj.y,but] = get_mouse_info(obj.h_fig,obj.h_axes);
      fprintf('Button Down: x = %.3f, y = %.3f, but = %d\n', obj.x, obj.y, but); % DEBUG ONLY
      rbbox;
    end
    %% button_up
    function button_up(obj,h_obj,event)
      % Run user defined button up callback
      if ~isempty(obj.fh_button_up)
        status = obj.fh_button_up(obj,h_obj,event);
        if status == 0
          return;
        end
      end
      
      % Get x,y position of user button release
      [x,y,but] = get_mouse_info(obj.h_fig,obj.h_axes);
      fprintf('Button Up: x = %.3f, y = %.3f, but = %d\n', x, y, but); % DEBUG ONLY
      
      layer_idx = get(obj.gui.layerTXT,'value');
                  
      if obj.zoom_mode
        zoom_button_up(x,y,but,struct('x',obj.x,'y',obj.y, ...
          'h_axes',obj.h_axes,'xlims',[1 size(obj.data,2)],'ylims',[1 size(obj.data,1)]));
      else
        if but == 2 || but == 3
          if obj.x == x
            obj.select_mask(round(x)) = true;
          else
            obj.shift_pressed
            if ~obj.shift_pressed
              obj.select_mask(:) = false;
              obj.update_slice;
            end

            obj.select_mask = obj.select_mask | (obj.layer(layer_idx).x(:,obj.slice) >= min(x,obj.x) ...
              & obj.layer(layer_idx).x(:,obj.slice) <= max(x,obj.x) ...
              & obj.layer(layer_idx).y(:,obj.slice) >= min(y,obj.y) ...
              & obj.layer(layer_idx).y(:,obj.slice) <= max(y,obj.y));
%             if x > obj.x;
%               for k = 1:size(obj.data,2);
%                 if k > obj.x && k < x;
%                   obj.select_mask(k,1) = true;
%                 end
%               end
%             
%             else
%               for k = 1:size(obj.data,2);
%                 if k < obj.x && k > x;
%                   obj.select_mask(k,1) = true;
%                 end
%               end
%             end
%             if y > obj.y;
%               for k = 1:size(obj.data,2);
%                 layer_idx = 1;
%                 if obj.layer(layer_idx).y(k,obj.slice) < obj.y && obj.layer(layer_idx).y(k,obj.slice) > y;
%                   obj.select_mask(k,1) = false;
%                 end
%               end
%             else
%               for k = 1:size(obj.data,2);
%                 if obj.layer(layer_idx).y(k,obj.slice) > obj.y && obj.layer(layer_idx).y(k,obj.slice) < y;
%                   obj.mask_select(k,1) = false;
%                 end
%               end
%             end
          end
        else
           layer_idx = get(obj.gui.layerPM,'value');
          obj.layer(layer_idx).y(round(x),obj.slice) = y;
        end
        obj.update_slice();
       end
      end
    
    %% button_motion
    function button_motion(obj,hObj,event)
      % Run user defined button up callback
      if ~isempty(obj.fh_button_motion)
        status = obj.fh_button_motion(obj,h_obj,event);
        if status == 0
          return;
        end
      end
      
      [x,y,but] = get_mouse_info(obj.h_fig,obj.h_axes);
      set(obj.h_control_plot,'x',obj.slice,'y',y);
    end
   
    %% button_scroll
    function button_scroll(obj,h_obj,event)
      zoom_button_scroll(event,struct('h_fig',obj.h_fig, ...
        'h_axes',obj.h_axes,'xlims',[1 size(obj.data,2)],'ylims',[1 size(obj.data,1)]));
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
      
      if ~isempty(obj.fh_key_press)
        status = obj.fh_key_press(src,event);
        if status == 0
          return;
        end
      end
      
      % Check to make sure that a key was pressed and not
      % just a modifier (e.g. shift, ctrl, alt)
      if ~isempty(event.Key)
        
        % see event.Modifier for modifiers
        switch event.Key
          
          case 'f1'
            obj.help_menu()
            
          case 'z'
            if obj.ctrl_pressed
              %% zoom reset
              axis tight;
            else
              % toggle zoom mode
              obj.zoom_mode = ~obj.zoom_mode;
              if obj.zoom_mode
                set(obj.h_fig,'pointer','custom');
              else
                set(obj.h_fig,'pointer','arrow');
              end
            end
            
          case 'downarrow' % Down-arrow: Pan down
            zoom_arrow(event,struct('h_axes',obj.h_axes, ...
              'xlims',[1 size(obj.data,2)],'ylims',[1 size(obj.data,1)]));
            
          case 'uparrow' % Up-arrow: Pan up
            zoom_arrow(event,struct('h_axes',obj.h_axes, ...
              'xlims',[1 size(obj.data,2)],'ylims',[1 size(obj.data,1)]));
            
          case 'rightarrow' % Right arrow: Pan right
            zoom_arrow(event,struct('h_axes',obj.h_axes, ...
              'xlims',[1 size(obj.data,2)],'ylims',[1 size(obj.data,1)]));
            
          case 'leftarrow' % Left arrow: Pan left
            zoom_arrow(event,struct('h_axes',obj.h_axes, ...
              'xlims',[1 size(obj.data,2)],'ylims',[1 size(obj.data,1)]));
            
          case 'n'
            obj.slice = obj.slice + 1;
            obj.update_slice();
            
          case 'p'
            obj.slice = obj.slice - 1;
            obj.update_slice();
            
          case 'delete'
            layer_idx = 1;
            for k = 1:64;
              if obj.select_mask(k,1) == 1;
                obj.layer(layer_idx).y(k,obj.slice) = NaN;
              end
            end
            obj.update_slice()
            obj.select_mask = logical(zeros(size(obj.data,2),1));
            
          otherwise
            
        end
        
      end
      
      if ~isempty(obj.fh_key_press)
        obj.fh_key_press(src,event)
      end
    end
    
    %% key_release
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
  
    %% Update slice
    function update_slice(obj)
      if obj.slice <= 0
        obj.slice = 1;
      end
      if obj.slice > size(obj.data,3)
        obj.slice = size(obj.data,3);
      end
      
      set(obj.h_image,'CData',obj.data(:,:,obj.slice));
      
      for layer_idx = 1:numel(obj.layer)
        set(obj.layer(layer_idx).h_plot, ...
          'XData', obj.layer(layer_idx).x(:,obj.slice), ...
          'YData', obj.layer(layer_idx).y(:,obj.slice));
        title(sprintf('Slice:%d',obj.slice),'parent',obj.h_axes)
        
      end
      layer_idx = get(obj.gui.layerPM,'value');
      x_select = obj.layer(layer_idx).x(:,obj.slice);
      y_select = obj.layer(layer_idx).y(:,obj.slice);
      set(obj.gui.h_select_plot,'XData',x_select(obj.select_mask), ...
        'YData',y_select(obj.select_mask));
      
      [x,y,but] = get_mouse_info(obj.h_fig,obj.h_axes);
      set(obj.h_control_plot,'x',obj.slice,'y',y);
    end
    
    %% Help
    function help_menu(obj)
      fprintf('Key Short Cuts\n');
      
      fprintf('? Mode\n');
      fprintf('scroll: zoom in/out at point\n');
      
      fprintf('Zoom Mode\n');
      fprintf('left-click and drag: zoom to selection\n');
      fprintf('left-click: zoom in at point\n');
      fprintf('right-click: zoom out at point\n');
      fprintf('scroll: zoom in/out at point\n');
    end
  end
  
end

