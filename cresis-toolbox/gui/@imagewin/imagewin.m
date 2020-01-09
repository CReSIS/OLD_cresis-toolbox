classdef imagewin < handle
  % imagewin: A class for controlling histogram and detrending image
  % handles. Works with linear power image data.
  %
  % mdata = load('C:\tmp\rds\2011_Greenland_P3\CSARP_mvdr\20110407_06\Data_20110407_06_002.mat');
  % figure(1); clf;
  % h_image = imagesc(mdata.Data);
  % img = imagewin('1: Image Params', h_image, false);
  
  properties
    h_fig
    h_gui
    img_title
    img
    hide_only
    CData
    
    img_fn % default filename (stores last filename saved or opened, initially empty)
    img_fn_save_formats % Cell vector for ui_put_file
    img_fn_open_formats % Cell vector for ui_get_file
    % img_fn_save_fh and img_fn_open_fh format:
    %   input: current filename string
    %   output: filename string
    img_fn_save_fh % Function handle to get image to save default filename (optional)
    img_fn_open_fh % Function handle to get image to open default filename (optional)
    % save_fn and open_fn format:
    %   input: filename string
    %   output: scalar. If >= 1, function will continue regular imagewin
    %     save/open function.
    save_fh % Function handle to save files (optional)
    open_fh % Function handle to open files (optional)
    
    save_mat_fh % Function handle for custom save
  end
  
  properties (SetAccess = private, GetAccess = private)
  end
    
  events
    open_event % Signalled when a user opens a file
    save_event % Signalled when a user saves a file
  end
  
  methods
    function obj = imagewin(img_title, img, hide_only)
      % img_title = string containing figure titlebar title
      % img = handle to matlab image object
      % hide_only = logical, when true, closing window just hides the
      %   window
      
      %%% Pre Initialization %%%
      % Any code not using output argument (obj)
      
      %%% Post Initialization %%%
      % Any code, including access to object
      obj.img = -1;
      obj.img_title = img_title;
      obj.hide_only = hide_only;
      
      % Set Save/Open properties
      obj.img_fn = '';
      obj.img_fn_save_formats = {'*.png', 'All png files (*.png)'; ...
        '*.*', 'All Files (*.*)'};
      obj.img_fn_open_formats = {'*.png', 'Png files (*.png)'; ...
        '*.jpg;*.jpeg', 'JPEG Files (*.jpg;*.jpeg)'; ...
        '*.jp2;*.jpf;*.jpx;*.j2c;*.j2k', 'JPEG-2000 Files (*.jp2;*.jpf;*.jpx;*.j2c;*.j2k)'; ...
        '*.tif;*.tiff', 'TIFF Files (*.tiff;*.tif)'; ...
        '*.*',  'All Files (*.*)'};
      obj.img_fn_save_fh = [];
      obj.img_fn_open_fh = [];
      obj.save_fh = [];
      obj.open_fh = [];
      obj.save_mat_fh = [];
      
      obj.create_ui();
      set_img(obj, img);
    end
    
    function delete(obj)
      try
        delete(obj.h_gui.min_slider);
      end
      try
        delete(obj.h_gui.max_slider);
      end
      try
        delete(obj.h_gui.hist_slider);
      end
      try
        delete(obj.h_fig);
      end
    end
    
    create_ui(obj);

    function set_auto_caxis(obj,auto_caxis_setting)
      set(obj.h_gui.caxis_autoCB,'Value',auto_caxis_setting);
      if auto_caxis_setting
        % If switching to auto setting, then update the caxis settings
        % since this may cause the caxis limits to change
        obj.update_caxis();
      end
    end

    function toggle_visibility(obj,visibility_flag)
      if exist('visibility_flag','var')
        cur_visibility = ~visibility_flag;
      else
        cur_visibility = strcmpi(get(obj.h_fig,'Visible'),'on')
      end
      if cur_visibility
        set(obj.h_fig,'Visible','off');
      else
        set(obj.h_fig,'Visible','on');
        figure(obj.h_fig);
      end
    end
    
    function set_img(obj, img)
      if obj.img == -1
        reset_vals = true;
      else
        reset_vals = false;
      end
      obj.img = img;
      if img == -1
        return;
      end
      obj.CData = get(obj.img,'CData');
      obj.update_filter();
      obj.set_limits(reset_vals);
      obj.update_caxis();
    end
    
    function set_cdata(obj, CData)
      obj.CData = CData;
      obj.update_filter();
      obj.set_limits(false);
      obj.update_caxis();
      
      if ~ishandle(obj.img)
        return;
      end
      C = get(obj.img,'CData');
      if isempty(C)
        return;
      end
      clims = [finitemin(C(:)) finitemax(C(:))];
      
      if numel(clims) == 2 && all(isfinite(clims))
        if get(obj.h_gui.caxis_autoCB,'Value')
          caxis(get(obj.img,'parent'),clims);
          obj.h_gui.min_slider.set_value(clims(1));
          obj.h_gui.max_slider.set_value(clims(end));
        else
          % Increase clims if necessary and keep slider values where they
          % are at.
          clims(1) = min(clims(1),obj.h_gui.min_slider.get_value());
          clims(2) = max(clims(2),obj.h_gui.max_slider.get_value());
        end
        
        obj.h_gui.min_slider.set_value_range(clims);
        obj.h_gui.max_slider.set_value_range(clims);
        
        cmap = get(obj.h_gui.colormapPM,'Value');
        if cmap == 1
          base_cmap = 1-gray(256);
        else
          base_cmap = jet(256);
        end
        hist_exp = obj.h_gui.hist_slider.get_value();
        cmap = interp1(linspace(0,1,256).', base_cmap, linspace(0,1,256).^(10.^(hist_exp/10)));
        colormap(get(obj.img,'parent'),cmap);
      end
    end
   
    function clims = get_limits(obj)
      clims = [0 1];
      if ishandle(obj.img)
        CData = get(obj.img,'CData');
        min_val = min(CData(isfinite(CData(:))));
        if ~isempty(min_val)
          clims(1) = min_val;
          clims(2) = max(CData(isfinite(CData(:))));
        end
      end
    end
    
    function set_limits(obj,reset_vals)
      clims = obj.get_limits();
      cur_min = obj.h_gui.min_slider.get_min();
      cur_max = obj.h_gui.min_slider.get_max();
      if reset_vals
        cur_min = inf;
        cur_max = -inf;
      end
      if clims(1) > cur_min
        clims(1) = cur_min;
      elseif clims(1) < cur_min
        obj.h_gui.min_slider.set_min(clims(1));
        obj.h_gui.max_slider.set_min(clims(1));
      end
      if clims(2) < cur_max
        clims(2) = cur_max;
      elseif clims(2) > cur_max
        obj.h_gui.min_slider.set_max(clims(2));
        obj.h_gui.max_slider.set_max(clims(2));
      end
      if reset_vals
        obj.h_gui.min_slider.set_value(clims(1));
        obj.h_gui.max_slider.set_value(clims(end));
      end
    end
    
    close_win(obj,varargin);
    
    % Call this update_caxis in your zoom/pan functions so that auto-caxis works
    update_caxis(obj,varargin);
    
    % Update for filter selection changes
    update_filter(obj,varargin);

    % Save pushbutton callback
    save_callback(obj,h_obj,event);
    save_mat_callback(obj,h_obj,event);

    % Open pushbutton callback
    open_callback(obj,h_obj,event);
    
  end
  
end
