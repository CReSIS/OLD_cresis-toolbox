% Class vector_editor
%
% Class which creates a GUI for editing geographic vector data.
% Useful for flight or mission planning.
%
% obj = vector_editor(ct_filename_gis([],'antarctica\Landsat-7\Antarctica_LIMA_480m.tif'));
% obj = vector_editor(ct_filename_gis([],'greenland\Landsat-7\mzl7geo_90m_lzw.tif'));

classdef (HandleCompatible = true) vector_editor < handle
  properties
    flines % Struct array of flight line vector data
    plotonly % graphics_handler class
    cur % Index to selected flight line
    geotiff_fn % geotiff filename
    proj % Projection information (empty is geodetic)
    
    selection % Structure of selection information
    % selection.button_down coordinates
    % selection.h_wnt_plot
    
    h_fig % Figure handle
    h_axes % Axes handle
    h_image % Image handle
    h_gui % Structure of graphics handles
    zoom_mode % logical indicating the current zoom mode
    zoom_mode_x % double containing x-position of button down
    zoom_mode_y % double containing y-position of button down
    
    save_fn % Filename to save lines to
    export_fn % Filename to export to
    open_fn_dir % Path to open files from
    special_mode % Structure for special selection mode
  end
  
  methods
    function obj = vector_editor(geotiff_fn,save_fn)
      if ~exist('geotiff_fn','var')
        geotiff_fn = '';
      end
      obj.geotiff_fn = geotiff_fn;
      if ~exist('save_fn','var')
        save_fn = '';
      end
      obj.save_fn = save_fn;
      
      obj.flines = struct([]);
      
      save_fn_dir = fileparts(obj.save_fn);
      obj.save_fn = obj.save_fn;
      obj.export_fn = fullfile(save_fn_dir,'export');
      obj.open_fn_dir = save_fn_dir;
      
      % geotiff_fn = Geotiff to load, leave empty for no raster image
      obj.special_mode.type = false;
      obj.h_fig = figure;
      h_fig_pos = get(obj.h_fig,'Position');
      set(obj.h_fig,'Position',h_fig_pos + [0 -150 150 150]);
      set(obj.h_fig,'NumberTitle','off');
      if isa(obj.h_fig,'double')
          set(obj.h_fig,'Name',sprintf('Vec %d',obj.h_fig));
      else
          set(obj.h_fig,'Name',sprintf('Vec %d',obj.h_fig.Number));
      end
      set(obj.h_fig,'MenuBar','none');
      set(obj.h_fig,'ToolBar','none');
      
      % h_gui.h_panel
      % h_gui.h_axes
      % h_gui.h_table
      % h_gui.h_panel_table
      % h_gui.h_flines
      %% Create widgets of main table
      obj.h_gui.h_lpanel = uipanel('Parent',obj.h_fig);
      set(obj.h_gui.h_lpanel,'Title','Vector Editor');
      set(obj.h_gui.h_lpanel,'TitlePosition','CenterTop');
      set(obj.h_gui.h_lpanel,'HighlightColor',[0.8 0.8 0.8]);
      set(obj.h_gui.h_lpanel,'ShadowColor',[0.6 0.6 0.6]);
      obj.h_gui.h_rpanel = uipanel('Parent',obj.h_fig);
      set(obj.h_gui.h_rpanel,'HighlightColor',[0.8 0.8 0.8]);
      set(obj.h_gui.h_rpanel,'ShadowColor',[0.6 0.6 0.6]);
      
      obj.h_axes = axes('Parent',obj.h_gui.h_rpanel);
      title(obj.save_fn,'interpreter','none','parent',obj.h_axes);
      obj.plotonly = graphics_handler(obj.h_axes,'none');
      set(obj.plotonly.h_fig,'NumberTitle','off');
      if isa(obj.h_fig,'double')
          set(obj.plotonly.h_fig,'Name',sprintf('VecGH %d %d',obj.h_fig,obj.plotonly.h_fig));
      else
          set(obj.plotonly.h_fig,'Name',sprintf('VecGH %d %d',obj.h_fig.Number,obj.plotonly.h_fig.Number));
      end
      axes(obj.h_axes);
      
      obj.selection.h_wpnt_plot = plot(xlim,ylim,'rx','Parent',obj.h_axes,'LineWidth',3,'MarkerSize',15);
      set(obj.selection.h_wpnt_plot,'XData',[]);
      set(obj.selection.h_wpnt_plot,'YData',[]);
      
      obj.h_gui.statusText = uicontrol('parent',obj.h_gui.h_rpanel);
      set(obj.h_gui.statusText,'Style','text');
      set(obj.h_gui.statusText,'HorizontalAlignment','left');
      set(obj.h_gui.statusText,'String','Creating user interface');
      set(obj.h_gui.statusText,'Position',[1 1 120 20]);
      set(obj.h_gui.statusText,'Units','points');
      
      %% Setup main table
      obj.h_gui.h_table.ui = obj.h_fig;
      obj.h_gui.h_table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
      obj.h_gui.h_table.height_margin = NaN*zeros(30,30);
      obj.h_gui.h_table.false_height = NaN*zeros(30,30);
      row = 1; col = 1;
      obj.h_gui.h_table.handles{row,col}   = obj.h_gui.h_lpanel;
      obj.h_gui.h_table.width(row,col)     = 180;
      obj.h_gui.h_table.height(row,col)    = inf;
      obj.h_gui.h_table.false_height(row,col) = 0;
      col = col + 1;
      obj.h_gui.h_table.handles{row,col}   = obj.h_gui.h_rpanel;
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
      
      %% Setup right panel table
      obj.h_gui.h_rpanel_table.ui = obj.h_gui.h_rpanel;
      obj.h_gui.h_rpanel_table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
      obj.h_gui.h_rpanel_table.height_margin = NaN*zeros(30,30);
      obj.h_gui.h_rpanel_table.false_height = NaN*zeros(30,30);
      row = 1; col = 1;
      obj.h_gui.h_rpanel_table.handles{row,col}   = obj.h_axes;
      obj.h_gui.h_rpanel_table.width(row,col)     = inf;
      obj.h_gui.h_rpanel_table.height(row,col)    = inf;
      obj.h_gui.h_rpanel_table.width_margin(row,col) = 30;
      obj.h_gui.h_rpanel_table.height_margin(row,col) = 30;
      obj.h_gui.h_rpanel_table.false_height(row,col) = 0;
      row = row + 1; col = 1;
      obj.h_gui.h_rpanel_table.handles{row,col}   = obj.h_gui.statusText;
      obj.h_gui.h_rpanel_table.width(row,col)     = inf;
      obj.h_gui.h_rpanel_table.height(row,col)    = 20;
      obj.h_gui.h_rpanel_table.false_height(row,col) = 0;
      obj.h_gui.h_rpanel_table.width_margin(row,col) = 0;
      obj.h_gui.h_rpanel_table.height_margin(row,col) = 0;
      obj.h_gui.h_rpanel_table.width_margin ...
        = obj.h_gui.h_rpanel_table.width_margin(1:row,1:col);
      obj.h_gui.h_rpanel_table.height_margin ...
        = obj.h_gui.h_rpanel_table.height_margin(1:row,1:col);
      obj.h_gui.h_rpanel_table.false_height ...
        = obj.h_gui.h_rpanel_table.false_height(1:row,1:col);
      clear row col
      table_draw(obj.h_gui.h_rpanel_table);
      
      %% Setup left subpanel GUI objects
      obj.h_gui.flines.moveupPB = uicontrol('Parent',obj.h_gui.h_lpanel);
      set(obj.h_gui.flines.moveupPB,'Style','PushButton');
      set(obj.h_gui.flines.moveupPB,'String','^');
      set(obj.h_gui.flines.moveupPB,'Callback',@obj.moveupPB_callback);
      
      obj.h_gui.flines.movedownPB = uicontrol('Parent',obj.h_gui.h_lpanel);
      set(obj.h_gui.flines.movedownPB,'Style','PushButton');
      set(obj.h_gui.flines.movedownPB,'String','v');
      set(obj.h_gui.flines.movedownPB,'Callback',@obj.movedownPB_callback);
      
      obj.h_gui.flines.openPB = uicontrol('Parent',obj.h_gui.h_lpanel);
      set(obj.h_gui.flines.openPB,'Style','PushButton');
      set(obj.h_gui.flines.openPB,'String','Open');
      set(obj.h_gui.flines.openPB,'Callback',@obj.openPB_callback);
      
      obj.h_gui.flines.newPB = uicontrol('Parent',obj.h_gui.h_lpanel);
      set(obj.h_gui.flines.newPB,'Style','PushButton');
      set(obj.h_gui.flines.newPB,'String','New');
      set(obj.h_gui.flines.newPB,'Callback',@obj.newPB_callback);
      
      obj.h_gui.flines.openPB = uicontrol('Parent',obj.h_gui.h_lpanel);
      set(obj.h_gui.flines.openPB,'Style','PushButton');
      set(obj.h_gui.flines.openPB,'String','Open');
      set(obj.h_gui.flines.openPB,'Callback',@obj.openPB_callback);
      
      obj.h_gui.flines.savePB = uicontrol('Parent',obj.h_gui.h_lpanel);
      set(obj.h_gui.flines.savePB,'Style','PushButton');
      set(obj.h_gui.flines.savePB,'String','Save');
      set(obj.h_gui.flines.savePB,'Callback',@obj.savePB_callback);
      
      obj.h_gui.flines.saveasPB = uicontrol('Parent',obj.h_gui.h_lpanel);
      set(obj.h_gui.flines.saveasPB,'Style','PushButton');
      set(obj.h_gui.flines.saveasPB,'String','Save As');
      set(obj.h_gui.flines.saveasPB,'Callback',@obj.saveasPB_callback);
      
      obj.h_gui.toolPM = uicontrol('Parent',obj.h_gui.h_lpanel);
      set(obj.h_gui.toolPM,'Style','PopupMenu');
      menu_string{1} = 'Select';
      menu_string{2} = 'Insert';
      set(obj.h_gui.toolPM,'String',menu_string);
      set(obj.h_gui.toolPM,'Value',1);
      set(obj.h_gui.toolPM,'HorizontalAlignment','Center');
      set(obj.h_gui.toolPM,'Value',1);
      set(obj.h_gui.toolPM,'Callback',@obj.toolPM_callback);
      
      obj.h_gui.exportPB = uicontrol('Parent',obj.h_gui.h_lpanel);
      set(obj.h_gui.exportPB,'Style','PushButton');
      set(obj.h_gui.exportPB,'String','Export');
      set(obj.h_gui.exportPB,'Callback',@obj.exportPB_callback);
      
      obj.h_gui.wpnts_or_flinesCB = uicontrol('Parent',obj.h_gui.h_lpanel);
      set(obj.h_gui.wpnts_or_flinesCB,'Style','checkbox');
      set(obj.h_gui.wpnts_or_flinesCB,'HorizontalAlignment','Center');
      set(obj.h_gui.wpnts_or_flinesCB,'FontName','fixed');
      set(obj.h_gui.wpnts_or_flinesCB,'String','Line-Select');
      set(obj.h_gui.wpnts_or_flinesCB,'Value',0);
      
      obj.h_gui.multi_or_singleCB = uicontrol('Parent',obj.h_gui.h_lpanel);
      set(obj.h_gui.multi_or_singleCB,'Style','checkbox');
      set(obj.h_gui.multi_or_singleCB,'HorizontalAlignment','Center');
      set(obj.h_gui.multi_or_singleCB,'FontName','fixed');
      set(obj.h_gui.multi_or_singleCB,'String','Multi-Select');
      set(obj.h_gui.multi_or_singleCB,'Value',0);
      
      obj.h_gui.zoom_selectCB = uicontrol('Parent',obj.h_gui.h_lpanel);
      set(obj.h_gui.zoom_selectCB,'Style','checkbox');
      set(obj.h_gui.zoom_selectCB,'HorizontalAlignment','Center');
      set(obj.h_gui.zoom_selectCB,'FontName','fixed');
      set(obj.h_gui.zoom_selectCB,'String','Zoom-Select');
      set(obj.h_gui.zoom_selectCB,'Value',1);
      
      obj.h_gui.wpnts.insert_beforePB = uicontrol('Parent',obj.h_gui.h_lpanel);
      set(obj.h_gui.wpnts.insert_beforePB,'Style','PushButton');
      set(obj.h_gui.wpnts.insert_beforePB,'String','Insert-Before');
      set(obj.h_gui.wpnts.insert_beforePB,'Callback',@obj.insert_beforePB_callback);
      
      obj.h_gui.wpnts.insert_afterPB = uicontrol('Parent',obj.h_gui.h_lpanel);
      set(obj.h_gui.wpnts.insert_afterPB,'Style','PushButton');
      set(obj.h_gui.wpnts.insert_afterPB,'String','Insert-After');
      set(obj.h_gui.wpnts.insert_afterPB,'Callback',@obj.insert_afterPB_callback);
      
      obj.h_gui.wpnts.moveupPB = uicontrol('Parent',obj.h_gui.h_lpanel);
      set(obj.h_gui.wpnts.moveupPB,'Style','PushButton');
      set(obj.h_gui.wpnts.moveupPB,'String','^');
      set(obj.h_gui.wpnts.moveupPB,'Callback',@obj.moveupPB_callback);
      
      obj.h_gui.wpnts.movedownPB = uicontrol('Parent',obj.h_gui.h_lpanel);
      set(obj.h_gui.wpnts.movedownPB,'Style','PushButton');
      set(obj.h_gui.wpnts.movedownPB,'String','v');
      set(obj.h_gui.wpnts.movedownPB,'Callback',@obj.movedownPB_callback);
      
      obj.h_gui.geo_sortPB = uicontrol('Parent',obj.h_gui.h_lpanel);
      set(obj.h_gui.geo_sortPB,'Style','PushButton');
      set(obj.h_gui.geo_sortPB,'String','Geometry Sort');
      set(obj.h_gui.geo_sortPB,'Callback',@obj.geo_sortPB_callback);
      
      obj.h_gui.name_sortPB = uicontrol('Parent',obj.h_gui.h_lpanel);
      set(obj.h_gui.name_sortPB,'Style','PushButton');
      set(obj.h_gui.name_sortPB,'String','Name Sort');
      set(obj.h_gui.name_sortPB,'Callback',@obj.name_sortPB_callback);
      
      %% Setup left subpanel table
      obj.h_gui.h_subpanel_table.ui = [];
      obj.h_gui.h_subpanel_table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
      obj.h_gui.h_subpanel_table.height_margin = NaN*zeros(30,30);
      
      row = 1; col = 1;
      obj.h_gui.h_subpanel_table.handles{row,col}   = obj.h_gui.flines.moveupPB;
      obj.h_gui.h_subpanel_table.width(row,col)     = inf;
      obj.h_gui.h_subpanel_table.height(row,col)    = 25;
      col = col + 1;
      obj.h_gui.h_subpanel_table.handles{row,col}   = obj.h_gui.flines.movedownPB;
      obj.h_gui.h_subpanel_table.width(row,col)     = inf;
      obj.h_gui.h_subpanel_table.height(row,col)    = 25;
      col = 1; row = row + 1;
      obj.h_gui.h_subpanel_table.handles{row,col}   = obj.h_gui.flines.newPB;
      obj.h_gui.h_subpanel_table.width(row,col)     = inf;
      obj.h_gui.h_subpanel_table.height(row,col)    = 25;
      col = col + 1;
      obj.h_gui.h_subpanel_table.handles{row,col}   = obj.h_gui.flines.openPB;
      obj.h_gui.h_subpanel_table.width(row,col)     = inf;
      obj.h_gui.h_subpanel_table.height(row,col)    = 25;
      col = 1; row = row + 1;
      obj.h_gui.h_subpanel_table.handles{row,col}   = obj.h_gui.flines.savePB;
      obj.h_gui.h_subpanel_table.width(row,col)     = inf;
      obj.h_gui.h_subpanel_table.height(row,col)    = 25;
      col = col + 1;
      obj.h_gui.h_subpanel_table.handles{row,col}   = obj.h_gui.flines.saveasPB;
      obj.h_gui.h_subpanel_table.width(row,col)     = inf;
      obj.h_gui.h_subpanel_table.height(row,col)    = 25;
      col = 1; row = row + 1;
      obj.h_gui.h_subpanel_table.handles{row,col}   = obj.h_gui.toolPM;
      obj.h_gui.h_subpanel_table.width(row,col)     = inf;
      obj.h_gui.h_subpanel_table.height(row,col)    = 25;
      col = col + 1;
      obj.h_gui.h_subpanel_table.handles{row,col}   = obj.h_gui.exportPB;
      obj.h_gui.h_subpanel_table.width(row,col)     = inf;
      obj.h_gui.h_subpanel_table.height(row,col)    = 25;
      col = 1; row = row + 1;
      obj.h_gui.h_subpanel_table.handles{row,col}   = obj.h_gui.wpnts_or_flinesCB;
      obj.h_gui.h_subpanel_table.width(row,col)     = inf;
      obj.h_gui.h_subpanel_table.height(row,col)    = 25;
      col = col + 1;
      obj.h_gui.h_subpanel_table.handles{row,col}   = obj.h_gui.multi_or_singleCB;
      obj.h_gui.h_subpanel_table.width(row,col)     = inf;
      obj.h_gui.h_subpanel_table.height(row,col)    = 25;
      col = col + 1;
      obj.h_gui.h_subpanel_table.handles{row,col}   = obj.h_gui.zoom_selectCB;
      obj.h_gui.h_subpanel_table.width(row,col)     = inf;
      obj.h_gui.h_subpanel_table.height(row,col)    = 25;
      col = 1; row = row + 1;
      obj.h_gui.h_subpanel_table.handles{row,col}   = obj.h_gui.wpnts.insert_beforePB;
      obj.h_gui.h_subpanel_table.width(row,col)     = inf;
      obj.h_gui.h_subpanel_table.height(row,col)    = 25;
      col = col + 1;
      obj.h_gui.h_subpanel_table.handles{row,col}   = obj.h_gui.wpnts.insert_afterPB;
      obj.h_gui.h_subpanel_table.width(row,col)     = inf;
      obj.h_gui.h_subpanel_table.height(row,col)    = 25;
      col = 1; row = row + 1;
      obj.h_gui.h_subpanel_table.handles{row,col}   = obj.h_gui.wpnts.moveupPB;
      obj.h_gui.h_subpanel_table.width(row,col)     = inf;
      obj.h_gui.h_subpanel_table.height(row,col)    = 25;
      col = col + 1;
      obj.h_gui.h_subpanel_table.handles{row,col}   = obj.h_gui.wpnts.movedownPB;
      obj.h_gui.h_subpanel_table.width(row,col)     = inf;
      obj.h_gui.h_subpanel_table.height(row,col)    = 25;
      col = 1; row = row + 1;
      obj.h_gui.h_subpanel_table.handles{row,col}   = obj.h_gui.geo_sortPB;
      obj.h_gui.h_subpanel_table.width(row,col)     = inf;
      obj.h_gui.h_subpanel_table.height(row,col)    = 25;
      col = col + 1;
      obj.h_gui.h_subpanel_table.handles{row,col}   = obj.h_gui.name_sortPB;
      obj.h_gui.h_subpanel_table.width(row,col)     = inf;
      obj.h_gui.h_subpanel_table.height(row,col)    = 25;
      
      obj.h_gui.h_subpanel_table.width_margin ...
        = obj.h_gui.h_subpanel_table.width_margin(1:row,1:col);
      obj.h_gui.h_subpanel_table.height_margin ...
        = obj.h_gui.h_subpanel_table.height_margin(1:row,1:col);
      clear row col
      
      %% Setup left panel GUI Objects
      obj.h_gui.flines.listLB = uicontrol('Parent',obj.h_gui.h_lpanel);
      set(obj.h_gui.flines.listLB,'Style','listbox');
      set(obj.h_gui.flines.listLB,'HorizontalAlignment','Center');
      set(obj.h_gui.flines.listLB,'FontName','fixed');
      set(obj.h_gui.flines.listLB,'Value',[]);
      set(obj.h_gui.flines.listLB,'Callback',@obj.flinesLB_callback);
      set(obj.h_gui.flines.listLB,'Max',1e9);
      
      %% Source list box context menu
      obj.h_gui.flines.listCM = uicontextmenu;
      % Define the context menu items and install their callbacks
      obj.h_gui.flines.listCM_item1 = uimenu(obj.h_gui.flines.listCM, 'Label', 'Copy', 'Callback', @obj.flinesLB_menu_callback);
      obj.h_gui.flines.listCM_item1 = uimenu(obj.h_gui.flines.listCM, 'Label', 'Delete', 'Callback', @obj.flinesLB_menu_callback);
      obj.h_gui.flines.listCM_item1 = uimenu(obj.h_gui.flines.listCM, 'Label', 'Geometry Sort', 'Callback', @obj.flinesLB_menu_callback);
      obj.h_gui.flines.listCM_item1 = uimenu(obj.h_gui.flines.listCM, 'Label', 'Merge', 'Callback', @obj.flinesLB_menu_callback);
      obj.h_gui.flines.listCM_item1 = uimenu(obj.h_gui.flines.listCM, 'Label', 'Move', 'Callback', @obj.flinesLB_menu_callback);
      obj.h_gui.flines.listCM_item1 = uimenu(obj.h_gui.flines.listCM, 'Label', 'Paste', 'Callback', @obj.flinesLB_menu_callback);
      obj.h_gui.flines.listCM_item1 = uimenu(obj.h_gui.flines.listCM, 'Label', 'Paste Special', 'Callback', @obj.flinesLB_menu_callback);
      obj.h_gui.flines.listCM_item1 = uimenu(obj.h_gui.flines.listCM, 'Label', 'Rename', 'Callback', @obj.flinesLB_menu_callback);
      obj.h_gui.flines.listCM_item1 = uimenu(obj.h_gui.flines.listCM, 'Label', 'Reverse Flines', 'Callback', @obj.flinesLB_menu_callback);
      obj.h_gui.flines.listCM_item1 = uimenu(obj.h_gui.flines.listCM, 'Label', 'Reverse Wpnts', 'Callback', @obj.flinesLB_menu_callback);
      obj.h_gui.flines.listCM_item1 = uimenu(obj.h_gui.flines.listCM, 'Label', 'Rotate', 'Callback', @obj.flinesLB_menu_callback);
      set(obj.h_gui.flines.listLB,'uicontextmenu',obj.h_gui.flines.listCM)
      
      obj.h_gui.wpnts.listLB = uicontrol('Parent',obj.h_gui.h_lpanel);
      set(obj.h_gui.wpnts.listLB,'Style','listbox');
      set(obj.h_gui.wpnts.listLB,'HorizontalAlignment','Center');
      set(obj.h_gui.wpnts.listLB,'FontName','fixed');
      set(obj.h_gui.wpnts.listLB,'Value',[]);
      set(obj.h_gui.wpnts.listLB,'Callback',@obj.wpntsLB_callback);
      set(obj.h_gui.wpnts.listLB,'Max',1e9);
      
      %% Source list box context menu
      obj.h_gui.wpnts.listCM = uicontextmenu;
      % Define the context menu items and install their callbacks
      obj.h_gui.wpnts.listCM_item1 = uimenu(obj.h_gui.wpnts.listCM, 'Label', 'Change Length', 'Callback', @obj.wpntsLB_menu_callback);
      obj.h_gui.wpnts.listCM_item1 = uimenu(obj.h_gui.wpnts.listCM, 'Label', 'Edit Geographic', 'Callback', @obj.wpntsLB_menu_callback);
      obj.h_gui.wpnts.listCM_item1 = uimenu(obj.h_gui.wpnts.listCM, 'Label', 'Edit Map', 'Callback', @obj.wpntsLB_menu_callback);
      obj.h_gui.wpnts.listCM_item1 = uimenu(obj.h_gui.wpnts.listCM, 'Label', 'Insert Geographic', 'Callback', @obj.wpntsLB_menu_callback);
      obj.h_gui.wpnts.listCM_item1 = uimenu(obj.h_gui.wpnts.listCM, 'Label', 'Insert Map', 'Callback', @obj.wpntsLB_menu_callback);
      obj.h_gui.wpnts.listCM_item1 = uimenu(obj.h_gui.wpnts.listCM, 'Label', 'Copy', 'Callback', @obj.wpntsLB_menu_callback);
      obj.h_gui.wpnts.listCM_item1 = uimenu(obj.h_gui.wpnts.listCM, 'Label', 'Delete', 'Callback', @obj.wpntsLB_menu_callback);
      obj.h_gui.wpnts.listCM_item1 = uimenu(obj.h_gui.wpnts.listCM, 'Label', 'Paste', 'Callback', @obj.wpntsLB_menu_callback);
      obj.h_gui.wpnts.listCM_item1 = uimenu(obj.h_gui.wpnts.listCM, 'Label', 'Rename', 'Callback', @obj.wpntsLB_menu_callback);
      obj.h_gui.wpnts.listCM_item1 = uimenu(obj.h_gui.wpnts.listCM, 'Label', 'Rename Auto', 'Callback', @obj.wpntsLB_menu_callback);
      obj.h_gui.wpnts.listCM_item1 = uimenu(obj.h_gui.wpnts.listCM, 'Label', 'Reverse Order', 'Callback', @obj.wpntsLB_menu_callback);
      obj.h_gui.wpnts.listCM_item1 = uimenu(obj.h_gui.wpnts.listCM, 'Label', 'Rotate', 'Callback', @obj.wpntsLB_menu_callback);
      set(obj.h_gui.wpnts.listLB,'uicontextmenu',obj.h_gui.wpnts.listCM)
      
      %% Setup left panel
      obj.h_gui.h_lpanel_table.ui = obj.h_gui.h_lpanel;
      obj.h_gui.h_lpanel_table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
      obj.h_gui.h_lpanel_table.height_margin = NaN*zeros(30,30);
      obj.h_gui.h_lpanel_table.false_width = NaN*zeros(30,30); % Just make these bigger than they have to be
      obj.h_gui.h_lpanel_table.offset = [0 10];
      row = 1; col = 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.flines.listLB;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = inf;
      obj.h_gui.h_lpanel_table.width_margin(row,col) = 3;
      obj.h_gui.h_lpanel_table.height_margin(row,col) = 3;
      obj.h_gui.h_lpanel_table.false_width(row,col) = 3;
      col = 1; row = row + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.h_subpanel_table;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 8*25;
      obj.h_gui.h_lpanel_table.width_margin(row,col) = 0;
      obj.h_gui.h_lpanel_table.height_margin(row,col) = 0;
      obj.h_gui.h_lpanel_table.false_width(row,col) = 3;
      col = 1; row = row + 1;
      obj.h_gui.h_lpanel_table.handles{row,col}   = obj.h_gui.wpnts.listLB;
      obj.h_gui.h_lpanel_table.width(row,col)     = inf;
      obj.h_gui.h_lpanel_table.height(row,col)    = 75;
      obj.h_gui.h_lpanel_table.width_margin(row,col) = 3;
      obj.h_gui.h_lpanel_table.height_margin(row,col) = 3;
      obj.h_gui.h_lpanel_table.false_width(row,col) = 3;
      
      obj.h_gui.h_lpanel_table.width_margin ...
        = obj.h_gui.h_lpanel_table.width_margin(1:row,1:col);
      obj.h_gui.h_lpanel_table.height_margin ...
        = obj.h_gui.h_lpanel_table.height_margin(1:row,1:col);
      obj.h_gui.h_lpanel_table.false_width ...
        = obj.h_gui.h_lpanel_table.false_width(1:row,1:col);
      clear row col
      table_draw(obj.h_gui.h_lpanel_table);
      
      %% Set up general handles
      set(obj.h_fig,'WindowButtonDownFcn',@obj.button_down);
      set(obj.h_fig,'WindowButtonUpFcn',@obj.button_up);
      set(obj.h_fig,'WindowScrollWheelFcn',@obj.button_scroll);
      set(obj.h_fig,'CloseRequestFcn',@obj.close_win);
      set(obj.h_fig,'WindowKeyPressFcn',@obj.key_press);
      
      % Set up zoom
      zoom_setup(obj.h_fig);
      obj.zoom_mode = true;
      set(obj.h_fig,'pointer','custom');
      
      %% Load file if passed in to constructor
      [save_fn_dir,save_fn_name,save_fn_ext] = fileparts(obj.save_fn);
      if strcmpi(save_fn_ext,'.mat')
        obj.openMatFile(obj.save_fn,0,1,1,1);
      else
        obj.update_geotiff(false);
      end
      
    end
    
    function delete(obj)
      % Delete the map figure handle
      try
        delete(obj.h_fig);
      end
      try
         delete(obj.plotonly);
      end
    end
    
    function close_win(obj,h_obj,event)
      delete(obj);
    end
    
    function update_geotiff(obj,update_graphics)
      if ~exist(obj.geotiff_fn,'file')
        warning('Valid geotiff file required.\n  Current file: (%s).', obj.geotiff_fn);
        
        [geotiff_fn, geotiff_fn_dir, filterindex] = uigetfile( ...
          {'*.tif','GeoTiff files (*.tif)'}, ...
          'Select Geotiff map file', ...
          obj.geotiff_fn, 'MultiSelect', 'off');
        
        if isequal(geotiff_fn,0) || isequal(geotiff_fn_dir,0)
          if isempty(obj.proj)
            % No projection, cannot continue
            error('Aborting... failed to open (%s)',geotiff_fn);
          else
            warning('Aborting... failed to open (%s)',geotiff_fn);
            return;
          end
        end
        
        obj.geotiff_fn = fullfile(geotiff_fn_dir, geotiff_fn);
      end
      
      obj.proj = geotiffinfo(obj.geotiff_fn);
      
      % Read the image
      fprintf('Reading the geotiff %s... may take a while\n', obj.geotiff_fn);
      [RGB, CMAP, R, tmp] = geotiffread(obj.geotiff_fn);
      fprintf('  Done loading geotiff\n');
      if size(RGB,3) == 3 && strcmp(class(RGB),'uint16') && max(RGB(:)) <= 255
        RGB = uint8(RGB);
      end
      R = R/1e3;
      
      if isa(RGB,'int16')
        RGB = double(RGB);
        RGB(RGB == 32767 | RGB == -32767 | RGB == -32768) = NaN;
        RGB = (RGB - min(RGB(:))) / (max(RGB(:)) - min(RGB(:)));
      elseif isa(RGB,'single')
        RGB = double(RGB);
      end
      
      % Store all the existing plotonly objects
      plotonly = handle2struct(obj.plotonly.handles);
      obj.plotonly.delete_handle(1:length(obj.plotonly.handles));
      
      if isobject(obj.h_image)
        delete(obj.h_image);
      end
      obj.h_image = mapshow(RGB,R,'Parent',obj.h_axes);
      xlabel('X (km)');
      ylabel('Y (km)');
      
      % Put all the existing plotonly objects back on the plot
      new_plotonly = struct2handle(plotonly,obj.h_axes);
      obj.plotonly.insert_handle(new_plotonly);
      
      % For each plotonly handle, update the projection
      for plotonly_idx=1:length(obj.plotonly.handles)
        userdata = get(obj.plotonly.handles(plotonly_idx),'userdata');
        [x,y] = projfwd(obj.proj,userdata.lat,userdata.lon);
        set(obj.plotonly.handles(plotonly_idx),'XData',x/1e3,'YData',y/1e3);
      end
      
      for pos = 1:length(obj.flines)
        [obj.flines(pos).x,obj.flines(pos).y] ...
          = projfwd(obj.proj,obj.flines(pos).lat,obj.flines(pos).lon);
      end
      
      if update_graphics
        obj.update_statusText();
        obj.update_flineGraphics();
      end
    end
    
    function button_up(obj,h_obj,event)
      [x,y,but] = get_mouse_info(obj.h_fig,obj.h_axes);
      
      % Make sure that click is on the right side panel
      mouse_pos = get(obj.h_fig,'CurrentPoint');
      
      % Check to make sure mouse clicked inside of obj.h_axes
      %   Since extends the full y-length, just check to the right of minimum x
      set(obj.h_gui.h_rpanel,'Units','normalized');
      uipanel_pos = get(obj.h_gui.h_rpanel,'Position');
      set(obj.h_gui.h_rpanel,'Units','Points');
      if mouse_pos(1) <= uipanel_pos(1)
        return
      end
      
      if obj.zoom_mode && obj.special_mode.type ~= 2
        xlims = get(obj.h_image,'XData'); xlims = sort(xlims([1 end]));
        ylims = get(obj.h_image,'YData'); ylims = sort(ylims([1 end]));
        zoom_button_up(x,y,but,struct('x',obj.zoom_mode_x,'y',obj.zoom_mode_y, ...
          'h_axes',obj.h_axes,'xlims',xlims,'ylims',ylims));
        return;
      end
      
      tools = get(obj.h_gui.toolPM,'String');
      tool = tools{get(obj.h_gui.toolPM,'Value')};
      if but == 1
        if obj.special_mode.type == 2
          return;
        end
        if strcmpi(tool,'select')
          fline_select = get(obj.h_gui.wpnts_or_flinesCB,'Value');
          multi_select = get(obj.h_gui.multi_or_singleCB,'Value');
          
          if multi_select
            %           obj.selection.button_down.x = x;
            %           obj.selection.button_down.y = y;
            %   xPoly = polyPts(:,1);
            %   yPoly = polyPts(:,2);
            %           [xF,yF] = projfwd(proj,LAT,LON);
            %           xF = xF/1e3;   yF = yF/1e3;
            %           IN = inpolygon(xF,yF,xPoly,yPoly);
            
          elseif fline_select
            % Find the closest flight line
            % interpm, reducem
            xlims = xlim(obj.h_axes);
            ylims = ylim(obj.h_axes);
            if x<xlims(1) || x>xlims(2) || y<ylims(1) || y>ylims(2)
              return
            end
            x = x*1e3;
            y = y*1e3;
            min_dist = inf;
            min_pos = [];
            for pos = 1:length(obj.flines)
              [lat_dense,lon_dense] = interpm(obj.flines(pos).lat,obj.flines(pos).lon,3e-2);
              [x_dense,y_dense] = projfwd(obj.proj,lat_dense,lon_dense);
              dist = min((x_dense - x).^2 + (y_dense - y).^2);
              if dist < min_dist
                min_dist = dist;
                min_pos = pos;
              end
            end
            if isempty(min_pos)
              warning('No way points loaded to select');
              return;
            end
            % Find the closest waypoint
            min_dist = inf;
            min_wpnt_idx = [];
            for wpnt_idx = 1:length(obj.flines(min_pos).x)
              dist = sqrt((obj.flines(min_pos).x(wpnt_idx) - x).^2 ...
                + (obj.flines(min_pos).y(wpnt_idx) - y).^2);
              if dist < min_dist
                min_dist = dist;
                min_wpnt_idx = wpnt_idx;
              end
            end
            set(obj.h_gui.flines.listLB,'Value',min_pos);
            %% Waypoint select
            set(obj.h_gui.wpnts.listLB,'Value',min_wpnt_idx);
            flinesLB_callback(obj,h_obj,event);
            
          else
            % Find the closest waypoint
            % interpm, reducem
            xlims = xlim;
            ylims = ylim;
            if x<xlims(1) || x>xlims(2) || y<ylims(1) || y>ylims(2)
              return
            end
            x = x*1e3;
            y = y*1e3;
            min_dist = inf;
            min_pos = [];
            for pos = 1:length(obj.flines)
              for wpnt_idx = 1:length(obj.flines(pos).x)
                dist = sqrt((obj.flines(pos).x(wpnt_idx) - x).^2 ...
                  + (obj.flines(pos).y(wpnt_idx) - y).^2);
                if dist < min_dist
                  min_dist = dist;
                  min_pos = pos;
                  min_wpnt_idx = wpnt_idx;
                end
              end
            end
            if isempty(min_pos)
              warning('No way points loaded to select');
              return;
            end
            set(obj.h_gui.flines.listLB,'Value',min_pos);
            %% Fill way points listbox with selection
            set(obj.h_gui.wpnts.listLB,'Value',[]);
            set(obj.h_gui.wpnts.listLB,'String', ...
              obj.flines(min_pos).wpnt_names);
            set(obj.h_gui.wpnts.listLB,'Value',min_wpnt_idx);
            wpntsLB_callback(obj,h_obj,event);
            obj.update_flineGraphics();
            
          end
        elseif strcmpi(tool,'insert')
            xlims = xlim(obj.h_axes);
            ylims = ylim(obj.h_axes);
            if x<xlims(1) || x>xlims(2) || y<ylims(1) || y>ylims(2)
              return
            end
          obj.insert_wpnt(x,y);
          obj.update_statusText();
          obj.wpntsLB_callback();
        end
      elseif but == 3
      end
    end
    
    function button_down(obj,h_obj,event)
      [x,y,but] = get_mouse_info(obj.h_fig,obj.h_axes);
      
      % Make sure that click is on the right side panel
      mouse_pos = get(obj.h_fig,'CurrentPoint');
      
      % Check to make sure mouse clicked inside of obj.h_axes
      %   Since extends the full y-length, just check to the right of minimum x
      set(obj.h_gui.h_rpanel,'Units','normalized');
      uipanel_pos = get(obj.h_gui.h_rpanel,'Position');
      set(obj.h_gui.h_rpanel,'Units','Points');
      if mouse_pos(1) <= uipanel_pos(1)
        return
      end
      
      if obj.zoom_mode && obj.special_mode.type ~= 2
        obj.zoom_mode_x = x;
        obj.zoom_mode_y = y;
        rbbox;
        return;
      end
      
      tools = get(obj.h_gui.toolPM,'String');
      tool = tools{get(obj.h_gui.toolPM,'Value')};
      if but == 1
        if strcmpi(tool,'select')
          multi_select = get(obj.h_gui.multi_or_singleCB,'Value');
          if multi_select
            obj.selection.button_down.x = x;
            obj.selection.button_down.y = y;
            rbbox;
          end
        end
      elseif but == 3
      end
    end
    
    function button_scroll(obj,h_obj,event)
      [x,y,but] = get_mouse_info(obj.h_fig,obj.h_axes);
      
      % Make sure that click is on the right side panel
      mouse_pos = get(obj.h_fig,'CurrentPoint');
      
      % Check to make sure mouse clicked inside of obj.h_axes
      %   Since extends the full y-length, just check to the right of minimum x
      set(obj.h_gui.h_rpanel,'Units','normalized');
      uipanel_pos = get(obj.h_gui.h_rpanel,'Position');
      set(obj.h_gui.h_rpanel,'Units','Points');
      if mouse_pos(1) <= uipanel_pos(1)
        return
      end
      
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
      
      xlim(obj.h_axes,xlims);
      ylim(obj.h_axes,ylims);
      
    end
    
    function key_press(obj,src,event)
      % see event.Modifier for modifiers

      current_object = gco;
      xlims = get(obj.h_image,'XData'); xlims = sort(xlims([1 end]));
      ylims = get(obj.h_image,'YData'); ylims = sort(ylims([1 end]));
      
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
          if any(current_object == obj.h_gui.flines.listLB) || any(current_object == obj.h_gui.wpnts.listLB)
            return
          end
          zoom_arrow(event,struct('h_axes',obj.h_axes, ...
            'xlims',xlims,'ylims',ylims));
          
        case 'uparrow' % Up-arrow: Pan up
          if any(current_object == obj.h_gui.flines.listLB) || any(current_object == obj.h_gui.wpnts.listLB)
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
    
    function flinesLB_callback(obj,h_obj,event)
      cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
      if iscell(cur_fline_selected)
        cur_fline_selected = cell2mat(cur_fline_selected);
      end
      if isempty(cur_fline_selected)
        % Empty way points listbox
        set(obj.h_gui.wpnts.listLB,'String',{});
        set(obj.h_gui.wpnts.listLB,'Value',[]);
        return;
      end
      
      %% Determine geographic limits of selected frames
      xlims = [inf -inf];
      ylims = [inf -inf];
      for pos = 1:length(obj.flines)
        if any(pos == cur_fline_selected)
          if min(obj.flines(pos).x) < xlims(1)
            xlims(1) = min(obj.flines(pos).x);
          end
          if max(obj.flines(pos).x) > xlims(2)
            xlims(2) = max(obj.flines(pos).x);
          end
          if min(obj.flines(pos).y) < ylims(1)
            ylims(1) = min(obj.flines(pos).y);
          end
          if max(obj.flines(pos).y) > ylims(2)
            ylims(2) = max(obj.flines(pos).y);
          end
        end
      end
      
      obj.update_flineGraphics();
      
      %% Set map axis
      if get(obj.h_gui.zoom_selectCB,'Value') && isfinite(xlims(1))
        axis(obj.h_axes,'normal');
        
        xlims = xlims/1e3;
        ylims = ylims/1e3;
        cur_xlim = xlim(obj.h_axes);
        cur_ylim = ylim(obj.h_axes);
        axes_position = get(obj.h_axes,'Position');
        aspect_ratio = axes_position(3) / axes_position(4);
        xlims(1) = min(xlims(1),cur_xlim(1));
        xlims(2) = max(xlims(2),cur_xlim(2));
        ylims(1) = min(ylims(1),cur_ylim(1));
        ylims(2) = max(ylims(2),cur_ylim(2));
        
        if diff(xlims) / diff(ylims) > aspect_ratio
          ylims = mean(ylims) + diff(xlims)/aspect_ratio/2 * [-1 1];
        else
          xlims = mean(xlims) + diff(ylims)*aspect_ratio/2 * [-1 1];
        end
        xlim(obj.h_axes,xlims);
        ylim(obj.h_axes,ylims);
      end
      
      cur_fline_selected = cur_fline_selected(1);
      
      %% Fill way points listbox with selection
      set(obj.h_gui.wpnts.listLB,'String', ...
        obj.flines(cur_fline_selected).wpnt_names);
      
      set(obj.h_gui.wpnts.listLB,'Value',[]);
      set(obj.h_gui.wpnts.listLB,'ListboxTop',1)
    end
    
    function wpntsLB_callback(obj,h_obj,event)
      cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
      if iscell(cur_fline_selected)
        cur_fline_selected = cell2mat(cur_fline_selected);
      end
      if isempty(cur_fline_selected)
        warning('Can not moveup: no flight line selected\n');
        return;
      end
      cur_fline_selected = cur_fline_selected(1);
      
      cur_wpnt_selected = get(obj.h_gui.wpnts.listLB,'Value');
      if iscell(cur_wpnt_selected)
        cur_wpnt_selected = cell2mat(cur_wpnt_selected);
      end
      
      set(obj.selection.h_wpnt_plot, ...
        'XData', obj.flines(cur_fline_selected).x(cur_wpnt_selected)/1e3, ...
        'YData', obj.flines(cur_fline_selected).y(cur_wpnt_selected)/1e3);
      
      %% Determine geographic limits of selected way points
      xlims = [inf -inf];
      ylims = [inf -inf];
      for wpnt = 1:length(obj.flines(cur_fline_selected).lat)
        if any(wpnt == cur_wpnt_selected)
          if obj.flines(cur_fline_selected).x(wpnt) < xlims(1)
            xlims(1) = obj.flines(cur_fline_selected).x(wpnt)-1;
          end
          if obj.flines(cur_fline_selected).x(wpnt) > xlims(2)
            xlims(2) = obj.flines(cur_fline_selected).x(wpnt)+1;
          end
          if obj.flines(cur_fline_selected).y(wpnt) < ylims(1)
            ylims(1) = obj.flines(cur_fline_selected).y(wpnt)-1;
          end
          if obj.flines(cur_fline_selected).y(wpnt) > ylims(2)
            ylims(2) = obj.flines(cur_fline_selected).y(wpnt)+1;
          end
        end
      end
      
      %% Set map axis
      if get(obj.h_gui.zoom_selectCB,'Value') && isfinite(xlims(1))
        axis(obj.h_axes,'normal');
        
        xlims = xlims/1e3;
        ylims = ylims/1e3;
        cur_xlim = xlim(obj.h_axes);
        cur_ylim = ylim(obj.h_axes);
        axes_position = get(obj.h_axes,'Position');
        aspect_ratio = axes_position(3) / axes_position(4);
        xlims(1) = min(xlims(1),cur_xlim(1));
        xlims(2) = max(xlims(2),cur_xlim(2));
        ylims(1) = min(ylims(1),cur_ylim(1));
        ylims(2) = max(ylims(2),cur_ylim(2));
        
        if diff(xlims) / diff(ylims) > aspect_ratio
          ylims = mean(ylims) + diff(xlims)/aspect_ratio/2 * [-1 1];
        else
          xlims = mean(xlims) + diff(ylims)*aspect_ratio/2 * [-1 1];
        end
        xlim(obj.h_axes,xlims);
        ylim(obj.h_axes,ylims);
      end
      
      
    end
    
    function newPB_callback(obj,h_obj,event)
      prompt = {'Enter new flight line name:'};
      dlg_title = 'Flight line name';
      num_lines = 1;
      def = {'new_fline'};
      answer = inputdlg(prompt,dlg_title,num_lines,def);
      
      if ~isempty(answer) && ~isempty(answer{1})
        answer = answer{1};
        fline.lat = [];
        fline.lon = [];
        fline.wpnt_names = {};
        fline.name = answer;
        obj.insert_fline(fline,[],true);
        % Select the flight line that was just created
        set(obj.h_gui.flines.listLB,'Value',length(obj.flines));
      end
    end
    
    function openPB_callback(obj,h_obj,event)
      [open_fns, open_fn_dir, filterindex] = uigetfile( ...
        {'*.mat','MAT-files (*.mat)'; ...
        '*.shp','Shape files (*.shp)'; ...
        '*.csv',  'CSV Files (*.csv)'; ...
        '*.txt',  'TXT Files (*.txt)'; ...
        '*.kml',  'KML Files (*.kml)'; ...
        '*.tif',  'Geotiff Files (*.tif)'; ...
        '*.*',  'All Files (*.*)'}, ...
        'Pick a file', ...
        obj.open_fn_dir, 'MultiSelect', 'on');
      
      if isequal(open_fns,0) || isequal(open_fn_dir,0)
        return;
      end
      
      obj.open_fn_dir = open_fn_dir;
      
      if ischar(open_fns)
        open_fns = {open_fns};
      end
      default_x_field = '1';
      default_y_field = '2';
      default_name_field = '';
      shape_file_coordinates = 'map';
      other_file_coordinates = 'geodetic';
      default_delim = '';
      default_num_header_lines = '0';
      default_plot_params = 'k.';
      default_map_only = 'false';
      default_unload_field = '1';
      default_geotiff_field = '1';
      default_plotonly_field = '1';
      default_flightlines_field = '1';
      for file_idx = 1:length(open_fns)
        fn = fullfile(open_fn_dir, open_fns{file_idx});
        
        [~,~,ext] = fileparts(fn);
        
        if any(strcmpi(ext,{'.csv','.txt'}))
          fid = fopen(fn,'r');
          header_line = fgets(fid);
          if any(header_line == ',')
            default_delim = ',';
          end
          fprintf('%s', header_line);
          line_str = fgets(fid);
          fprintf('%s', line_str);
          
          id_field = [];
          if strcmpi(ext,'.csv')
            default_delim = ',';
            default_num_header_lines = '1';
            
            fields = regexpi(header_line,default_delim,'split');
            for field_idx=1:length(fields)
              fields{field_idx} = fields{field_idx}(fields{field_idx} ~= 13 & fields{field_idx} ~= 10);
            end
            lat_field = find(~cellfun(@isempty,regexpi(fields,'.*lat.*')),1);
            if ~isempty(lat_field)
              default_y_field = sprintf('%d',lat_field);
            end
            
            lon_field = find(~cellfun(@isempty,regexpi(fields,'.*lon.*')),1);
            if ~isempty(lon_field)
              default_x_field = sprintf('%d',lon_field);
            end
            
            id_field = find(~cellfun(@isempty,regexpi(fields,'.*name.*')),1);
            if isempty(id_field)
              id_field = find(~cellfun(@isempty,regexpi(fields,'.*id.*')),1);
              if ~isempty(id_field)
                default_name_field = sprintf('%d',id_field);
              end
            else
              default_name_field = sprintf('%d',id_field);
            end
          end
          
          prompt = {'X/Lon Field #:','Y/Lat Field #', ...
            'Name/ID Field # (blank for none)', '(G)eodetic/(M)ap coordinates:', ...
            'Delimiter', '# of Header Lines', 'Matlab plot parameter', ...
            'Plot map only? (true/false)'};
          dlg_title = 'Define Read Parameters:';
          num_lines = 1;
          def = {default_x_field, default_y_field, ...
            default_name_field, other_file_coordinates, ...
            default_delim, default_num_header_lines, default_plot_params, ...
            default_map_only};
          answer = inputdlg(prompt,dlg_title,num_lines,def);
          
          if length(answer) == 8 && ~isempty(answer{1}) && ~isempty(answer{2})
            default_x_field = answer{1};
            default_y_field = answer{2};
            default_name_field = answer{3};
            if isempty(answer{4}) || ~isempty(regexpi(answer{4},'g'))
              other_file_coordinates = 'Geodetic';
            else
              other_file_coordinates = 'Map';
            end
            default_delim = answer{5};
            default_num_header_lines = answer{6};
            default_plot_params = answer{7};
            default_map_only = answer{8};
            
            fseek(fid,0,-1);
            for idx = 1:eval(default_num_header_lines)
              line_str = fgets(fid);
            end
            line_str = fgets(fid);
            if isempty(default_delim)
              C = textscan(line_str,'%s');
            else
              C = textscan(line_str,'%s','Delimiter',default_delim);
            end
            num_fields = length(C{1});
            
            fseek(fid,0,-1);
            textscan_str = '';
            try
              x_idx = str2double(default_x_field);
            end
            try
              y_idx = str2double(default_y_field);
            end
            try
              name_idx = str2double(default_name_field);
            end
            x_found = false;
            y_found = false;
            name_found = false;
            for idx = 1:num_fields
              if idx == x_idx
                textscan_str = cat(2,textscan_str,'%f');
                x_found = true;
              elseif idx == y_idx
                textscan_str = cat(2,textscan_str,'%f');
                y_found = true;
              elseif idx == name_idx
                textscan_str = cat(2,textscan_str,'%s');
                name_found = true;
              else
                textscan_str = cat(2,textscan_str,'%s');
              end
            end
            if ~x_found || ~y_found
              warning('X or Y coordinate missing');
              return;
            end
            if isempty(default_delim)
              C = textscan(fid, textscan_str, ...
                'HeaderLines',eval(default_num_header_lines));
            else
              C = textscan(fid, textscan_str,'Delimiter',default_delim, ...
                'HeaderLines',eval(default_num_header_lines));
            end
            fclose(fid);
            
            if strcmpi(other_file_coordinates,'map')
              x = C{x_idx};
              y = C{y_idx};
            else
              [x,y] = projfwd(obj.proj,C{y_idx},C{x_idx});
            end
            x = reshape(x,[1 length(x)]);
            y = reshape(y,[1 length(y)]);
            
            try
              map_only_flag = eval(default_map_only)
            catch
              map_only_flag = true;
            end
            if ~map_only_flag
              cur_fline = length(obj.flines)+1;
              [fline.lat,fline.lon] = projinv(obj.proj,x,y);
              if ~name_found
                fline.wpnt_names = {};
                for wpnt_idx = 1:length(fline.lat)
                  fline.wpnt_names{wpnt_idx} = sprintf('S_%d', wpnt_idx);
                end
              else
                fline.wpnt_names = C{name_idx};
              end
              fline.name = open_fns{file_idx};
              obj.insert_fline(fline,[],true);
            else
              hold(obj.h_axes, 'on');
              % Store the geodetic coordinates in the plot handle so we can
              % update the projection later if we need to.
              [userdata.lat,userdata.lon] = projinv(obj.proj,x,y);
              new_plotonly = plot(x/1e3,y/1e3,default_plot_params,'Parent',obj.h_axes,'UserData',userdata);
              obj.plotonly.insert_handle(new_plotonly);
            end
          end
          
        elseif strcmpi(ext,'.kml')
          kml = kml2struct(fn);
          
          % Load each line in... match each point in the line with
          % corresponding point in file and rename waypoint in line
          line_idxs = find(strcmp('Line',{kml.Geometry}));
          point_idxs = find(strcmp('Point',{kml.Geometry}));
          
          for line_idx = line_idxs
            fline = [];
            fline.wpnt_names = {};
            fline.lon = kml(line_idx).Lon;
            fline.lat = kml(line_idx).Lat;
            for point_idx = 1:length(kml(line_idx).Lon)
              match_idx = find(kml(line_idx).Lon(point_idx) == cell2mat({kml(point_idxs).Lon}) ...
                & kml(line_idx).Lat(point_idx) == cell2mat({kml(point_idxs).Lat}));
              if ~isempty(match_idx)
                match_idx = match_idx(1);
                fline.wpnt_names{point_idx} = kml(point_idxs(match_idx)).Name;
              else
                fline.wpnt_names{point_idx} = sprintf('K_%d',point_idx);
              end
            end
            fline.name = kml(line_idx).Name;
            obj.insert_fline(fline,[],true);
          end
          
        elseif strcmpi(ext,'.shp')
          % Shape file lines sorted alphabetically?
          % filename_name (name is 1.1, 1.2t, etc)
          S = shaperead(fn);
          
          fields = fieldnames(S);
          for field_idx = 1:length(fields)
            fprintf('%s\n', fields{field_idx})
          end
          id_field = find(~cellfun(@isempty,regexpi(fields,'.*id.*')),1)
          if isempty(id_field)
            id_field = '';
          else
            id_field = fields{id_field};
          end
          
          prompt = {'X/Longitude:','Y/Latitude:','Name:','(Map) or (Geo)graphic?'};
          dlg_title = 'Define Shapefield Fields';
          num_lines = 1;
          def = {'X','Y',id_field,shape_file_coordinates};
          answer = inputdlg(prompt,dlg_title,num_lines,def);
          
          if length(answer) == 4 && ~isempty(answer{1}) && ~isempty(answer{2})
            if isempty(answer{4}) || ~isempty(regexpi(answer{4},'g'))
              shape_file_coordinates = 'Geometry';
            else
              shape_file_coordinates = 'Map';
            end
            for idx = 1:length(S)
              if strcmpi(shape_file_coordinates,'map')
                [fline.lat fline.lon] ...
                  = projinv(obj.proj, S(idx).(answer{1}), S(idx).(answer{2}));
              else
                fline.lat = S(idx).(answer{2});
                fline.lon = S(idx).(answer{1});
              end
              fline.lat = fline.lat(1:end-1);
              fline.lon = fline.lon(1:end-1);
              fline.wpnt_names = {};
              for wpnt_idx = 1:length(fline.lat)
                fline.wpnt_names{wpnt_idx} = sprintf('S_%d', wpnt_idx);
              end
              if ~isempty(answer{3})
                fline_name_fields = regexpi(answer{3},':','split');
                if isnumeric(S(idx).(fline_name_fields{1}))
                  fline.name = sprintf('%d',S(idx).(fline_name_fields{1}));
                else
                  fline.name = S(idx).(fline_name_fields{1});
                end
                for field_idx = 2:length(fline_name_fields)
                  if isnumeric(S(idx).(fline_name_fields{field_idx}))
                    fline.name = sprintf('%s_%d', fline.name, S(idx).(fline_name_fields{field_idx}));
                  else
                    fline.name = sprintf('%s_%s', fline.name, S(idx).(fline_name_fields{field_idx}));
                  end
                end
              else
                fline.name = sprintf('%s_%d', open_fns{file_idx}, idx);
              end
              obj.insert_fline(fline,[],true);
            end
          end
          
        elseif strcmpi(ext,'.mat')
          prompt = {'Unload all first','Load geotiff','Load plotonly', ...
            'Load flightlines'};
          dlg_title = 'Define which values to load:';
          num_lines = 1;
          def = {default_unload_field, default_geotiff_field, default_plotonly_field, ...
            default_flightlines_field};
          answer = inputdlg(prompt,dlg_title,num_lines,def);
          
          if length(answer) == 4
            default_unload_field = answer{1};
            default_geotiff_field = answer{2};
            default_plotonly_field = answer{3};
            default_flightlines_field = answer{4};
            try
              unload_field = eval(default_unload_field);
              if (unload_field)
                unload_field = true;
              else
                unload_field = false;
              end
            catch
              unload_field = true;
            end
            try
              load_geotiff_field = eval(default_geotiff_field);
              if (load_geotiff_field)
                load_geotiff_field = true;
              else
                load_geotiff_field = false;
              end
            catch
              load_geotiff_field = true;
            end
            try
              load_plotonly_field = eval(default_plotonly_field);
              if (load_plotonly_field)
                load_plotonly_field = true;
              else
                load_plotonly_field = false;
              end
            catch
              load_plotonly_field = true;
            end
            try
              load_flightlines_field = eval(default_flightlines_field);
              if (load_flightlines_field)
                load_flightlines_field = true;
              else
                load_flightlines_field = false;
              end
            catch
              load_flightlines_field = true;
            end
            
            obj.openMatFile(fn,unload_field,load_geotiff_field,load_plotonly_field,load_flightlines_field);
          end
          
        elseif strcmpi(ext,'.tif')
          % Update geotiff
          obj.geotiff_fn = fn;
          obj.update_geotiff(true);
        end
      end
    end
    
    function openMatFile(obj,fn,unload_field,load_geotiff_field,load_plotonly_field,load_flightlines_field)
      if ~exist(fn,'file')
        warning('Could not find file %s\n', obj.save_fn);
        return;
      end
      new_data = load(fn);
      obj.save_fn = fn;
      title(obj.save_fn,'interpreter','none','parent',obj.h_axes);
      % Unload all data
      if unload_field
        delete_flines(obj,ones(size(obj.flines)));
        obj.plotonly.delete_handle(1:length(obj.plotonly.handles));
      end
      % Update geotiff
      if isempty(obj.proj) || (load_geotiff_field && isfield(new_data,'geotiff_fn') ...
          && exist(new_data.geotiff_fn,'file') ...
          && ~strcmpi(obj.geotiff_fn,new_data.geotiff_fn))
        obj.geotiff_fn = new_data.geotiff_fn;
        obj.update_geotiff(true);
      end
      if load_plotonly_field
        % Insert all plotonly graphics
        new_plotonly = struct2handle(new_data.plotonly,obj.h_axes);
        
        obj.plotonly.insert_handle(new_plotonly);
      end
      if load_flightlines_field
        % Insert all flight lines
        for pos = 1:length(new_data.flines)
          obj.insert_fline(new_data.flines(pos),[],0);
        end
      end
    end
    
    function savePB_callback(obj,h_obj,event)
      if isempty(obj.save_fn) || exist(obj.save_fn,'dir')
        saveasPB_callback(obj,h_obj,event)
      else
        fprintf('Saving flight lines to %s\n',obj.save_fn);
        geotiff_fn = obj.geotiff_fn;
        flines = obj.flines;
        plotonly = handle2struct(obj.plotonly.handles);
        save(obj.save_fn, 'flines', 'geotiff_fn', 'plotonly')
      end
    end
    
    function saveasPB_callback(obj,h_obj,event)
      [filename, pathname] = uiputfile( ...
        {'*.mat', 'All MATLAB Files (*.mat)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Save as', obj.save_fn);
      if isempty(filename) || filename(1) == 0
        return;
      end
      
      obj.save_fn = fullfile(pathname, filename);
      title(obj.save_fn,'interpreter','none','parent',obj.h_axes);
      fprintf('Saving flight lines to %s\n',obj.save_fn);
      geotiff_fn = obj.geotiff_fn;
      flines = obj.flines;
      plotonly = handle2struct(obj.plotonly.handles);
      save(obj.save_fn, 'flines', 'geotiff_fn', 'plotonly')
    end
    
    function insert_fline(obj,fline,pos,selected)
      % Required Inputs: fline.name, fline.lat, fline.lon, fline.wpnt_names
      if isempty(pos)
        pos = length(obj.flines) + 1;
      end
      
      cur_fline_names = get(obj.h_gui.flines.listLB,'String');
      cur_fline_names = [cur_fline_names(1:pos-1);
        {fline.name}; cur_fline_names(pos:end)];
      set(obj.h_gui.flines.listLB,'String',cur_fline_names);
      
      cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
      if iscell(cur_fline_selected)
        cur_fline_selected = cell2mat(cur_fline_selected);
      end
      cur_fline_selected = sort([cur_fline_selected pos]);
      set(obj.h_gui.flines.listLB,'Value',cur_fline_selected);
      
      % Map fields and make sure lat,lon,wpnt_names are 1xN
      obj.flines(pos+1:end+1) = obj.flines(pos:end);
      obj.flines(pos).lat = fline.lat(:).';
      obj.flines(pos).lon = fline.lon(:).';
      obj.flines(pos).wpnt_names = fline.wpnt_names(:).';
      obj.flines(pos).name = fline.name;
      
      hold(obj.h_axes,'on')
      [obj.flines(pos).x obj.flines(pos).y] ...
        = projfwd(obj.proj, obj.flines(pos).lat, obj.flines(pos).lon);
      if length(obj.flines(pos).x) >= 1
        obj.flines(pos).handle ...
          = plot(obj.flines(pos).x/1e3,obj.flines(pos).y/1e3,'b.-','Parent',obj.h_axes);
        obj.flines(pos).start_handle ...
          = plot(obj.flines(pos).x(1)/1e3,obj.flines(pos).y(1)/1e3,'bo','Parent',obj.h_axes);
      else
        obj.flines(pos).handle = plot(0,0,'b.-','Parent',obj.h_axes);
        uistack(obj.flines(pos).handle, 'top');
        set(obj.flines(pos).handle,'XData',[]);
        set(obj.flines(pos).handle,'YData',[]);
        obj.flines(pos).start_handle = plot(0,0,'bo','Parent',obj.h_axes);
        set(obj.flines(pos).start_handle,'XData',[]);
        set(obj.flines(pos).start_handle,'YData',[]);
      end
      hold(obj.h_axes,'off')
      
      obj.flines(pos).along_track ...
        = geodetic_to_along_track(obj.flines(pos).lat, ...
        obj.flines(pos).lon);
      
      obj.flinesLB_callback();
    end
    
    function insert_wpnt(obj,x,y,name)
      cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
      if iscell(cur_fline_selected)
        cur_fline_selected = cell2mat(cur_fline_selected);
      end
      if isempty(cur_fline_selected)
        warning('Can not insert: no flight line selected\n');
        return;
      end
      cur_fline_selected = cur_fline_selected(1);
      
      cur_wpnt_selected = get(obj.h_gui.wpnts.listLB,'Value');
      if iscell(cur_wpnt_selected)
        cur_wpnt_selected = cell2mat(cur_wpnt_selected);
      end
      if isempty(cur_wpnt_selected)
        % Insert before first element if nothing selected
        pos = 1;
      else
        % Insert after currently selected element
        pos = cur_wpnt_selected(1) + 1;
      end
      obj.flines(cur_fline_selected).x = [obj.flines(cur_fline_selected).x(1:pos-1) 0 obj.flines(cur_fline_selected).x(pos:end)];
      obj.flines(cur_fline_selected).y = [obj.flines(cur_fline_selected).y(1:pos-1) 0 obj.flines(cur_fline_selected).y(pos:end)];
      obj.flines(cur_fline_selected).wpnt_names = [obj.flines(cur_fline_selected).wpnt_names(1:pos-1), {''}, obj.flines(cur_fline_selected).wpnt_names(pos:end)];
      obj.flines(cur_fline_selected).x(pos) = x*1e3;
      obj.flines(cur_fline_selected).y(pos) = y*1e3;
      done = false;
      if exist('name','var') ...
          && all(~strcmpi(name,obj.flines(cur_fline_selected).wpnt_names))
        wpnt_name = name;
        done = true;
      end
      way_pnt_idx = 1;
      while ~done
        if ~exist('name','var')
          wpnt_name = sprintf('%d', way_pnt_idx);
        else
          wpnt_name = sprintf('%s_%d', name, way_pnt_idx);
        end
        if all(~strcmpi(wpnt_name,obj.flines(cur_fline_selected).wpnt_names))
          done = true;
        else
          way_pnt_idx = way_pnt_idx + 1;
        end
      end
      obj.flines(cur_fline_selected).wpnt_names(pos) = {wpnt_name};
      
      [obj.flines(cur_fline_selected).lat obj.flines(cur_fline_selected).lon] ...
        = projinv(obj.proj,obj.flines(cur_fline_selected).x,obj.flines(cur_fline_selected).y);
      obj.flines(cur_fline_selected).along_track ...
        = geodetic_to_along_track(obj.flines(cur_fline_selected).lat, ...
        obj.flines(cur_fline_selected).lon);
      
      cur_wpnt_names = get(obj.h_gui.wpnts.listLB,'String');
      cur_wpnt_names = [cur_wpnt_names(1:pos-1);
        obj.flines(cur_fline_selected).wpnt_names(pos); cur_wpnt_names(pos:end)];
      set(obj.h_gui.wpnts.listLB,'String',cur_wpnt_names);
      
      set(obj.h_gui.wpnts.listLB,'Value',pos);
      
      set(obj.flines(cur_fline_selected).handle,'XData',obj.flines(cur_fline_selected).x/1e3,'YData',obj.flines(cur_fline_selected).y/1e3);
    end
    
    function reset_copy(obj)
      obj.selection.copy.wpnts = [];
      obj.selection.copy.flines = [];
    end
    
    function wpntsLB_menu_callback(obj,h_obj,event)
      command = get(h_obj,'Label');
      
      if strcmpi(command,'Change Length')
        cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
        if iscell(cur_fline_selected)
          cur_fline_selected = cell2mat(cur_fline_selected);
        end
        if isempty(cur_fline_selected)
          warning('No flight line selected');
          return;
        end
        cur_fline_selected = cur_fline_selected(1);
        
        if length(obj.flines(cur_fline_selected).lat) < 1
          warning('Can not shorten this line, only one waypoint');
          return;
        end
        
        cur_wpnt_selected = get(obj.h_gui.wpnts.listLB,'Value');
        if iscell(cur_wpnt_selected)
          cur_wpnt_selected = cell2mat(cur_wpnt_selected);
        end
        if isempty(cur_wpnt_selected)
          warning('No waypoint selected');
          return;
        end
        cur_wpnt_selected = cur_wpnt_selected(1);
        
        if cur_wpnt_selected == 1
          ref_wpnt = 1;
        else
          ref_wpnt = -1;
        end
        
        % Get the length to the previous waypoint
        current_length = geodetic_to_along_track( ...
          obj.flines(cur_fline_selected).lat([cur_wpnt_selected+ref_wpnt cur_wpnt_selected]), ...
          obj.flines(cur_fline_selected).lon([cur_wpnt_selected+ref_wpnt cur_wpnt_selected]));
        current_length = current_length(2);
        
        done = false;
        while ~done
          prompt = {'Length (km):','Relative reference waypoint:'};
          def = {sprintf('%g',current_length/1e3),sprintf('%d',ref_wpnt)};
          dlg_title = 'Change Length';
          num_lines = 1;
          answer = inputdlg(prompt,dlg_title,num_lines,def);
          
          if length(answer) == 2
            if ~isempty(answer{1}) && ~isempty(answer{2})
              done = true;
            else
              ref_pnt = str2double(answer{2});
              current_length = geodetic_to_along_track( ...
                obj.flines(cur_fline_selected).lat([cur_wpnt_selected+ref_wpnt cur_wpnt_selected]), ...
                obj.flines(cur_fline_selected).lon([cur_wpnt_selected+ref_wpnt cur_wpnt_selected]));
              current_length = current_length(2);
            end
          else
            return;
          end
        end
        
        if length(answer) == 2 && ~isempty(answer{1}) && ~isempty(answer{2})
          ref_pnt = str2double(answer{2});
          if answer{1}(1) == 'r'
            current_length = geodetic_to_along_track( ...
              obj.flines(cur_fline_selected).lat([cur_wpnt_selected+ref_wpnt cur_wpnt_selected]), ...
              obj.flines(cur_fline_selected).lon([cur_wpnt_selected+ref_wpnt cur_wpnt_selected]));
            current_length = current_length(2);
            new_length = current_length + str2double(answer{1}) * 1e3;
          else
            new_length = str2double(answer{1}) * 1e3;
          end
          
          unit_vector = [diff(obj.flines(cur_fline_selected).x([cur_wpnt_selected+ref_pnt cur_wpnt_selected]));
            diff(obj.flines(cur_fline_selected).y([cur_wpnt_selected+ref_pnt cur_wpnt_selected]))];
          unit_vector = unit_vector ./ sqrt(dot(unit_vector,unit_vector));
          
          x0 = obj.flines(cur_fline_selected).x(cur_wpnt_selected+ref_pnt);
          y0 = obj.flines(cur_fline_selected).y(cur_wpnt_selected+ref_pnt);
          obj.flines(cur_fline_selected).x(cur_wpnt_selected) ...
            = x0 + unit_vector(1)*new_length;
          obj.flines(cur_fline_selected).y(cur_wpnt_selected) ...
            = y0 + unit_vector(2)*new_length;
          [obj.flines(cur_fline_selected).lat(cur_wpnt_selected) ...
            obj.flines(cur_fline_selected).lon(cur_wpnt_selected)] ...
            = projinv(obj.proj, obj.flines(cur_fline_selected).x(cur_wpnt_selected), ...
            obj.flines(cur_fline_selected).y(cur_wpnt_selected));
          
          obj.update_flineGraphics();
        end
        
      elseif any(strcmpi(command,{'Edit Geographic','Edit Map'}))
        cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
        if iscell(cur_fline_selected)
          cur_fline_selected = cell2mat(cur_fline_selected);
        end
        if isempty(cur_fline_selected)
          warning('No flight line selected');
          return;
        end
        cur_fline_selected = cur_fline_selected(1);
        
        cur_wpnt_selected = get(obj.h_gui.wpnts.listLB,'Value');
        if iscell(cur_wpnt_selected)
          cur_wpnt_selected = cell2mat(cur_wpnt_selected);
        end
        if isempty(cur_wpnt_selected)
          warning('No waypoint selected');
          return;
        end
        cur_wpnt_selected = cur_wpnt_selected(1);
        
        if strcmpi(command,{'Edit Geographic'})
          prompt = {'Latitude:','Longitude:','Name:'};
          def = {sprintf('%.12g',obj.flines(cur_fline_selected).lat(cur_wpnt_selected)), ...
            sprintf('%.12g',obj.flines(cur_fline_selected).lon(cur_wpnt_selected)), ...
            obj.flines(cur_fline_selected).wpnt_names{cur_wpnt_selected}};
        else
          prompt = {'X:','Y:','Name:'};
          def = {sprintf('%.12g',obj.flines(cur_fline_selected).x(cur_wpnt_selected)), ...
            sprintf('%.12g',obj.flines(cur_fline_selected).y(cur_wpnt_selected)), ...
            obj.flines(cur_fline_selected).wpnt_names{cur_wpnt_selected}};
        end
        dlg_title = 'Edit waypoint';
        num_lines = 1;
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        
        if length(answer) == 3 && ~isempty(answer{1}) && ~isempty(answer{2})
          if strcmpi(command,{'Edit Geographic'})
            lat = eval(answer{1});
            lon = eval(answer{2});
            [x y] ...
              = projfwd(obj.proj, lat, lon);
          else
            x = eval(answer{1});
            y = eval(answer{2});
            [lat lon] ...
              = projfwd(obj.proj, x, y);
          end
          obj.flines(cur_fline_selected).x(cur_wpnt_selected) = x;
          obj.flines(cur_fline_selected).y(cur_wpnt_selected) = y;
          obj.flines(cur_fline_selected).wpnt_names(cur_wpnt_selected) = answer(3);
          obj.flines(cur_fline_selected).lat(cur_wpnt_selected) = lat;
          obj.flines(cur_fline_selected).lon(cur_wpnt_selected) = lon;
          obj.flines(cur_fline_selected).along_track ...
            = geodetic_to_along_track(obj.flines(cur_fline_selected).lat, ...
            obj.flines(cur_fline_selected).lon);
          set(obj.flines(cur_fline_selected).handle,'XData',obj.flines(cur_fline_selected).x/1e3,'YData',obj.flines(cur_fline_selected).y/1e3);
          
          obj.wpntsLB_callback();
          obj.update_statusText();
        end
        
      elseif strcmpi(command,'Insert Geographic')
        obj.reset_copy()
        prompt = {'Latitude:','Longitude:','Name:'};
        dlg_title = 'New waypoint';
        num_lines = 1;
        def = {'','',''};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        
        if length(answer) == 3 && ~isempty(answer{1}) && ~isempty(answer{2})
          [x y] ...
            = projfwd(obj.proj, eval(answer{1}), eval(answer{2}));
          obj.insert_wpnt(x/1e3,y/1e3,answer{3});
        end
        obj.wpntsLB_callback();
        obj.update_statusText();
        
      elseif strcmpi(command,'Insert Map')
        obj.reset_copy()
        prompt = {'X (km):','Y (km):','Name:'};
        dlg_title = 'New waypoint';
        num_lines = 1;
        def = {'','',''};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        
        if length(answer) == 3 && ~isempty(answer{1}) && ~isempty(answer{2})
          obj.insert_wpnt(eval(answer{1}),eval(answer{2}),answer{3});
        end
        obj.wpntsLB_callback();
        obj.update_statusText();
        
      elseif strcmpi(command,'Copy')
        cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
        if iscell(cur_fline_selected)
          cur_fline_selected = cell2mat(cur_fline_selected);
        end
        if isempty(cur_fline_selected)
          warning('No flight line selected');
          return;
        end
        cur_fline_selected = cur_fline_selected(1);
        obj.selection.copy.flines = cur_fline_selected;
        
        cur_wpnt_selected = get(obj.h_gui.wpnts.listLB,'Value');
        if iscell(cur_wpnt_selected)
          cur_wpnt_selected = cell2mat(cur_wpnt_selected);
        end
        if isempty(cur_wpnt_selected)
          warning('No waypoint selected');
          return;
        end
        obj.selection.copy.wpnts = cur_wpnt_selected;
        
      elseif strcmpi(command,'Cut')
      elseif strcmpi(command,'Paste')
        if isempty(obj.selection.copy.flines)
          warning('No flight line selected');
          return;
        end
        cur_fline_selected = obj.selection.copy.flines(1);
        fline = obj.flines(cur_fline_selected);
        for cur_wpnt_selected = obj.selection.copy.wpnts
          x = fline.x(cur_wpnt_selected);
          y = fline.y(cur_wpnt_selected);
          name = fline.wpnt_names(cur_wpnt_selected);
          for idx = 1:length(x)
            obj.insert_wpnt(x(idx)/1e3,y(idx)/1e3,name{idx});
          end
        end
        obj.wpntsLB_callback();
        obj.update_statusText();
        
      elseif strcmpi(command,'Delete')
        obj.reset_copy()
        delete_pointPB_callback(obj,h_obj,event)
        
      elseif strcmpi(command,'Rename') || strcmpi(command,'Rename Auto')
        cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
        if iscell(cur_fline_selected)
          cur_fline_selected = cell2mat(cur_fline_selected);
        end
        if isempty(cur_fline_selected)
          warning('No flight line selected');
          return;
        end
        cur_fline_selected = cur_fline_selected(1);
        
        cur_wpnt_selected = get(obj.h_gui.wpnts.listLB,'Value');
        if iscell(cur_wpnt_selected)
          cur_wpnt_selected = cell2mat(cur_wpnt_selected);
        end
        if isempty(cur_wpnt_selected)
          warning('No waypoint selected');
          return;
        end
        
        cur_names = get(obj.h_gui.wpnts.listLB,'String');
        if strcmpi(command,'Rename')
          prompt = {'Rename waypoint:'};
          dlg_title = 'Waypoint rename';
          def = cur_names(cur_wpnt_selected(1));
        else
          prompt = {'Auto rename waypoint (base):','Repeat if same: '};
          dlg_title = 'Waypoint auto rename';
          def = [cur_names(cur_wpnt_selected(1)) {'true'}];
        end
        num_lines = 1;
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        
        if length(answer)>=1 && ~isempty(answer{1})
          if strcmpi(command,'Rename')
            answer = answer{1};
            cur_names{cur_wpnt_selected} = answer;
            set(obj.h_gui.wpnts.listLB,'String',cur_names);
            obj.flines(cur_fline_selected).wpnt_names{cur_wpnt_selected} = answer;
          else
            if length(answer)>=2 && ~isempty(answer{2})
              try
                repeat_if_same = eval(answer{2});
              catch
                repeat_if_same = true;
              end
            end
            answer = answer{1};
            wpnt_idx_rename = 1;
            for wpnt_idx = 1:length(cur_wpnt_selected)
              wpnt = cur_wpnt_selected(wpnt_idx);
              if repeat_if_same
                wpnt_is_same = false;
                for same_wpnt_idx = 1:wpnt_idx-1
                  same_wpnt = cur_wpnt_selected(same_wpnt_idx);
                  if (obj.flines(cur_fline_selected).lat(wpnt) ...
                      == obj.flines(cur_fline_selected).lat(same_wpnt) ...
                      && obj.flines(cur_fline_selected).lon(wpnt) ...
                      == obj.flines(cur_fline_selected).lon(same_wpnt))
                    wpnt_is_same = true;
                    break;
                  end
                end
              else
                wpnt_is_same = false;
              end
              if wpnt_is_same
                cur_names{wpnt} = obj.flines(cur_fline_selected).wpnt_names{same_wpnt};
              else
                cur_names{wpnt} = sprintf('%s%d', answer, wpnt_idx_rename);
                wpnt_idx_rename = wpnt_idx_rename + 1;
              end
              obj.flines(cur_fline_selected).wpnt_names{wpnt} = cur_names{wpnt};
            end
            set(obj.h_gui.wpnts.listLB,'String',cur_names);
          end
        end
        
      elseif strcmpi(command,'Reverse Order')
        obj.reset_copy()
        cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
        if iscell(cur_fline_selected)
          cur_fline_selected = cell2mat(cur_fline_selected);
        end
        if isempty(cur_fline_selected)
          warning('No flight line selected');
          return;
        end
        cur_fline_selected = cur_fline_selected(1);
        
        cur_wpnt_selected = get(obj.h_gui.wpnts.listLB,'Value');
        if iscell(cur_wpnt_selected)
          cur_wpnt_selected = cell2mat(cur_wpnt_selected);
        end
        if isempty(cur_wpnt_selected)
          cur_wpnt_selected = 1:length(obj.flines(cur_fline_selected));
        end
        
        obj.flines(cur_fline_selected).lat(cur_wpnt_selected) = obj.flines(cur_fline_selected).lat(cur_wpnt_selected(end:-1:1));
        obj.flines(cur_fline_selected).lon(cur_wpnt_selected) = obj.flines(cur_fline_selected).lon(cur_wpnt_selected(end:-1:1));
        obj.flines(cur_fline_selected).wpnt_names(cur_wpnt_selected) = obj.flines(cur_fline_selected).wpnt_names(cur_wpnt_selected(end:-1:1));
        obj.flines(cur_fline_selected).x(cur_wpnt_selected) = obj.flines(cur_fline_selected).x(cur_wpnt_selected(end:-1:1));
        obj.flines(cur_fline_selected).y(cur_wpnt_selected) = obj.flines(cur_fline_selected).y(cur_wpnt_selected(end:-1:1));
        
        obj.flines(cur_fline_selected).along_track ...
          = geodetic_to_along_track(obj.flines(cur_fline_selected).lat, ...
          obj.flines(cur_fline_selected).lon);
        
        cur_wpnt_names = get(obj.h_gui.wpnts.listLB,'String');
        cur_wpnt_names(cur_wpnt_selected) = cur_wpnt_names(cur_wpnt_selected(end:-1:1));
        set(obj.h_gui.wpnts.listLB,'String',cur_wpnt_names);
        
        obj.update_flineGraphics();
        
      elseif strcmpi(command,'Rotate')
        cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
        obj.rotate_waypoints(cur_fline_selected, false);
      end
    end
    
    function update_flineGraphics(obj)
      cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
      if iscell(cur_fline_selected)
        cur_fline_selected = cell2mat(cur_fline_selected);
      end
      
      cur_wpnt_selected = get(obj.h_gui.wpnts.listLB,'Value');
      if iscell(cur_wpnt_selected)
        cur_wpnt_selected = cell2mat(cur_wpnt_selected);
      end
      
      axes_children = get(obj.h_axes,'children');
      children_mask = zeros(size(axes_children));
      for pos = 1:length(obj.flines)
        if ~isempty(cur_fline_selected) && pos == cur_fline_selected(1)
          cur_wpnt_selected = cur_wpnt_selected(cur_wpnt_selected <= length(obj.flines(pos).x));
          
          % Update 'X' marks
          set(obj.selection.h_wpnt_plot,'XData',obj.flines(pos).x(cur_wpnt_selected)/1e3,'YData',obj.flines(pos).y(cur_wpnt_selected)/1e3);
          children_mask(axes_children == obj.selection.h_wpnt_plot) = -3;
          
          % Update way points listbox names
          cur_wpnt_names = cell(size(obj.flines(pos).wpnt_names));
          [cur_wpnt_names{:}] = deal(obj.flines(pos).wpnt_names{:});
          set(obj.h_gui.wpnts.listLB,'String',cur_wpnt_names);
        end
        
        % Update flight lines
        set(obj.flines(pos).handle,'XData',obj.flines(pos).x/1e3,'YData',obj.flines(pos).y/1e3);
        if length(obj.flines(pos).x) > 0
          set(obj.flines(pos).start_handle,'XData',obj.flines(pos).x(1)/1e3,'YData',obj.flines(pos).y(1)/1e3);
        else
          set(obj.flines(pos).start_handle,'XData',[],'YData',[]);
        end
        
        if any(pos == cur_fline_selected)
          set(obj.flines(pos).handle,'Color','Red');
          set(obj.flines(pos).start_handle,'Color','Red');
          children_mask(axes_children == obj.flines(pos).handle) = -1;
          children_mask(axes_children == obj.flines(pos).start_handle) = -2;
        else
          set(obj.flines(pos).handle,'Color','Blue');
          set(obj.flines(pos).start_handle,'Color','Blue');
        end
      end
      [~,children_sort_idx] = sort(children_mask);
      set(obj.h_axes,'children',axes_children(children_sort_idx));
      
      % Update flines listbox names
      cur_fline_names = cell(size(obj.flines));
      if ~isempty(obj.flines)
        [cur_fline_names{:}] = deal(obj.flines(:).name);
      end
      set(obj.h_gui.flines.listLB,'String',cur_fline_names);
      
      % Set ListTopBox
      cur_top_value = get(obj.h_gui.wpnts.listLB,'ListboxTop');
      if all(cur_wpnt_selected - cur_top_value > 5) ...
          || all(cur_wpnt_selected - cur_top_value < 0)
        if ~isempty(cur_wpnt_selected)
          set(obj.h_gui.wpnts.listLB,'ListboxTop',min(cur_wpnt_selected))
        else
          set(obj.h_gui.wpnts.listLB,'ListboxTop',1)
        end
      end
      obj.update_statusText();
    end
    
    function flinesLB_menu_callback(obj,h_obj,event)
      command = get(h_obj,'Label');
      
      if strcmpi(command,'Copy')
        cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
        if iscell(cur_fline_selected)
          cur_fline_selected = cell2mat(cur_fline_selected);
        end
        if isempty(cur_fline_selected)
          warning('No flight line selected');
          return;
        end
        obj.selection.copy.flines = cur_fline_selected;
        obj.selection.copy.wpnts = [];
        
      elseif strcmpi(command,'Cut')
      elseif strcmpi(command,'Delete')
        deletePB_callback(obj,h_obj,event)
        
      elseif strcmpi(command,'Geometry Sort')
        cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
        if iscell(cur_fline_selected)
          cur_fline_selected = cell2mat(cur_fline_selected);
        end
        
        % Sort the selected flight line structures
        cur_pos = cur_fline_selected(1);
        sort_order = cur_fline_selected(1);
        used_mask = zeros(size(cur_fline_selected));
        used_mask(1) = 1;
        while length(sort_order) < length(cur_fline_selected)
          min_dist = inf;
          for pos = cur_fline_selected(~used_mask)
            % Search for next closest flight line (look at the start/stop
            % position of each line)
            dist_begin = sqrt(abs(obj.flines(pos).x(1) - obj.flines(cur_pos).x(end)).^2 ...
              + abs(obj.flines(pos).y(1) - obj.flines(cur_pos).y(end)).^2);
            dist_end = sqrt(abs(obj.flines(pos).x(end) - obj.flines(cur_pos).x(end)).^2 ...
              + abs(obj.flines(pos).y(end) - obj.flines(cur_pos).y(end)).^2);
            if dist_end < dist_begin && dist_end < min_dist
              % The end point of the next potential line is closer
              closest_pos = pos;
              min_dist = dist_end;
              min_is_end = true;
            elseif dist_begin < min_dist
              % The begin point of the next potential line is closer
              closest_pos = pos;
              min_dist = dist_begin;
              min_is_end = false;
            end
          end
          cur_pos = closest_pos;
          used_mask(closest_pos == cur_fline_selected) = 1;
          sort_order = cat(2,sort_order,closest_pos);
          if min_is_end
            % Reverse the waypoint order in the flight line since this
            % flight line will be flown from end to beginning
            obj.reverse_waypoints(cur_pos);
          end
        end
        
        % Update our "database" of flight lines
        obj.flines(cur_fline_selected) = obj.flines(sort_order);
        
        % Update the listbox
        flines = get(obj.h_gui.flines.listLB,'String');
        flines(cur_fline_selected) = flines(sort_order);
        set(obj.h_gui.flines.listLB,'String',flines);
        
      elseif strcmpi(command,'Merge')
        cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
        if iscell(cur_fline_selected)
          cur_fline_selected = cell2mat(cur_fline_selected);
        end
        if isempty(cur_fline_selected)
          warning('No flight line selected');
          return;
        end
        obj.selection.copy.flines = cur_fline_selected;
        
        prompt = {'Enter new flight line name:'};
        dlg_title = 'Flight line name';
        num_lines = 1;
        def = {'new_fline'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        
        if ~isempty(answer) && ~isempty(answer{1})
          answer = answer{1};
          fline.lat = [];
          fline.lon = [];
          fline.wpnt_names = {};
          
          for pos = obj.selection.copy.flines
            fline.lat = [fline.lat reshape(obj.flines(pos).lat,[1 length(obj.flines(pos).lat)])];
            fline.lon = [fline.lon reshape(obj.flines(pos).lon,[1 length(obj.flines(pos).lon)])];
            fline.wpnt_names = [fline.wpnt_names reshape(obj.flines(pos).wpnt_names,[1 length(obj.flines(pos).wpnt_names)])];
          end
          fline.name = answer;
          obj.insert_fline(fline,[],true);
          % Select the flight line that was just created
          set(obj.h_gui.flines.listLB,'Value',length(obj.flines));
        end
        
      elseif strcmpi(command,'Move')
        % Next Button Up: Get Position To Paste At (special mode = 10)
        % Next Button Up: End Position To Paste At (special mode = 11)
        obj.special_mode.type = 2;
        poly_handle = impoly;
        position = wait(poly_handle);
        [polyPts] = getPosition(poly_handle);
        xPoly = polyPts(:,1);
        yPoly = polyPts(:,2);
        delete(poly_handle);
        obj.special_mode.type = 0;
        
        cur_names = get(obj.h_gui.flines.listLB,'String');
        
        for idx = 2:length(xPoly)
          % At each vertex offset from the first entered by the user,
          % copy all the selected flight lines
          for pos = obj.selection.copy.flines
            x = obj.flines(pos).x + (xPoly(idx) - xPoly(1))*1e3;
            y = obj.flines(pos).y + (yPoly(idx) - yPoly(1))*1e3;
            
            [fline.lat,fline.lon] = projinv(obj.proj,x,y);
            for wpnt_idx = 1:length(fline.lat)
              fline.wpnt_names{wpnt_idx} = sprintf('M_%d', wpnt_idx);
            end
            fline.name = sprintf('%s_%d', obj.flines(pos).name, idx-1);
            obj.insert_fline(fline,[]);
          end
        end
        
      elseif strcmpi(command,'Paste')
        cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
        if iscell(cur_fline_selected)
          cur_fline_selected = cell2mat(cur_fline_selected);
        end
        if isempty(cur_fline_selected)
          warning('No flight line selected');
          return;
        end
        cur_fline_selected = cur_fline_selected(1);
        
        flines = obj.flines(obj.selection.copy.flines);
        for idx = 1:length(obj.selection.copy.flines)
          fline = flines(idx);
          if ~isempty(obj.selection.copy.wpnts)
            fline.lat = fline.lat(obj.selection.copy.wpnts);
            fline.lon = fline.lon(obj.selection.copy.wpnts);
            fline.wpnt_names = fline.wpnt_names(obj.selection.copy.wpnts);
          end
          obj.insert_fline(fline,cur_fline_selected+idx);
        end
        
      elseif strcmpi(command,'Paste Special')
        % Next Button Up: Get Position To Paste At (special mode = 10)
        % Next Button Up: End Position To Paste At (special mode = 11)
        obj.special_mode.type = 2;
        poly_handle = impoly;
        position = wait(poly_handle);
        [polyPts] = getPosition(poly_handle);
        xPoly = polyPts(:,1);
        yPoly = polyPts(:,2);
        delete(poly_handle);
        obj.special_mode.type = 0;
        
        cur_names = get(obj.h_gui.flines.listLB,'String');
        prompt = {'Aircraft Spacing (km)'};
        dlg_title = 'Paste Special';
        num_lines = 1;
        def = {'1'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        
        if ~isempty(answer) && ~isempty(answer{1})
          % Paste away (If "R" is fline, "R_0", "R_1", etc will be lines
          answer{1} = str2double(answer{1});
          
          [lat,lon] = projinv(obj.proj,xPoly*1e3,yPoly*1e3);
          [lat,lon] = interpm(lat,lon,1/100);
          along_track = geodetic_to_along_track(lat,lon);
          along_track_axis = 0:answer{1}*1e3:along_track(end);
          
          obj.selection.copy.flines = obj.selection.copy.flines(1);
          
          lat = interp1(along_track,lat,along_track_axis);
          lon = interp1(along_track,lon,along_track_axis);
          [x0,y0] = projfwd(obj.proj,lat,lon);
          for idx = 1:length(along_track_axis)
            x = obj.flines(obj.selection.copy.flines).x + x0(idx) - obj.flines(obj.selection.copy.flines).x(1);
            y = obj.flines(obj.selection.copy.flines).y + y0(idx) - obj.flines(obj.selection.copy.flines).y(1);
            [fline.lat,fline.lon] = projinv(obj.proj,x,y);
            for wpnt_idx = 1:length(fline.lat)
              fline.wpnt_names{wpnt_idx} = sprintf('S_%d', wpnt_idx);
            end
            fline.name = sprintf('%s_%d', obj.flines(obj.selection.copy.flines).name, idx);
            obj.insert_fline(fline,obj.selection.copy.flines+idx);
          end
        end
        
      elseif strcmpi(command,'Rename')
        cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
        if iscell(cur_fline_selected)
          cur_fline_selected = cell2mat(cur_fline_selected);
        end
        if isempty(cur_fline_selected)
          warning('No flight line selected');
          return;
        end
        cur_fline_selected = cur_fline_selected(1);
        
        cur_names = get(obj.h_gui.flines.listLB,'String');
        prompt = {'Rename flight line:'};
        dlg_title = 'Flight line rename';
        num_lines = 1;
        def = cur_names(cur_fline_selected);
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        
        if ~isempty(answer) && ~isempty(answer{1})
          answer = answer{1};
          cur_names(cur_fline_selected) = {answer};
          set(obj.h_gui.flines.listLB,'String',cur_names);
          obj.flines(cur_fline_selected).name = answer;
        end
        
      elseif strcmpi(command,'Reverse Flines')
        cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
        if iscell(cur_fline_selected)
          cur_fline_selected = cell2mat(cur_fline_selected);
        end
        if isempty(cur_fline_selected)
          warning('No flight line selected');
          return;
        end
        
        % Reverse the flight line order now
        obj.flines(cur_fline_selected) = obj.flines(cur_fline_selected(end:-1:1));
        
        obj.update_flineGraphics();
        
      elseif strcmpi(command,'Reverse Wpnts')
        cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
        obj.reverse_waypoints(cur_fline_selected);
        
      elseif strcmpi(command,'Rotate')
        cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
        obj.rotate_waypoints(cur_fline_selected, true);
      end
    end
    
    function rotate_waypoints(obj,cur_fline_selected,all_wpnts_en)
      cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
      if iscell(cur_fline_selected)
        cur_fline_selected = cell2mat(cur_fline_selected);
      end
      if isempty(cur_fline_selected)
        warning('No flight line selected');
        return;
      end
      cur_fline_selected = cur_fline_selected(1);
      
      if all_wpnts_en
        cur_wpnt_selected = 1:length(obj.flines(cur_fline_selected).x);
      else
        cur_wpnt_selected = get(obj.h_gui.wpnts.listLB,'Value');
        if iscell(cur_wpnt_selected)
          cur_wpnt_selected = cell2mat(cur_wpnt_selected);
        end
        if isempty(cur_wpnt_selected)
          warning('No waypoint selected');
          return;
        end
      end
      
      prompt = {'Rotation Angle (deg):'};
      def = {'90'};
      dlg_title = 'Rotate Waypoints';
      num_lines = 1;
      answer = inputdlg(prompt,dlg_title,num_lines,def);
      
      if length(answer) == 1 && ~isempty(answer{1})
        rot = eval(answer{1});
      end
      % Rotate around the first waypoint in the line
      x0 = obj.flines(cur_fline_selected).x(cur_wpnt_selected(1));
      y0 = obj.flines(cur_fline_selected).y(cur_wpnt_selected(1));
      for wpnt = cur_wpnt_selected(2:end)
        x = obj.flines(cur_fline_selected).x(wpnt) - x0;
        y = obj.flines(cur_fline_selected).y(wpnt) - y0;
        obj.flines(cur_fline_selected).x(wpnt) = x0 + cosd(rot)*x - sind(rot)*y;
        obj.flines(cur_fline_selected).y(wpnt) = y0 + sind(rot)*x + cosd(rot)*y;
      end
     [obj.flines(cur_fline_selected).lat obj.flines(cur_fline_selected).lon] ...
        = projinv(obj.proj,obj.flines(cur_fline_selected).x,obj.flines(cur_fline_selected).y);
      obj.flines(cur_fline_selected).along_track ...
        = geodetic_to_along_track(obj.flines(cur_fline_selected).lat, ...
        obj.flines(cur_fline_selected).lon);
      set(obj.flines(cur_fline_selected).handle,'XData',obj.flines(cur_fline_selected).x/1e3,'YData',obj.flines(cur_fline_selected).y/1e3);
      
      obj.wpntsLB_callback();
      obj.update_statusText();
    end
    
    function reverse_waypoints(obj,cur_fline_selected)
      if iscell(cur_fline_selected)
        cur_fline_selected = cell2mat(cur_fline_selected);
      end
      if isempty(cur_fline_selected)
        warning('No flight line selected');
        return;
      end
      
      % Reverse the waypoints
      for pos = cur_fline_selected
        obj.flines(pos).lat = obj.flines(pos).lat(end:-1:1);
        obj.flines(pos).lon = obj.flines(pos).lon(end:-1:1);
        obj.flines(pos).wpnt_names = obj.flines(pos).wpnt_names(end:-1:1);
        obj.flines(pos).x = obj.flines(pos).x(end:-1:1);
        obj.flines(pos).y = obj.flines(pos).y(end:-1:1);
        
        obj.flines(pos).along_track ...
          = geodetic_to_along_track(obj.flines(pos).lat, ...
          obj.flines(pos).lon);
        
        cur_wpnt_names = get(obj.h_gui.wpnts.listLB,'String');
        cur_wpnt_names = cur_wpnt_names(end:-1:1);
        set(obj.h_gui.wpnts.listLB,'String',cur_wpnt_names);
      end
      
      obj.update_statusText();
      obj.update_flineGraphics();
    end
    
    function delete_flines(obj,delete_mask)
      if isempty(delete_mask)
        return
      end
      
      for pos = 1:length(obj.flines)
        if delete_mask(pos)
          delete(obj.flines(pos).handle);
          delete(obj.flines(pos).start_handle);
        end
      end
      obj.flines = obj.flines(~delete_mask);
      
      % Update listbox strings and value/selection
      cur_fline_names = get(obj.h_gui.flines.listLB,'String');
      cur_fline_names = cur_fline_names(~delete_mask);
      set(obj.h_gui.flines.listLB,'Value',[]);
      set(obj.h_gui.flines.listLB,'String',cur_fline_names);
      
      set(obj.h_gui.wpnts.listLB,'ListboxTop',1);
      set(obj.h_gui.wpnts.listLB,'Value',[]);
      set(obj.h_gui.wpnts.listLB,'String',{});
    end
    
    function name_sortPB_callback(obj,h_obj,event)
      cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
      if iscell(cur_fline_selected)
        cur_fline_selected = cell2mat(cur_fline_selected);
      end
      if isempty(cur_fline_selected)
        warning('No flight line selected');
        return;
      end
      
      % Name sort the waypoints (i.e. order them alphebetically)
      for pos = cur_fline_selected
        cur_wpnt_selected = get(obj.h_gui.wpnts.listLB,'Value');
        if iscell(cur_wpnt_selected)
          cur_wpnt_selected = cell2mat(cur_wpnt_selected);
        end
        if isempty(cur_wpnt_selected)
          warning('No waypoint selected');
          return;
        end
        
        % Sort the selected waypoints by their names
        [tmp,sort_order] = sort(obj.flines(pos).wpnt_names(cur_wpnt_selected));
        
        obj.flines(pos).lat(cur_wpnt_selected) = obj.flines(pos).lat(cur_wpnt_selected(sort_order));
        obj.flines(pos).lon(cur_wpnt_selected) = obj.flines(pos).lon(cur_wpnt_selected(sort_order));
        obj.flines(pos).wpnt_names(cur_wpnt_selected) = obj.flines(pos).wpnt_names(cur_wpnt_selected(sort_order));
        obj.flines(pos).x(cur_wpnt_selected) = obj.flines(pos).x(cur_wpnt_selected(sort_order));
        obj.flines(pos).y(cur_wpnt_selected) = obj.flines(pos).y(cur_wpnt_selected(sort_order));
        
        obj.flines(pos).along_track ...
          = geodetic_to_along_track(obj.flines(pos).lat, ...
          obj.flines(pos).lon);
        
        cur_wpnt_names = get(obj.h_gui.wpnts.listLB,'String');
        cur_wpnt_names = cur_wpnt_names(end:-1:1);
        set(obj.h_gui.wpnts.listLB,'String',cur_wpnt_names);
      end
      
      obj.update_statusText();
      obj.update_flineGraphics();
    end
    
    function geo_sortPB_callback(obj,h_obj,event)
      cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
      if iscell(cur_fline_selected)
        cur_fline_selected = cell2mat(cur_fline_selected);
      end
      if isempty(cur_fline_selected)
        warning('No flight line selected');
        return;
      end
      
      % Geometry sort the waypoints (i.e. order them starting from the
      % first according to the distance between them)
      for pos = cur_fline_selected
        cur_wpnt_selected = get(obj.h_gui.wpnts.listLB,'Value');
        if iscell(cur_wpnt_selected)
          cur_wpnt_selected = cell2mat(cur_wpnt_selected);
        end
        if isempty(cur_wpnt_selected)
          warning('No waypoint selected');
          return;
        end
        
        % Sort the selected waypoints
        cur_pos = cur_wpnt_selected(1);
        sort_order = cur_wpnt_selected(1);
        used_mask = zeros(size(cur_wpnt_selected));
        used_mask(1) = 1;
        while length(sort_order) < length(cur_wpnt_selected)
          min_dist = inf;
          for next_pos = cur_wpnt_selected(~used_mask)
            % Search for next closest waypoint
            dist = sqrt(abs(obj.flines(pos).x(next_pos) - obj.flines(pos).x(cur_pos)).^2 ...
              + abs(obj.flines(pos).y(next_pos) - obj.flines(pos).y(cur_pos)).^2);
            if dist < min_dist
              % The begin point of the next potential line is closer
              closest_pos = next_pos;
              min_dist = dist;
            end
          end
          cur_pos = closest_pos;
          used_mask(closest_pos == cur_wpnt_selected) = 1;
          sort_order = cat(2,sort_order,closest_pos);
        end
        
        obj.flines(pos).lat(cur_wpnt_selected) = obj.flines(pos).lat(cur_wpnt_selected(sort_order));
        obj.flines(pos).lon(cur_wpnt_selected) = obj.flines(pos).lon(cur_wpnt_selected(sort_order));
        obj.flines(pos).wpnt_names(cur_wpnt_selected) = obj.flines(pos).wpnt_names(cur_wpnt_selected(sort_order));
        obj.flines(pos).x(cur_wpnt_selected) = obj.flines(pos).x(cur_wpnt_selected(sort_order));
        obj.flines(pos).y(cur_wpnt_selected) = obj.flines(pos).y(cur_wpnt_selected(sort_order));
        
        obj.flines(pos).along_track ...
          = geodetic_to_along_track(obj.flines(pos).lat, ...
          obj.flines(pos).lon);
        
        cur_wpnt_names = get(obj.h_gui.wpnts.listLB,'String');
        cur_wpnt_names = cur_wpnt_names(end:-1:1);
        set(obj.h_gui.wpnts.listLB,'String',cur_wpnt_names);
      end
      
      obj.update_statusText();
      obj.update_flineGraphics();
    end
    
    function deletePB_callback(obj,h_obj,event)
      cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
      if iscell(cur_fline_selected)
        cur_fline_selected = cell2mat(cur_fline_selected);
      end
      
      % Delete flight line structures
      delete_mask = logical(zeros(size(obj.flines)));
      delete_mask(cur_fline_selected) = 1;
      obj.delete_flines(delete_mask);
    end
    
    function delete_pointPB_callback(obj,h_obj,event)
      cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
      if iscell(cur_fline_selected)
        cur_fline_selected = cell2mat(cur_fline_selected);
      end
      if isempty(cur_fline_selected)
        warning('Can not insert: no flight line selected\n');
        return;
      end
      cur_fline_selected = cur_fline_selected(1);
      
      wpnts = get(obj.h_gui.wpnts.listLB,'String');
      
      cur_wpnt_selected = get(obj.h_gui.wpnts.listLB,'Value');
      if iscell(cur_wpnt_selected)
        cur_wpnt_selected = cell2mat(cur_wpnt_selected);
      end
      delete_mask = logical(zeros(size(obj.flines(cur_fline_selected).x)));
      delete_mask(cur_wpnt_selected) = 1;
      obj.flines(cur_fline_selected).x = obj.flines(cur_fline_selected).x(~delete_mask);
      obj.flines(cur_fline_selected).y = obj.flines(cur_fline_selected).y(~delete_mask);
      obj.flines(cur_fline_selected).wpnt_names = obj.flines(cur_fline_selected).wpnt_names(~delete_mask);
      
      [obj.flines(cur_fline_selected).lat obj.flines(cur_fline_selected).lon] ...
        = projinv(obj.proj,obj.flines(cur_fline_selected).x,obj.flines(cur_fline_selected).y);
      obj.flines(cur_fline_selected).along_track ...
        = geodetic_to_along_track(obj.flines(cur_fline_selected).lat, ...
        obj.flines(cur_fline_selected).lon);
      
      wpnts = wpnts(~delete_mask);
      set(obj.h_gui.wpnts.listLB,'Value',[]);
      set(obj.h_gui.wpnts.listLB,'String',wpnts);
      
      set(obj.flines(cur_fline_selected).handle,'XData',obj.flines(cur_fline_selected).x/1e3,'YData',obj.flines(cur_fline_selected).y/1e3);
      
      obj.update_statusText();
      obj.wpntsLB_callback();
    end
    
    function toolPM_callback(obj,h_obj,event)
      
    end
    
    function exportPB_callback(obj,h_obj,event)
      % Try to create a default file path that makes sense
      [save_fn_dir,save_fn_name,save_fn_ext] = fileparts(obj.save_fn);
      if exist(save_fn_dir,'dir')
        export_fn_dir = save_fn_dir;
        if strcmpi(save_fn_ext,'.mat')
          export_fn_name = save_fn_name;
          export_fn_ext = '.jpg';
        else
          export_fn_name = 'export'
          export_fn_ext = '.jpg';
        end
        obj.export_fn = fullfile(export_fn_dir,[export_fn_name export_fn_ext]);
      end

      % Ask the user for the file path to use
      [filename, pathname] = uiputfile( ...
        {'*.jpg', 'All JPG Images (*.jpg)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Save as', obj.export_fn);
      
      if isequal(filename,0) || isequal(pathname,0)
        return;
      end
      
      obj.export_fn = fullfile(pathname, filename);
      fprintf('Exporting to %s\n',obj.export_fn);
      
      %% Create a jpg of map showing what ever is currently in view in the browser
      h_export_fig = figure;
      h_export_axes = axes;
      for child = get(obj.h_axes,'Children')'
        if strcmpi(get(child,'Type'),'Image')
          h_image = child;
          break;
        end
      end
      XData = get(h_image,'XData');
      YData = get(h_image,'YData');
      CData = get(h_image,'CData');
      h_export_image = imagesc(XData,YData,CData,'Parent',h_export_axes);
      xlabel('X (km)');
      ylabel('Y (km)');
      set(h_export_axes,'YDir','normal');
      axis(h_export_axes, axis(obj.h_axes));
      
      hold(h_export_axes,'on')
      for pos = 1:length(obj.flines)
        plot(obj.flines(pos).x/1e3,obj.flines(pos).y/1e3,'b.-','Parent',h_export_axes);
        if ~isempty(obj.flines(pos).x)
          plot(obj.flines(pos).x(1)/1e3,obj.flines(pos).y(1)/1e3,'bo','Parent',h_export_axes);
        end
        for wpnt_idx = 1:length(obj.flines(pos).x)
          text(obj.flines(pos).x(wpnt_idx)/1e3,obj.flines(pos).y(wpnt_idx)/1e3, ...
            obj.flines(pos).wpnt_names{wpnt_idx}, 'Color','black', ...
            'Parent',h_export_axes,'Interpreter','none','FontSize',10,'BackgroundColor','white','margin',0.5);
        end
      end
      hold(h_export_axes,'off')
      
      hold(h_export_axes,'on')
      plotonly = handle2struct(obj.plotonly.handles);
      struct2handle(plotonly,h_export_axes);
      hold(h_export_axes,'off')
      
      set(h_export_fig,'PaperPosition',[0.5 0.5 10 7.5]);
      set(h_export_fig,'PaperOrientation','Portrait');
      print(h_export_fig,'-djpeg','-r200',obj.export_fn);
      
      %% Create a Matlab figure file
      [export_fn_dir export_fn_name] = fileparts(obj.export_fn);
      export_fig_fn = fullfile(export_fn_dir,[export_fn_name '.fig']);
      clip_and_resample_image(h_export_image,h_export_axes,8);
      saveas(h_export_fig, export_fig_fn);
      fprintf('  %s\n', export_fig_fn);
      
      %% Create a separate CSV file for each flight line
      for pos = 1:length(obj.flines)
        export_csv_fn = fullfile(export_fn_dir,sprintf('%s_%s.csv', ...
          export_fn_name, obj.flines(pos).name));
        fprintf('  %s\n', export_csv_fn);
        
        [fid,msg] = fopen(export_csv_fn,'w');
        if fid < 0
          error('Could not open %s for writing: %s', export_csv_fn, msg);
        end
        
        fprintf(fid,'%15s,%15s,%15s,%15s\n','Lat_North','Lon_East','Name','Distance_km');
        along_track = geodetic_to_along_track(obj.flines(pos).lat,obj.flines(pos).lon,zeros(size(obj.flines(pos).lat)));
        for wpnt_idx = 1:length(obj.flines(pos).x)
          fprintf(fid,'%15.6f,%15.6f,%15s,%15.0f\n',obj.flines(pos).lat(wpnt_idx), ...
            obj.flines(pos).lon(wpnt_idx), obj.flines(pos).wpnt_names{wpnt_idx}, ...
            along_track(wpnt_idx)/1e3);
        end
        fclose(fid);
      end
        
      %% Create a separate TXT file for each flight line (AWI Format)
      for pos = 1:length(obj.flines)
        export_txt_fn = fullfile(export_fn_dir,sprintf('%s_AWI_%s.txt', ...
          export_fn_name, obj.flines(pos).name));
        fprintf('  %s\n', export_txt_fn);
        
        [fid,msg] = fopen(export_txt_fn,'w');
        if fid < 0
          error('Could not open %s for writing: %s', export_txt_fn, msg);
        end
        
        %% Create elevation profile
        elev = latlon2elevation_profile(obj.flines(pos).lat,obj.flines(pos).lon,true,1/100,1) + 500;
        
        fprintf(fid,'[Settings]\r\nRadiusEarthKm=6371.0\r\n[Track]\r\n;TrackPoints[Label]	Longitude[°]Latitude[°]	Altitude[m]\r\n');
        along_track = geodetic_to_along_track(obj.flines(pos).lat,obj.flines(pos).lon,zeros(size(obj.flines(pos).lat)));
        for wpnt_idx = 1:length(obj.flines(pos).x)
          fprintf(fid,'%s\t%.5f\t%.5f\t%.0f', obj.flines(pos).wpnt_names{wpnt_idx}, ...
            obj.flines(pos).lon(wpnt_idx),obj.flines(pos).lat(wpnt_idx),  ...
            elev(wpnt_idx));
          if wpnt_idx < length(obj.flines(pos).x)
            fprintf(fid,'\r\n');
          end
        end
        fclose(fid);
      end
      
      %% Create a separate CSV file for each flight line (Kenn Borek Air Format)
      for pos = 1:length(obj.flines)
        export_csv_fn = fullfile(export_fn_dir,sprintf('%s_KBA_%s.csv', ...
          export_fn_name, obj.flines(pos).name));
        fprintf('  %s\n', export_csv_fn);
        
        [fid,msg] = fopen(export_csv_fn,'w');
        if fid < 0
          error('Could not open %s for writing: %s', export_csv_fn, msg);
        end
        fprintf(fid,'%s,%s,%s,%s,%s,%s\n', ...
          'Name','Lat_North_deg','Lat_North_min', ...
          'Lon_East_deg','Lon_East_min','Distance_km');
        along_track = geodetic_to_along_track(obj.flines(pos).lat,obj.flines(pos).lon,zeros(size(obj.flines(pos).lat)));
        for wpnt_idx = 1:length(obj.flines(pos).x)
          lat_deg = fix(obj.flines(pos).lat(wpnt_idx));
          lat_min = abs((obj.flines(pos).lat(wpnt_idx) - lat_deg)*60);
          lon_deg = fix(obj.flines(pos).lon(wpnt_idx));
          lon_min = abs((obj.flines(pos).lon(wpnt_idx) - lon_deg)*60);
          fprintf(fid,'%s,%15.0f,%15.2f,%15.0f,%15.2f,%15.0f\n', ...
            obj.flines(pos).wpnt_names{wpnt_idx}, lat_deg, lat_min, ...
            lon_deg, lon_min, along_track(wpnt_idx)/1e3);
        end
        fclose(fid);
      end
      
      %% Create a KML file
      export_kml_fn = fullfile(export_fn_dir,sprintf('%s.kml', ...
        export_fn_name));
      fprintf('  %s\n', export_kml_fn);
      p = geoshape;
      name = {};
      color = {};
      for pos = 1:length(obj.flines)
        if ~isempty(obj.flines(pos).lat)
          name{pos} = obj.flines(pos).name;
          color{pos} = 'Blue';
          p = append(p, obj.flines(pos).lat, obj.flines(pos).lon, 'Name', obj.flines(pos).name);
        end
      end
      
      description = cell(size(name));
      kmlwrite(export_kml_fn, p, 'Color', color, 'Width', 2, ...
        'Description', description, 'Name', name);
      
      %% Create a Garmin Flight Plan file
      export_fpl_fn = fullfile(export_fn_dir,sprintf('%s.fpl', ...
        export_fn_name));
      for pos = 1:length(obj.flines)
        if ~isempty(obj.flines(pos).lat)
          fprintf('  %s\n', export_fpl_fn);
          write_garmin_fpl(export_fpl_fn, obj.flines(pos));
        end
      end
      
    end
    
    function insert_beforePB_callback(obj,h_obj,event)
      
    end
    
    function moveupPB_callback(obj,h_obj,event)
      obj.reset_copy();
      if obj.h_gui.flines.moveupPB == h_obj
        cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
        if iscell(cur_fline_selected)
          cur_fline_selected = cell2mat(cur_fline_selected);
        end
        if isempty(cur_fline_selected)
          warning('Can not moveup: no flight line selected\n');
          return;
        end
        cur_fline_selected = cur_fline_selected(1);
        
        flines = get(obj.h_gui.flines.listLB,'String');
        
        if cur_fline_selected > 1
          % Update our "database" of flight lines
          obj.flines(cur_fline_selected-1:cur_fline_selected) = obj.flines(cur_fline_selected:-1:cur_fline_selected-1);
          
          % Update the listbox
          flines(cur_fline_selected-1:cur_fline_selected) = flines(cur_fline_selected:-1:cur_fline_selected-1);
          set(obj.h_gui.flines.listLB,'String',flines);
          set(obj.h_gui.flines.listLB,'Value',cur_fline_selected-1);
        end
        
      elseif obj.h_gui.wpnts.moveupPB == h_obj
        cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
        if iscell(cur_fline_selected)
          cur_fline_selected = cell2mat(cur_fline_selected);
        end
        if isempty(cur_fline_selected)
          warning('Can not moveup: no flight line selected\n');
          return;
        end
        cur_fline_selected = cur_fline_selected(1);
        
        wpnts = get(obj.h_gui.wpnts.listLB,'String');
        
        cur_wpnt_selected = get(obj.h_gui.wpnts.listLB,'Value');
        if iscell(cur_wpnt_selected)
          cur_wpnt_selected = cell2mat(cur_wpnt_selected);
        end
        cur_wpnt_selected = cur_wpnt_selected(1);
        
        if cur_wpnt_selected > 1
          % Update our "database" of flight lines
          obj.flines(cur_fline_selected).x(cur_wpnt_selected-1:cur_wpnt_selected) = obj.flines(cur_fline_selected).x(cur_wpnt_selected:-1:cur_wpnt_selected-1);
          obj.flines(cur_fline_selected).y(cur_wpnt_selected-1:cur_wpnt_selected) = obj.flines(cur_fline_selected).y(cur_wpnt_selected:-1:cur_wpnt_selected-1);
          obj.flines(cur_fline_selected).wpnt_names(cur_wpnt_selected-1:cur_wpnt_selected) = obj.flines(cur_fline_selected).wpnt_names(cur_wpnt_selected:-1:cur_wpnt_selected-1);
          
          [obj.flines(cur_fline_selected).lat obj.flines(cur_fline_selected).lon] ...
            = projinv(obj.proj,obj.flines(cur_fline_selected).x,obj.flines(cur_fline_selected).y);
          
          % Update the listbox
          wpnts(cur_wpnt_selected-1:cur_wpnt_selected) = wpnts(cur_wpnt_selected:-1:cur_wpnt_selected-1);
          set(obj.h_gui.wpnts.listLB,'String',wpnts);
          set(obj.h_gui.wpnts.listLB,'Value',cur_wpnt_selected-1);
          
          % Update the map
          set(obj.flines(cur_fline_selected).handle,'XData',obj.flines(cur_fline_selected).x/1e3,'YData',obj.flines(cur_fline_selected).y/1e3);
        end
      end
    end
    
    function movedownPB_callback(obj,h_obj,event)
      obj.reset_copy();
      if obj.h_gui.flines.movedownPB == h_obj
        cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
        if iscell(cur_fline_selected)
          cur_fline_selected = cell2mat(cur_fline_selected);
        end
        if isempty(cur_fline_selected)
          warning('Can not moveup: no flight line selected\n');
          return;
        end
        cur_fline_selected = cur_fline_selected(1);
        
        flines = get(obj.h_gui.flines.listLB,'String');
        
        if cur_fline_selected < length(flines)
          % Update our "database" of flight lines
          obj.flines(cur_fline_selected:cur_fline_selected+1) = obj.flines(cur_fline_selected+1:-1:cur_fline_selected);
          
          % Update the listbox
          flines(cur_fline_selected:cur_fline_selected+1) = flines(cur_fline_selected+1:-1:cur_fline_selected);
          set(obj.h_gui.flines.listLB,'String',flines);
          set(obj.h_gui.flines.listLB,'Value',cur_fline_selected+1);
        end
        
      elseif obj.h_gui.wpnts.movedownPB == h_obj
        cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
        if iscell(cur_fline_selected)
          cur_fline_selected = cell2mat(cur_fline_selected);
        end
        if isempty(cur_fline_selected)
          warning('Can not moveup: no flight line selected\n');
          return;
        end
        cur_fline_selected = cur_fline_selected(1);
        
        wpnts = get(obj.h_gui.wpnts.listLB,'String');
        
        cur_wpnt_selected = get(obj.h_gui.wpnts.listLB,'Value');
        if iscell(cur_wpnt_selected)
          cur_wpnt_selected = cell2mat(cur_wpnt_selected);
        end
        cur_wpnt_selected = cur_wpnt_selected(1);
        
        if cur_wpnt_selected < length(wpnts)
          % Update our "database" of flight lines
          obj.flines(cur_fline_selected).x(cur_wpnt_selected:cur_wpnt_selected+1) ...
            = obj.flines(cur_fline_selected).x(cur_wpnt_selected+1:-1:cur_wpnt_selected);
          obj.flines(cur_fline_selected).y(cur_wpnt_selected:cur_wpnt_selected+1) ...
            = obj.flines(cur_fline_selected).y(cur_wpnt_selected+1:-1:cur_wpnt_selected);
          obj.flines(cur_fline_selected).wpnt_names(cur_wpnt_selected:cur_wpnt_selected+1) ...
            = obj.flines(cur_fline_selected).wpnt_names(cur_wpnt_selected+1:-1:cur_wpnt_selected);
          
          [obj.flines(cur_fline_selected).lat obj.flines(cur_fline_selected).lon] ...
            = projinv(obj.proj,obj.flines(cur_fline_selected).x,obj.flines(cur_fline_selected).y);
          
          % Update the listbox
          wpnts(cur_wpnt_selected:cur_wpnt_selected+1) = wpnts(cur_wpnt_selected+1:-1:cur_wpnt_selected);
          set(obj.h_gui.wpnts.listLB,'String',wpnts);
          set(obj.h_gui.wpnts.listLB,'Value',cur_wpnt_selected+1);
          
          % Update the map
          set(obj.flines(cur_fline_selected).handle,'XData',obj.flines(cur_fline_selected).x/1e3,'YData',obj.flines(cur_fline_selected).y/1e3);
        end
      end
    end
    
    function insert_afterPB_callback(obj,h_obj,event)
      if ~obj.special_mode.type
        cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
        if iscell(cur_fline_selected)
          cur_fline_selected = cell2mat(cur_fline_selected);
        end
        if isempty(cur_fline_selected)
          warning('Can not insert: no flight line selected\n');
          return;
        end
        
        % Get current selection
        obj.special_mode.flines = get(obj.h_gui.flines.listLB,'Value');
        obj.special_mode.wpnts = get(obj.h_gui.wpnts.listLB,'Value');
        % Deselect all
        set(obj.h_gui.flines.listLB,'Value',[]);
        set(obj.h_gui.wpnts.listLB,'Value',[]);
        % Enter special selection mode
        obj.special_mode.type = true;
        
      else
        % Click insert after again completes action
        obj.reset_copy();
        obj.special_mode.type = false;
        
        % Insert all waypoints that were selected
        cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
        if iscell(cur_fline_selected)
          cur_fline_selected = cell2mat(cur_fline_selected);
        end
        cur_fline_selected = cur_fline_selected(1);
        
        cur_wpnt_selected = get(obj.h_gui.wpnts.listLB,'Value');
        if iscell(cur_wpnt_selected)
          cur_wpnt_selected = cell2mat(cur_wpnt_selected);
        end
        
        set(obj.h_gui.flines.listLB,'Value',obj.special_mode.flines);
        obj.flinesLB_callback();
        set(obj.h_gui.wpnts.listLB,'Value',obj.special_mode.wpnts);
        obj.wpntsLB_callback();
        if isempty(cur_wpnt_selected)
          warning('No points selected to insert\n');
          return;
        else
          x = obj.flines(cur_fline_selected).x(cur_wpnt_selected);
          y = obj.flines(cur_fline_selected).y(cur_wpnt_selected);
          name = obj.flines(cur_fline_selected).wpnt_names(cur_wpnt_selected);
          for idx = 1:length(x)
            obj.insert_wpnt(x(idx)/1e3,y(idx)/1e3,name{idx});
          end
        end
      end
      obj.wpntsLB_callback();
      obj.update_statusText();
      
    end
    
    function update_statusText(obj)
      cur_fline_selected = get(obj.h_gui.flines.listLB,'Value');
      if iscell(cur_fline_selected)
        cur_fline_selected = cell2mat(cur_fline_selected);
      end
      if isempty(cur_fline_selected)
        % Empty way points listbox
        set(obj.h_gui.wpnts.listLB,'String',{});
        set(obj.h_gui.wpnts.listLB,'Value',[]);
        return;
      end
      
      total_length = 0;
      for pos = 1:length(obj.flines)
        if any(pos == cur_fline_selected)
          if ~isempty(obj.flines(pos).along_track)
            total_length = total_length + obj.flines(pos).along_track(end);
          end
        end
      end
      set(obj.h_gui.statusText,'String',sprintf('%.0f km\n', total_length/1e3));
    end
    
  end
end



