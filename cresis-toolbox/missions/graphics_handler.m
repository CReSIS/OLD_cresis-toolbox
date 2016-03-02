% Class graphics_handler
%
% Class which handles graphics on an axes. Provides a figure window with a
% listbox for selecting and deleting graphics.
%
% h_fig = figure; clf;
% h_axes = axes;
% h_plot = plot(1:10,'parent',h_axes);
% obj = graphics_handler(h_axes,'delete');
% obj.insert_handle(h_plot,'Test');
% hold(h_axes,'on');
% h_plot = plot(3:-1:1,'r','parent',h_axes);
% obj.insert_handle(h_plot,'Test2');
% obj.delete_handle('Test');

classdef (HandleCompatible = true) graphics_handler < handle
  properties
    handles % list of graphics handles (1xN array of handles)
    handle_names % list of graphics handle names (1xN cell array of strings)
    h_axes
    
    cur_handle_selected
    attribute_list_original
    close_action
    
    h_fig % Figure handle
    h_gui % Structure of graphics handles
  end
  
  events
    hide_figure
  end
  
  methods
    function obj = graphics_handler(h_axes,close_action)
      obj.h_axes = h_axes;
      obj.handles = [];
      obj.handle_names = {};
      obj.cur_handle_selected = [];
      obj.attribute_list_original = {};
      obj.close_action = close_action; % 'delete','hide','none'
      
      obj.h_fig = figure;
      set(obj.h_fig,'MenuBar','none');
      set(obj.h_fig,'ToolBar','none');
      set(obj.h_fig,'DockControls','off');
      h_fig_pos = get(obj.h_fig,'Position');
      set(obj.h_fig,'Position',[h_fig_pos(1:2) 200 200]);
      
      % h_gui.h_panel
      % h_gui.h_axes
      % h_gui.h_table
      % h_gui.h_panel_table
      % h_gui.h_flines
      %% Create widgets of main table
      
      obj.h_gui.handleText = uicontrol('parent',obj.h_fig);
      set(obj.h_gui.handleText,'Style','text');
      set(obj.h_gui.handleText,'HorizontalAlignment','left');
      set(obj.h_gui.handleText,'String','Graphics handle list');
      set(obj.h_gui.handleText,'Position',[1 1 120 15]);
      set(obj.h_gui.handleText,'Units','points');
      
      obj.h_gui.handlesLB = uicontrol('Parent',obj.h_fig);
      set(obj.h_gui.handlesLB,'Style','listbox');
      set(obj.h_gui.handlesLB,'HorizontalAlignment','Center');
      set(obj.h_gui.handlesLB,'FontName','fixed');
      set(obj.h_gui.handlesLB,'Value',[]);
      set(obj.h_gui.handlesLB,'Callback',@obj.handlesLB_callback);
      set(obj.h_gui.handlesLB,'Max',1e9);
      
            %% Source list box context menu
      obj.h_gui.handlesCM = uicontextmenu;
      % Define the context menu items and install their callbacks
      obj.h_gui.handlesCM_item1 = uimenu(obj.h_gui.handlesCM, 'Label', 'Delete', 'Callback', @obj.handlesLB_menu_callback);
      set(obj.h_gui.handlesLB,'uicontextmenu',obj.h_gui.handlesCM)
      
      obj.h_gui.selectText = uicontrol('parent',obj.h_fig);
      set(obj.h_gui.selectText,'Style','text');
      set(obj.h_gui.selectText,'HorizontalAlignment','left');
      set(obj.h_gui.selectText,'String','Selection properties');
      set(obj.h_gui.selectText,'Position',[1 1 120 15]);
      set(obj.h_gui.selectText,'Units','points');
      
      obj.h_gui.selectLE = uicontrol('parent',obj.h_fig);
      set(obj.h_gui.selectLE,'Style','edit');
      set(obj.h_gui.selectLE,'HorizontalAlignment','left');
      set(obj.h_gui.selectLE,'String','{''LineWidth'',2,''MarkerSize'',20}');
      set(obj.h_gui.selectLE,'Position',[1 1 120 20]);
      set(obj.h_gui.selectLE,'Units','points');
      
      %% Setup main table
      obj.h_gui.h_table.ui = obj.h_fig;
      obj.h_gui.h_table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
      obj.h_gui.h_table.height_margin = NaN*zeros(30,30);
      obj.h_gui.h_table.false_height = NaN*zeros(30,30);
      row = 1; col = 1;
      obj.h_gui.h_table.handles{row,col}   = obj.h_gui.handleText;
      obj.h_gui.h_table.width(row,col)     = inf;
      obj.h_gui.h_table.height(row,col)    = 15;
      obj.h_gui.h_table.false_height(row,col) = 0;
      row = row + 1;
      obj.h_gui.h_table.handles{row,col}   = obj.h_gui.handlesLB;
      obj.h_gui.h_table.width(row,col)     = inf;
      obj.h_gui.h_table.height(row,col)    = inf;
      obj.h_gui.h_table.false_height(row,col) = 0;
      row = row + 1;
      obj.h_gui.h_table.handles{row,col}   = obj.h_gui.selectText;
      obj.h_gui.h_table.width(row,col)     = inf;
      obj.h_gui.h_table.height(row,col)    = 15;
      obj.h_gui.h_table.false_height(row,col) = 0;
      row = row + 1;
      obj.h_gui.h_table.handles{row,col}   = obj.h_gui.selectLE;
      obj.h_gui.h_table.width(row,col)     = inf;
      obj.h_gui.h_table.height(row,col)    = 20;
      obj.h_gui.h_table.false_height(row,col) = 0;
      
      obj.h_gui.h_table.width_margin ...
        = obj.h_gui.h_table.width_margin(1:row,1:col);
      obj.h_gui.h_table.height_margin ...
        = obj.h_gui.h_table.height_margin(1:row,1:col);
      obj.h_gui.h_table.false_height ...
        = obj.h_gui.h_table.false_height(1:row,1:col);
      clear row col
      table_draw(obj.h_gui.h_table);
      
      %% Set up general handles
      set(obj.h_fig,'CloseRequestFcn',@obj.close_win);
      
    end
    
    function delete(obj)
      % Delete the figure handle
      try
        delete(obj.h_fig);
      end
      % Delete the children
      for idx = 1:length(obj.handles)
        try
          delete(obj.handles(idx));
        end
      end
    end
    
    function close_win(obj,h_obj,event)
      if strcmpi(obj.close_action,'delete')
        delete(obj);
      elseif strcmpi(obj.close_action,'hide')
        set_visibility(obj,'off');
        notify(obj,'hide_figure');
      end
    end
    
    function set_visibility(obj,state)
      set(obj.h_fig,'Visible',state);
    end
    
    function delete_handle(obj,handle_refs)
      delete_match_mask = zeros(size(obj.handles));
      selection_mask = zeros(size(obj.handles));
      for handle_idx = 1:length(handle_refs)
        if iscell(handle_refs)
          handle_ref = handle_refs{handle_idx};
        else
          handle_ref = handle_refs(handle_idx);
        end
        % handle_ref: name of handle (string) or index (integer) of handle to delete
        if ischar(handle_ref)
          match_idx = find(strcmpi(handle_ref,obj.handle_names));
        else
          match_idx = handle_ref;
        end
        if ~isempty(match_idx) && match_idx >= 1 && match_idx <= length(obj.handles);
          try
            delete(obj.handles(match_idx));
          end
          delete_match_mask(match_idx) = 1;
        end
      end
      selection_idxs = get(obj.h_gui.handlesLB,'Value');
      selection_mask(selection_idxs) = 1;
      selection_mask = selection_mask & delete_match_mask;
      selection_mask = selection_mask(~delete_match_mask);
      obj.handles = obj.handles(~delete_match_mask);
      obj.handle_names = obj.handle_names(~delete_match_mask);
      set(obj.h_gui.handlesLB,'String',obj.handle_names,'Value',find(selection_mask),'ListboxTop',1);
    end
    
    function insert_handle(obj,handles,names)
      for handle_idx = 1:length(handles)
        handle = handles(handle_idx);
        if exist('names','var')
          name = names{handle_idx};
        else
          name = '';
        end
        if ishandle(handle)
          obj.handles = [obj.handles handle];
          if isempty(name)
            name = sprintf('%d', length(obj.handles));
          end
          obj.handle_names = [obj.handle_names name];
          set(obj.h_gui.handlesLB,'String',obj.handle_names);
        else
          warning('Invalid handle, not added to list.');
        end
      end
    end
    
    function handlesLB_callback(obj,h_obj,event)
      cur_handle_selected = get(obj.h_gui.handlesLB,'Value');
      if iscell(cur_handle_selected)
        cur_handle_selected = cell2mat(cur_handle_selected);
      end
      if ~isempty(obj.cur_handle_selected)
        try
          for idx = 1:2:length(obj.attribute_list_original)
            set(obj.handles(obj.cur_handle_selected), ...
              obj.attribute_list_original{idx},obj.attribute_list_original{idx+1});
          end
        end
      end
      if isempty(cur_handle_selected)
        return;
      end
      % Only highlight the first entry
      cur_handle_selected = cur_handle_selected(1);
      try
        attribute_list = eval(get(obj.h_gui.selectLE,'String'));
        obj.cur_handle_selected = cur_handle_selected;
        for idx = 1:2:length(attribute_list)
          obj.attribute_list_original{idx} = attribute_list{idx};
          obj.attribute_list_original{idx+1} = get(obj.handles(cur_handle_selected),attribute_list{idx});
          set(obj.handles(cur_handle_selected),attribute_list{idx},attribute_list{idx+1});
        end
      end
    end
    
    function handlesLB_menu_callback(obj,h_obj,event)
      cur_handle_selected = get(obj.h_gui.handlesLB,'Value');
      if iscell(cur_handle_selected)
        cur_handle_selected = cell2mat(cur_handle_selected);
      end
      delete_handle(obj,cur_handle_selected);
    end
  
  end
end



