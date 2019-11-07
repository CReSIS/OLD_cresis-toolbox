classdef popupmenu_edit < handle
% popupmenu_edit: A class for creating a popup menu which allows items to
% be edited. Press enter to add a new item. Press backspace to edit an
% item. Pressing enter adds the item. If the string is empty then no item
% is added or the old item is deleted.
%
% parent: parent GUI handle for each of the GUI objects
% list_values: cell vector of strings to populate popup menu with, default
%   is empty and popup menu starts as an edit box to add an item
% value: scalar integer index into list_values vector
%
% See also: popupmenu_edit, run_popupmenu_edit
%
% Author: John Paden

  properties
    h_valuePM
    h_valueCM
    value
    list_values % cell array of strings to show in popup menu
    enable_mode % true to enable
    style_mode % 'edit' or 'popupmenu'
  end
  
  properties (SetAccess = private, GetAccess = private)
  end
   
  events
    selection_changed
  end
  
  methods
    function obj = popupmenu_edit(parent,list_values,value)
      % List boxes start enabled by default
      obj.enable_mode = true;
      if exist('list_values','var')
        obj.list_values = list_values;
      else
        obj.list_values = {};
      end
      if exist('value','var')
        obj.value = value;
      else
        obj.value = [];
      end
      if isempty(obj.list_values)
        obj.style_mode = 0;
        obj.value = [];
      else
        obj.style_mode = 1;
        if isempty(obj.value)
          obj.value = 1;
        end
        obj.value = round(obj.value);
        if isempty(obj.value) || obj.value < 1 || obj.value > length(obj.list_values)
          obj.value = [];
        end
      end
      
      obj.h_valuePM = uicontrol('Parent',parent);
      if obj.style_mode
        set(obj.h_valuePM,'Style','popupmenu');
      else
        set(obj.h_valuePM,'Style','edit');
      end
      set(obj.h_valuePM,'String',obj.list_values);
      set(obj.h_valuePM,'Value',obj.value);
      set(obj.h_valuePM,'HorizontalAlignment','Center');
      set(obj.h_valuePM,'Callback',@obj.valuePM_callback);
      set(obj.h_valuePM,'KeyPressFcn',@obj.key_press);

      obj.h_valueCM = uicontextmenu('Parent',parent);
      % Define the context menu items and install their callbacks
      uimenu(obj.h_valueCM, 'Label', 'Add', 'Callback', @obj.valueCM_callback);
      uimenu(obj.h_valueCM, 'Label', 'Edit', 'Callback', @obj.valueCM_callback);
      set(obj.h_valuePM,'uicontextmenu',obj.h_valueCM)

    end
    
    function delete(obj)
      try
        delete(obj.h_valuePM);
      end
      try
        delete(obj.h_valueCM);
      end
    end
    
    function valueCM_callback(obj,status,event)
      fprintf('valueCM_callback\n');
      if strcmpi(event.Source.Label,'Add')
        if obj.style_mode == 1
          obj.value = [];
          set(obj.h_valuePM,'Style','edit');
          set(obj.h_valuePM,'String',{''});
          set(obj.h_valuePM,'uicontextmenu',[]);
          obj.style_mode = 0;
        end
      elseif strcmpi(event.Source.Label,'Edit')
        if obj.style_mode == 1
          obj.value = get(obj.h_valuePM,'Value');
          set(obj.h_valuePM,'Style','edit');
          set(obj.h_valuePM,'String',obj.list_values(obj.value));
          set(obj.h_valuePM,'uicontextmenu',[]);
          obj.style_mode = 0;
        end
      end
    end
    
    function valuePM_callback(obj,status,event)
      %fprintf('valuePM_callback\n');
      if obj.style_mode == 0
        new_str = get(obj.h_valuePM,'String');
        if ~isempty(new_str{1})
          if isempty(obj.value)
          obj.list_values{end+1} = new_str{1};
          obj.value = length(obj.list_values);
          else
            obj.list_values{obj.value} = new_str{1};
          end
          obj.update_PM;
        end
      end
    end
    
    function key_press(obj,status,event)
      %fprintf('keyboard_callback %s\n', event.Key);
      if ~isempty(event.Key)
        switch event.Key
          case 'return'
            %fprintf('  return\n');
            if obj.style_mode == 1
              obj.value = [];
              set(obj.h_valuePM,'Style','edit');
              set(obj.h_valuePM,'String',{''});
              set(obj.h_valuePM,'uicontextmenu',[]);
              obj.style_mode = 0;
            end
          case 'backspace'
            %fprintf('  backspace\n');
            if obj.style_mode == 1
              obj.value = get(obj.h_valuePM,'Value');
              obj.list_values(obj.value) = [];
              obj.update_PM;
            end
        end
      end
    end
    
    function update_PM(obj)
      if isempty(obj.list_values)
        obj.value = [];
        set(obj.h_valuePM,'Style','edit');
        set(obj.h_valuePM,'String',{''});
        set(obj.h_valuePM,'uicontextmenu',[]);
        obj.style_mode = 0;
      else
        set(obj.h_valuePM,'String',obj.list_values);
        if isempty(obj.value)
          obj.value = 1;
        end
        obj.value = min(obj.value,length(obj.list_values));
        set(obj.h_valuePM,'Value',obj.value);
        set(obj.h_valuePM,'Style','popupmenu');
        set(obj.h_valuePM,'uicontextmenu',obj.h_valueCM);
        obj.style_mode = 1;
      end
    end
    
    function val = get(obj,prop)
      if strcmpi(prop,'String')
        val = get_list(obj);
      elseif strcmpi(prop,'Value')
        val = get_value(obj);
      else
        keyboard
      end
    end
    
    function val = set(obj,prop,value)
      if strcmpi(prop,'String')
        set_list(obj,value);
      elseif strcmpi(prop,'Enable')
        set_enable(obj,value)
      elseif strcmpi(prop,'Value')
        set_value(obj,value);
      else
        keyboard
      end
    end
    
    function val = get_value(obj)
      val = obj.value;
    end
    
    function set_value(obj, val)
      obj.value = val;
      obj.update_PM;
    end
    
    function val = get_selected_string(obj)
      if isempty(obj.list_values)
        val = [];
      else
        val = obj.list_values{obj.value};
      end
    end
    
    function set_selected_string(obj, selected)
      obj.value = find(strcmp(selected,obj.list_values),1);
      obj.update_PM;
    end
    
    function vals = get_list(obj, vals)
      vals = obj.list_values;
    end
    
    function set_list(obj, vals)
      obj.list_values = vals;
      obj.update_PM;
    end
    
    function set_enable(obj,enable_mode)
      if ischar(enable_mode)
        if strcmpi(enable_mode,'on')
          enable_mode = true;
        else
          enable_mode = false;
        end
      end
      if obj.enable_mode ~= enable_mode
        obj.enable_mode = enable_mode;
        if obj.enable_mode
          set(obj.h_valuePM,'Enable','on');
          set(obj.h_valuePM,'uicontextmenu',obj.h_valueCM);
        else
          set(obj.h_valuePM,'Enable','off');
          set(obj.h_valuePM,'uicontextmenu',[]);
        end
      end
    end
    
  end
  
end
