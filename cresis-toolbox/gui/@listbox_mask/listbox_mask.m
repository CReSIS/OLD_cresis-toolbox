classdef listbox_mask < handle
% listbox_mask: A class for creating an integrated label (optional)
% listbox, and line edit with a search field.  Whatever is typed
% in the search field causes the contents of the listbox to change
% to what matches this field. A secondary mask can optionally be
% applied which restricts the fields further.  Supports a separate
% sorting index as well.
%
% parent: parent GUI handle for each of the GUI objects
% name: label (no label if left blank)
% list_values: cell vector of strings
% init_value: scalar integer index into list_values vector
%
% figure(1); clf;
% myLBM = listbox_mask(1,'Test',{'1','2','3','2b'},[1 1 1 1],[1 2 4 3],1);
% %% Create the table
% obj.h_gui.table.ui=1;
% obj.h_gui.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
% obj.h_gui.table.height_margin = NaN*zeros(30,30);
% obj.h_gui.table.false_width = NaN*zeros(30,30);
% obj.h_gui.table.false_height = NaN*zeros(30,30);
% obj.h_gui.table.offset = [0 0];
% 
% row = 1; col = 1;
% obj.h_gui.table.handles{row,col}   = myLBM.h_text;
% obj.h_gui.table.width(row,col)     = inf;
% obj.h_gui.table.height(row,col)    = 20;
% obj.h_gui.table.width_margin(row,col) = 1;
% obj.h_gui.table.height_margin(row,col) = 1;
% 
% row = row + 1; col = 1;
% obj.h_gui.table.handles{row,col}   = myLBM.h_list;
% obj.h_gui.table.width(row,col)     = inf;
% obj.h_gui.table.height(row,col)    = inf;
% obj.h_gui.table.width_margin(row,col) = 1;
% obj.h_gui.table.height_margin(row,col) = 1;
% 
% row = row + 1; col = 1;
% obj.h_gui.table.handles{row,col}   = myLBM.h_LE;
% obj.h_gui.table.width(row,col)     = inf;
% obj.h_gui.table.height(row,col)    = 20;
% obj.h_gui.table.width_margin(row,col) = 1;
% obj.h_gui.table.height_margin(row,col) = 1;
% 
% table_draw(obj.h_gui.table);

  properties
    h_text
    h_list
    h_LE
    list_values
    cur_value
    cur_mask
    search_mask
    sort_idxs
    double_click
  end
  
  properties (SetAccess = private, GetAccess = private)
  end
   
  events
    list_changed
  end
  
  methods
    function obj = listbox_mask(parent,name,list_values,init_mask,sort_idxs,init_value)
      %%% Pre Initialization %%%
      % Any code not using output argument (obj)
      
      %%% Post Initialization %%%
      % Any code, including access to object
      
      double_click = false;
      
      if isempty(name)
        obj.h_text = 0;
      else
        obj.h_text = uicontrol('Parent',parent);
        set(obj.h_text,'Style','Text');
        set(obj.h_text,'String',name);
      end
    
      obj.list_values = list_values;
      if ~exist('init_value','var') || isempty(init_value)
        obj.cur_value = 1;
      else
        obj.cur_value = init_value;
      end
      if ~exist('init_mask','var') || isempty(init_mask)
        obj.cur_mask = logical(ones(size(obj.list_values)));
      else
        obj.cur_mask = logical(init_mask);
      end
      obj.cur_mask = reshape(obj.cur_mask,[1 length(obj.cur_mask)]);
      if isempty(obj.cur_mask)
        obj.cur_mask = [];
      end
      if ~exist('sort_idxs','var') || isempty(sort_idxs)
        obj.sort_idxs = 1:length(obj.list_values);
      else
        obj.sort_idxs = sort_idxs;
      end
      
      obj.search_mask = ones(size(obj.list_values));
      
      obj.h_list = uicontrol('Parent',parent);
      set(obj.h_list,'Style','listbox');
      set(obj.h_list,'Callback',@obj.list_callback);
      
      obj.h_LE = uicontrol('Parent',parent);
      set(obj.h_LE,'Style','Edit');
      set(obj.h_LE,'String','.*');
      set(obj.h_LE,'Callback',@obj.list_callback);
      set(obj.h_LE,'TooltipString',sprintf('See "help regexpi".\n  ".*" selects all\n  "0519" selects all with the string 0519 occuring anywhere\n  "^2012" selects all starting with 2012\n  "(?!^2016.*)^.*" selects all besides 2016'));

      obj.update_listbox();
    end
    
    function delete(obj)
    end

    function list_callback(obj, h_obj, event)
      list_string = get(obj.h_list,'String');
      cur_entry = get(obj.h_list,'Value');
      obj.double_click = false;
      if h_obj == obj.h_list
        % User selected a value in the list
        if strcmpi(get(gcf,'SelectionType'),'Open')
          obj.double_click = true;
        else
          set(obj.h_list,'Max',1);
          good_mask = obj.search_mask & obj.cur_mask;
          sorted_mask = good_mask(obj.sort_idxs);
          cur_value = find(sorted_mask,cur_entry);
          obj.cur_value = obj.sort_idxs(cur_value(end));
        end
      elseif h_obj == obj.h_LE
        % User entered a new search criteria
        obj.search_mask = cellfun(@(x) ~isempty(x),regexpi(obj.list_values,get(obj.h_LE,'String')));
        obj.update_listbox();
      end
      notify(obj,'list_changed');
    end

    function update_listbox(obj)
      good_mask = obj.search_mask & obj.cur_mask;
      sorted_mask = good_mask(obj.sort_idxs);
      
      list_values_sorted = obj.list_values(obj.sort_idxs);
      
      set(obj.h_list,'Value',1);
      set(obj.h_list,'String',list_values_sorted(sorted_mask));
      cur_val_sorted = find(obj.sort_idxs == obj.cur_value);
      if sorted_mask(cur_val_sorted)
        cur_val_listbox = sum(sorted_mask(1:cur_val_sorted));
        set(obj.h_list,'Value',cur_val_listbox);
        set(obj.h_list,'Max',1);
      else
        set(obj.h_list,'Max',2);
        set(obj.h_list,'Value',[]);
      end
    end
    
    function set_value(obj, new_val)
      if isempty(new_val)
        obj.cur_value = 1;
      else
        obj.cur_value = new_val;
      end
      obj.update_listbox();
    end
    
    function set_mask(obj,new_mask)
      if isempty(new_mask)
        obj.cur_mask = logical(ones(size(obj.list_values)));
      else
        obj.cur_mask = logical(new_mask);
      end
      obj.cur_mask = reshape(obj.cur_mask,[1 length(obj.cur_mask)]);
      if isempty(obj.cur_mask)
        obj.cur_mask = [];
      end
      obj.update_listbox();
    end
    
    function combined_mask = get_combined_mask(obj)
      combined_mask = obj.cur_mask & obj.search_mask;
    end
    
    function set_sort(obj,new_sort_idxs)
      if isempty(new_sort_idxs)
        obj.sort_idxs = 1:length(obj.list_values);
      else
        obj.sort_idxs = new_sort_idxs;
      end
      obj.update_listbox();
    end
       
    function set_list_values(obj,list_values,init_mask,sort_idxs,init_value)
      obj.list_values = list_values;
      if ~exist('init_value','var') || isempty(init_value)
        obj.cur_value = 1;
      else
        obj.cur_value = init_value;
      end
      if ~exist('init_mask','var') || isempty(init_mask)
        obj.cur_mask = logical(ones(size(obj.list_values)));
      else
        obj.cur_mask = logical(init_mask);
      end
      obj.cur_mask = reshape(obj.cur_mask,[1 length(obj.cur_mask)]);
      if isempty(obj.cur_mask)
        obj.cur_mask = [];
      end
      if ~exist('sort_idxs','var') || isempty(sort_idxs)
        obj.sort_idxs = 1:length(obj.list_values);
      else
        obj.sort_idxs = sort_idxs;
      end
      
      obj.search_mask = cellfun(@(x) ~isempty(x),regexpi(obj.list_values,get(obj.h_LE,'String')));
      
      obj.update_listbox();
    end
    
    function reset(obj, enable_notify)
      % User entered a new search criteria
      set(obj.h_LE,'String','.*');
      obj.search_mask = cellfun(@(x) ~isempty(x),regexpi(obj.list_values,get(obj.h_LE,'String')));
      if exist('enable_notify','var') || isempty(enable_notify) || enable_notify
        obj.update_listbox();
        notify(obj,'list_changed');
      end
    end
    
  end
  
end


