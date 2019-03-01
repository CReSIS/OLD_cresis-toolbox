% Class coverage_maps_master_fig_files_loader
%
% Class which creates a GUI for editing the master fig file
%
% obj = coverage_maps_master_fig_files_loader_new_script('C:\tmp\coverage_maps\Greenland_All_Seasons.fig');
% obj = coverage_maps_master_fig_files_loader(1);

classdef (HandleCompatible = true) coverage_maps_master_fig_files_loader_new_script < handle
  properties
    h_fig % Figure handle
    h_ctrl % Control figure handle
    h_children % Children of h_fig
    h_geotiff % Geotiff handle
    h_gui % h_ctrl GUI handles struct
    
    season_names % cell vector of strings containing season names
    delete_figure_on_obj_delete % logical (delete h_fig on obj deletion)
  end
  
  methods
    function obj = coverage_maps_master_fig_files_loader_new_script(master_fig)

      obj.delete_figure_on_obj_delete = false;
      if ischar(master_fig)
        % User passed in a filename
        obj.h_fig = open(master_fig);
        obj.delete_figure_on_obj_delete = true;
      elseif ishandle(master_fig);
        % User passed in a handle
        obj.h_fig = master_fig;
      else
        error('Argument must be a figure handle or filename');
      end
      
      %% Parse master figure UserData
      obj.season_names = flipud(get(obj.h_fig,'UserData'));
      for season_idx = 1:length(obj.season_names)
        obj.season_names{season_idx} = sprintf('%d: %s', season_idx, obj.season_names{season_idx});
      end
      obj.h_children = get(obj.h_fig,'Children');
      obj.h_children = get(obj.h_children(2),'Children');
      obj.h_geotiff = obj.h_children(end);
      % Sort children handles into 3xN matrix where N is the length of
      % season_names and row 1 is good data, row 2 is moderate quality data,
      % and row 3 is bad data (no bottom)
      obj.h_children = reshape(obj.h_children(1:end-5),[4 length(obj.season_names)]);

      %% Create control figure
      obj.h_ctrl = figure;
      h_ctrl_pos = get(obj.h_ctrl,'Position');
      set(obj.h_ctrl,'Position',[h_ctrl_pos(1:2) 200 400]);
      
      obj.h_gui.seasonLB = uicontrol('Parent',obj.h_ctrl);
      set(obj.h_gui.seasonLB,'Style','listbox');
      set(obj.h_gui.seasonLB,'HorizontalAlignment','Center');
      set(obj.h_gui.seasonLB,'FontName','fixed');
      set(obj.h_gui.seasonLB,'String',obj.season_names);
      set(obj.h_gui.seasonLB,'Value',1:length(obj.season_names));
      set(obj.h_gui.seasonLB,'Callback',@obj.seasonsLB_callback);
      set(obj.h_gui.seasonLB,'Max',1e9);
      
      obj.h_gui.viewPM = uicontrol('Parent',obj.h_ctrl);
      set(obj.h_gui.viewPM,'Style','PopupMenu');
      menu_string{1} = 'Red';
      menu_string{2} = 'Magenta';
      menu_string{3} = 'Yellow';
      menu_string{4} = 'Green';
      menu_string{5} = 'All';
      menu_string{6} = 'Clear';
      menu_string{7} = 'Clear All';
      set(obj.h_gui.viewPM,'String',menu_string);
      set(obj.h_gui.viewPM,'Value',1);
      set(obj.h_gui.viewPM,'HorizontalAlignment','Center');
      set(obj.h_gui.viewPM,'Value',1);
      set(obj.h_gui.viewPM,'Callback',@obj.viewPM_callback);
      
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
      obj.h_gui.h_table.handles{row,col}   = obj.h_gui.viewPM;
      obj.h_gui.h_table.width(row,col)     = inf;
      obj.h_gui.h_table.height(row,col)    = 25;
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
        delete(obj.h_ctrl);
      catch ME
        disp(ME)
      end
      if obj.delete_figure_on_obj_delete
        try
          delete(obj.h_fig);
        catch ME
          disp(ME)
        end
      end
    end
    
    function close_win(obj,h_obj,event)
      delete(obj);
    end
    
    function seasonsLB_callback(obj,h_obj,event)
      refresh_view(obj);
    end
    
    function viewPM_callback(obj,h_obj,event)
      refresh_view(obj);
    end
    
    function refresh_view(obj)
      seasons_selected = get(obj.h_gui.seasonLB,'Value');
      view_selected = get(obj.h_gui.viewPM,'Value');
      selected_mask = logical(zeros(size(obj.season_names)));
      selected_mask(seasons_selected) = 1;
      %       menu_string{1} = 'Green';
      %       menu_string{2} = 'Yellow';
      %       menu_string{3} = 'Red';
      %       menu_string{4} = 'Magenta';
      %       menu_string{5} = 'All';
      % RED
      %       if view_selected == 1
      %         temp = (obj.h_children(1,selected_mask));
      %         flag = 1;
      %         for i=1:length(temp)
      %           if((strcmp(temp(i).('Visible'),'off'))==1)
      %             set(obj.h_children(1,selected_mask),'Visible','on');
      %             flag = 0;
      %             break;
      %           end
      %         end
      %         if(flag)
      %           set(obj.h_children(1,selected_mask),'Visible','off');
      %         end
      if( view_selected <= 4)
        
        temp = (obj.h_children(view_selected,selected_mask));
        set(obj.h_children(view_selected,selected_mask),'Visible','on');
        set(obj.h_children(view_selected,~selected_mask),'Visible','off');
        %
        %         temp2 = (obj.h_children(view_selected,~selected_mask));
         
        %Find out how many seasons are already switched on for that color
%         all_masks = selected_mask + ~selected_mask;
%         temp_view_selected = (obj.h_children(view_selected,all_masks));
%         counter = 0;
%         for i=1:length(temp_view_selected)
%           if((strcmp(temp_view_selected(i).('Visible'),'on'))==1)
%             counter = counter+1;
%           end
%         end
%         
%         flag = 1;
%         for i=1:length(temp)
%           if((strcmp(temp(i).('Visible'),'off'))==1)
%             set(obj.h_children(view_selected,selected_mask),'Visible','on');
%             set(obj.h_children(view_selected,~selected_mask),'Visible','off');
% 
%             flag = 0;
%             break;
%           end
%         end
%         if (flag && counter == length(all_masks))
%             set(obj.h_children(view_selected,selected_mask),'Visible','on');
%             set(obj.h_children(view_selected,~selected_mask),'Visible','off');
%         elseif (flag && length(temp)>1)
%             set(obj.h_children(view_selected,selected_mask),'Visible','on');
%             set(obj.h_children(view_selected,~selected_mask),'Visible','off');
%         end
%         if(flag && length(temp)>1)
%           set(obj.h_children(view_selected,selected_mask),'Visible','off');
% %         elseif (flag && length(temp)>1 && length(temp)>sum(selected_mask(:)==1))
% %           set(obj.h_children(view_selected,~selected_mask),'Visible','off');
%         elseif (flag && length(temp)==1)
%           if(flag)
%             set(obj.h_children(view_selected,selected_mask),'Visible','off');
%           end
%         end

        
      else
        if view_selected == 5
          set(obj.h_children(1:4,selected_mask),'Visible','on');
                  set(obj.h_children(1:4,~selected_mask),'Visible','off');
        elseif view_selected == 6
                  set(obj.h_children(4,selected_mask),'Visible','on');
          set(obj.h_children(1:4, selected_mask),'Visible','off');
        elseif view_selected == 7
          %         set(obj.h_children(4,selected_mask),'Visible','on');
          selected_mask(1:length(selected_mask)) = 1;
          set(obj.h_children(1:4, selected_mask),'Visible','off');
        end
      end

     
    end
  end
  
end

