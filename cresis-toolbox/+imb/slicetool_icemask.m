classdef (HandleCompatible = true) slicetool_icemask < imb.slicetool
  % Slice_browser tool which calls detect.cpp (HMM)
  
  properties
      ice
      cmd
      ice_layer
      sb
  end
  
  properties (SetAccess = private, GetAccess = private)
  end
  
  events
    IceChange
  end
  
  methods
    function obj = slicetool_icemask()
      obj.create_option_ui();
      obj.tool_name = 'Ice Mask';
      obj.tool_menu_name = '(I)ce Mask';
      obj.tool_shortcut = 'i';
    end
    
    function cmd = apply_PB_callback(obj,sb)
      obj.sb = sb;
      layer_name = {sb.layer.name};
      
      % Read and set sb.layer
      % Read sb.data
      % Read sb.slice
      % Read sb.layer_idx
      fprintf('Apply %s to layer %d slice %d\n', obj.tool_name, sb.layer_idx, sb.slice);
      
%       
    end
    
    function set_custom_data(obj,custom_data)
      param.mdata.twtt = custom_data.mdata.twtt;
      param.mdata.theta = custom_data.mdata.theta;
      param.mdata.ice_mask = custom_data.mdata.ice_mask;
      param.mdata.param_combine = custom_data.mdata.param_combine;
%       rmfield(param,'mdata');
      param.DEM = custom_data.DEM;
      param.R = custom_data.R;
      param.ice_mask_fn = custom_data.ice_mask_fn;
      param.ice_mask = custom_data.ice_mask;
      param.proj = custom_data.proj;
      obj.sb = custom_data.sb;
      addlistener(obj.sb,'SliceChange',@obj.runChangeSlice);
      addlistener(obj.sb.undo_stack,'synchronize_event',@obj.run_undo_sync);
      obj.sb.save_callback = @obj.save;
      
      obj.ice = imb.ice_mask_edit(param);
      addlistener(obj.ice,'IceChange',@obj.run_ice_change);
      addlistener(obj.ice,'SliceChange',@obj.run_slice_change);
      addlistener(obj.ice,'Undo',@obj.run_sb_undo);
      addlistener(obj.ice,'Redo',@obj.run_sb_redo);
      obj.ice.local_undo_flag = 0;
    end
    
    function create_option_ui(obj)
    end
    
    function add_listener(obj,src)
      addlistener(src,'SliceChange',@obj.runChangeSlice);
      addlistener(src.undo_stack,'synchronize_event',@obj.run_undo_sync);
    end
    
    function evnts = get_events(obj)
      evnts.src = obj;
      evnts.evnts{1} = 'IceChange';
    end
    
    function runChangeSlice(obj,src,~)
      obj.ice.change_slice(src.slice);
    end
    
    function cmd = overlay_intersects(obj,cmd)
      
      thetas = [];
      slices = [];
      vals = [];
      val_olds = [];
      cmd_op_idxs = [];
      % get theta and slice idxs for changes
      for cmd_idx = 1:length(cmd)
        if strcmp(cmd{cmd_idx}.type,'standard') && ...
            cmd{cmd_idx}.redo.layer == obj.sb.layer(cmd{cmd_idx}.redo.layer).mask_layer
          thetas = [thetas , cmd{cmd_idx}.redo.x];
          slices = [slices , cmd{cmd_idx}.redo.slice];
          vals = [vals , cmd{cmd_idx}.redo.y];
          val_olds = [val_olds , cmd{cmd_idx}.undo.y];
          cmd_op_idxs = [cmd_op_idxs , cmd_idx];
        end
      end

      if isempty(thetas)
        return
      end

      cmd_ice = obj.ice.edit_twtt(thetas,slices,vals,val_olds);
      cmd_append = obj.form_sb_cmd(cmd_ice,cmd{1}.redo.layer);
      
      % check for same slices
      del_cmd_idx = [];
      for i = 1:length(cmd_append)
        if strcmp(cmd_append{i}.type,'standard') && ...
            ismember(cmd_append{i}.redo.slice,slices)
          cmd{cmd_append{i}.redo.slice==slices} = cmd_append{i};
          del_cmd_idx = [del_cmd_idx,i];
        end
      end
      cmd_append(del_cmd_idx) = [];
      
      cmd = [cmd,cmd_append];
    end 
    
    
    function cmd = push_request(obj,cmd)
      cmd_ice = {};
      for cmd_idx = 1:length(cmd)
        if strcmp(cmd{cmd_idx}.type,'ice_mask')
          return;
        end
      end
      
      if isempty(cmd_ice)
        for cmd_idx = 1:length(cmd)
          if strcmp(cmd{cmd_idx}.type,'standard') && ...
              cmd{cmd_idx}.redo.layer == obj.sb.layer(cmd{cmd_idx}.redo.layer).mask_layer            
            
              cmd = obj.overlay_intersects(cmd);
              return
          end
        end
      end
      
    end
    
    
    function run_undo_sync(obj,src,~)
      obj.ice.undo_sync(src,[]);
    end
    
    
    function cmd = form_sb_cmd(obj,cmd,layer)
      
      for i = 1:length(cmd)
        if strcmp(cmd{i}.type,'ice_mask')
          ice_cmd_idx = i;
        end
      end
      
      if isempty(cmd{1}.redo.data_mask_idx)
        return;
      end
      
      [theta_idx,slice_idx] = ind2sub([size(obj.ice.intersections,2),size(obj.ice.intersections,3)],cmd{1}.redo.data_mask_idx);
      unique_slices = unique(slice_idx);
      
      for s = 1:length(unique_slices)
        slice = unique_slices(s);
        
        cmd_idx = length(cmd)+1;
        
        s_idx = find(slice_idx==slice);
        theta = theta_idx(s_idx);
%           ind = sub2ind([size(obj.intersections,2),size(obj.intersections,3)],theta,s_idx);

        cmd{cmd_idx}.undo.layer = layer;
        cmd{cmd_idx}.redo.layer = layer;
        cmd{cmd_idx}.undo.x = theta;
        cmd{cmd_idx}.redo.x = theta;
        cmd{cmd_idx}.undo.y = logical(obj.ice.mdata.ice_mask(theta,slice));
        cmd{cmd_idx}.redo.y = logical(cmd{1}.redo.data_mask(s_idx));
        cmd{cmd_idx}.undo.slice = slice;
        cmd{cmd_idx}.redo.slice = slice;
        cmd{cmd_idx}.type = 'standard';
      end
    end
    
    
    function run_ice_change(obj,cmd,~)
      cmd = obj.form_sb_cmd(obj.ice.cmd,obj.sb.layer(1).mask_layer);
      obj.sb.push(cmd);
    end
    
    
    function run_slice_change(obj,src,~)
      obj.sb.slice = src.slice;
      obj.sb.update_slice();
    end
    
    
    function run_sb_redo(obj,src,~)
      obj.sb.undo_stack.redo();
    end
    
    
    function run_sb_undo(obj,src,~)
      obj.sb.undo_stack.pop();
    end
    
    function save(obj)
      obj.ice.save();
    end
    
  end
  
end