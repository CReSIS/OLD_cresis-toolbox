classdef (HandleCompatible = true) slicetool_threshold < imb.slicetool
   
    properties
        sb
        slice
        theta
        img
        Time
    end
    
    properties (SetAccess = private, GetAccess = private)
    end
    
    events
    end
    
    methods
        function obj = slicetool_threshold()
            obj.create_option_ui();
            obj.tool_name = 'Threshold';
            obj.tool_menu_name = '(T)hreshold';
            obj.tool_shortcut = 't';
        end
        
        function cmd = apply_PB_callback(obj,sb)
            
            fprintf('Apply %s to layer %d slice %d\n', obj.tool_name, sb.layer_idx, sb.slice);
            control_idx = sb.layer(sb.layer_idx).control_layer;
            surf_idx = sb.layer(sb.layer_idx).surf_layer;
            rline = sb.slice;
            twtt_sur_all = tomo.threshold(obj.theta,obj.img, obj.Time,sb.slice);
            
            % Create cmd for layer change
            cmd = [];
            cmd{1}.undo.slice = sb.slice;
            cmd{1}.redo.slice = sb.slice;
            cmd{1}.undo.layer = sb.layer_idx;
            cmd{1}.redo.layer = sb.layer_idx;
            cmd{1}.undo.x = 1:size(sb.layer(sb.layer_idx).y,1);
            cmd{1}.undo.y = sb.layer(sb.layer_idx).y(:,sb.slice);
            cmd{1}.redo.x = 1:size(sb.layer(sb.layer_idx).y,1);
            cmd{1}.redo.y =  interp1(obj.Time,1:length(obj.Time),twtt_sur_all(:,sb.slice));
            cmd{1}.type = 'standard';           
            fprintf('Applied Threshold\n');
            
        end
        
        function set_custom_data(obj,custom_data)
            obj.custom_data.ice_mask = custom_data.ice_mask;
            obj.theta = custom_data.theta ;
            obj.img = custom_data.img ;
            obj.Time = custom_data.Time;
        end
     
        function create_option_ui(obj)
        end
        
    end
    
end