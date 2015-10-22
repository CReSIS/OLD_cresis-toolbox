function synch_layers(obj)
% Synchronize the pointer from echowin class with the pointer from undo_stack
% calss


for echowin_idx = 1:length(obj.echowin_list)
  if obj.echowin_list(echowin_idx).isvalid && ...
      obj.echowin_list(echowin_idx).undo_stack_ptr ~= obj.pointer
    cur_echowin = obj.echowin_list(echowin_idx);
    echo_pointer = cur_echowin.undo_stack_ptr;
    cur_frame_idxs = cur_echowin.eg.frame_idxs;
    while echo_pointer > obj.pointer
      %% =======================================================================
      % Undo
      %========================================================================
      % Pop the last change off of the undo stack
      % And make sure the undo is within current frames
      undo = obj.stack{echo_pointer};
      if (undo.frame_idxs(1) >= cur_frame_idxs(1) && undo.frame_idxs(1) <= cur_frame_idxs(end)) ||...
          (undo.frame_idxs(end) >= cur_frame_idxs(1) && undo.frame_idxs(end) <= cur_frame_idxs(end))
        % pull x_old_gps from undo, and make sure only pnts in current
        % frame is remained
        x_old_gps = undo.x_old_gps;
        y_old_twtt = undo.y_old_twtt;
        quality_old = undo.quality_old;
        manu_auto_old = undo.manu_auto_old;
        x_old_idx = (x_old_gps>=cur_echowin.eg.image_gps_time(1)) & ...
          (x_old_gps<=cur_echowin.eg.image_gps_time(end));
        x_old_gps = x_old_gps(x_old_idx);
        y_old_twtt = y_old_twtt(x_old_idx);
        quality_old = quality_old(x_old_idx);
        manu_auto_old = manu_auto_old(x_old_idx);
        % pull new pnts from undo, and make sure only pnts in current
        % frame is remained
        x_new_gps = undo.x_new_gps;
        x_new_idx = (x_new_gps>=cur_echowin.eg.image_gps_time(1)) & ...
          (x_new_gps<=cur_echowin.eg.image_gps_time(end));
        x_new_gps = x_new_gps(x_new_idx);
        % pull current layer
        cur_layer = undo.cur_layer;
        %========================================================================
        % Uptdate layer data, layer quality, and layer M_A
        layer_data_y = cur_echowin.eg.layer.y{cur_layer};
        layer_data_x = cur_echowin.eg.layer.x{cur_layer};
        layerQual = cur_echowin.eg.layer.qual{cur_layer};
        layerM_A = cur_echowin.eg.layer.type{cur_layer};
        % remove new pnts
        if ~isempty(x_new_gps)
          for idx = 1:length(x_new_gps)
            x_idx = logical(x_new_gps(idx) - layer_data_x);
            % but don't remove manual pts if operation is auto!
%             if any(manu_auto_old==2) % right now tools set all pts to auto in an auto operation so this should be ok (i.e. no mixed cmds)
%               if any(layerM_A(~x_idx)==1)
%                 keep_mask = logical(layerM_A(~x_idx)-2);
%                 x_idx(~x_idx) = keep_mask;
%               end
%             end
            layer_data_x = layer_data_x(x_idx);
            layer_data_y = layer_data_y(x_idx);
            layerQual = layerQual(x_idx);
            layerM_A = layerM_A(x_idx);
          end
        end
        % Add old pnts
        if ~isempty(x_old_gps)
          [layer_data_x,data_idx] = sort([layer_data_x,x_old_gps]);
          tmp = [layer_data_y,y_old_twtt];
          layer_data_y = tmp(data_idx);
          tmp = [layerQual,quality_old];
          layerQual = tmp(data_idx);
          tmp = [layerM_A,manu_auto_old];
          layerM_A = tmp(data_idx);
        end
        %=========================================================================
        % Convert layer_data to current unit
        % convert layer_data_x
        layer_x_curUnit = interp1(cur_echowin.eg.image_gps_time,...
          cur_echowin.eg.image_xaxis,...
          layer_data_x,'linear');
        layer_x_curUnit_manual = layer_x_curUnit(layerM_A == 1);
        layer_x_curUnit_auto = layer_x_curUnit(layerM_A == 2);
        % convert layer_data_y
        yaxis_choice = get(cur_echowin.left_panel.yaxisPM,'Value');
        layer_y_curUnit = zeros(1,length(layer_data_y));
        if yaxis_choice == 1 % TWTT
          layer_y_curUnit = layer_data_y*1e6;
          layer_y_curUnit_manual = layer_y_curUnit(layerM_A == 1);
          layer_y_curUnit_auto = layer_y_curUnit(layerM_A == 2);
        elseif yaxis_choice == 2 % WGS_84 Elevation
          elevation = interp1(cur_echowin.eg.gps_time,...
            cur_echowin.eg.elevation,...
            layer_data_x,'linear');
          surface = interp1(cur_echowin.eg.gps_time,...
            cur_echowin.eg.surface,...
            layer_data_x,'linear');
          physical_constants;
          for idx = 1: length(layer_data_y)
            if layer_data_y(idx) <= surface(idx)
              layer_y_curUnit(idx) = elevation(idx) - layer_data_y(idx)*c/2;
            else
              layer_y_curUnit(idx) = elevation(idx) - surface(idx)*c/2 - ...
                (layer_data_y(idx) - surface(idx))*c/(sqrt(er_ice)*2);
            end
          end
          layer_y_curUnit_manual = layer_y_curUnit(layerM_A == 1);
          layer_y_curUnit_auto = layer_y_curUnit(layerM_A == 2);
        elseif yaxis_choice == 3 % Depth/Range
          surface = interp1(cur_echowin.eg.gps_time,...
            cur_echowin.eg.surface,...
            layer_data_x,'linear');
          physical_constants;
          for idx = 1: length(layer_data_y)
            if layer_data_y(idx) <= surface(idx)
              layer_y_curUnit(idx) = layer_data_y(idx)*c/2;
            else
              layer_y_curUnit(idx) = surface(idx)*c/2 + (layer_data_y(idx) - surface(idx))*c/(sqrt(er_ice)*2);
            end
          end
          layer_y_curUnit_manual = layer_y_curUnit(layerM_A == 1);
          layer_y_curUnit_auto = layer_y_curUnit(layerM_A == 2);
        elseif yaxis_choice == 4 % Range bin
          layer_y_curUnit = interp1(cur_echowin.eg.time,cur_echowin.eg.image_yaxis,...
            layer_data_y,'linear');
          layer_y_curUnit_manual = layer_y_curUnit(layerM_A == 1);
          layer_y_curUnit_auto = layer_y_curUnit(layerM_A == 2);
        end
        %=======================================================================
        % Update layer handle
        % manual pts handle
        set(cur_echowin.layer_h(2*cur_layer-1),{'XData','YData'},{layer_x_curUnit_manual,layer_y_curUnit_manual});
        % auto pts handle
        set(cur_echowin.layer_h(2*cur_layer),{'XData','YData'},{layer_x_curUnit_auto,layer_y_curUnit_auto});
        % Update quality handles
        good_idx_manual = layerM_A == 1 & ~logical(layerQual-1);
        good_idx_auto = layerM_A == 2 & ~logical(layerQual-1);
        good_x_manual = layer_x_curUnit(good_idx_manual);
        good_x_auto = layer_x_curUnit(good_idx_auto);
        good_y_manual = layer_y_curUnit(good_idx_manual);
        good_y_auto = layer_y_curUnit(good_idx_auto);
        moderate_idx_manual = layerM_A == 1 & ~logical(layerQual-2);
        moderate_idx_auto = layerM_A == 2 & ~logical(layerQual-2);
        moderate_x_manual = layer_x_curUnit(moderate_idx_manual);
        moderate_x_auto = layer_x_curUnit(moderate_idx_auto);
        moderate_y_manual = layer_y_curUnit(moderate_idx_manual);
        moderate_y_auto = layer_y_curUnit(moderate_idx_auto);
        bad_idx_manual = layerM_A == 1 & ~logical(layerQual-3);
        bad_idx_auto = layerM_A == 2 & ~logical(layerQual-3);
        bad_x_manual = layer_x_curUnit(bad_idx_manual);
        bad_x_auto = layer_x_curUnit(bad_idx_auto);
        bad_y_manual = layer_y_curUnit(bad_idx_manual);
        bad_y_auto = layer_y_curUnit(bad_idx_auto);
        set(cur_echowin.quality_h(6*(cur_layer-1)+1),{'XData','YData'},{good_x_manual,good_y_manual});
        set(cur_echowin.quality_h(6*(cur_layer-1)+2),{'XData','YData'},{good_x_auto,good_y_auto});
        set(cur_echowin.quality_h(6*(cur_layer-1)+3),{'XData','YData'},{moderate_x_manual,moderate_y_manual});
        set(cur_echowin.quality_h(6*(cur_layer-1)+4),{'XData','YData'},{moderate_x_auto,moderate_y_auto});
        set(cur_echowin.quality_h(6*(cur_layer-1)+5),{'XData','YData'},{bad_x_manual,bad_y_manual});
        set(cur_echowin.quality_h(6*(cur_layer-1)+6),{'XData','YData'},{bad_x_auto,bad_y_auto});
        %=========================================================================
        % Update gCtrl layer data vector
        cur_echowin.eg.layer.y{cur_layer} = layer_data_y;
        cur_echowin.eg.layer.x{cur_layer} = layer_data_x;
        cur_echowin.eg.layer.qual{cur_layer} = layerQual;
        cur_echowin.eg.layer.type{cur_layer} = layerM_A;
        cur_echowin.eg.layer.x_curUnit{cur_layer} = layer_x_curUnit;
        cur_echowin.eg.layer.y_curUnit{cur_layer} = layer_y_curUnit;
      end
      %=======================================================================
      % Update stack pointer
      echo_pointer = echo_pointer-1;
    end
    while echo_pointer < obj.pointer
      %% ========================================================================
      % Redo
      %=======================================================================
      % Update stack pointer
      echo_pointer = echo_pointer + 1;
      %========================================================================
      % Pop the last change off of the undo stack
      % And make sure the undo is within current frames
      undo = obj.stack{echo_pointer};
      if (undo.frame_idxs(1) >= cur_frame_idxs(1) && undo.frame_idxs(1) <= cur_frame_idxs(end)) ||...
          (undo.frame_idxs(end) >= cur_frame_idxs(1) && undo.frame_idxs(end) <= cur_frame_idxs(end))
        % pull x_old_gps from undo, and make sure only pnts in current
        % frame remain
        x_old_gps = undo.x_old_gps;
        x_old_idx = (x_old_gps>=cur_echowin.eg.image_gps_time(1)) & ...
          (x_old_gps<=cur_echowin.eg.image_gps_time(end));
        x_old_gps = x_old_gps(x_old_idx);
        % pull new pnts from undo, and make sure only pnts in current
        % frame is remained
        x_new_gps = undo.x_new_gps;
        x_new_idx = (x_new_gps>=cur_echowin.eg.image_gps_time(1)) & ...
          (x_new_gps<=cur_echowin.eg.image_gps_time(end));
        x_new_gps = x_new_gps(x_new_idx);
        y_new_twtt = undo.y_new_twtt;
        y_new_twtt = y_new_twtt(x_new_idx);
        quality_new = undo.quality_new;
        quality_new = quality_new(x_new_idx);
        manu_auto_new = undo.manu_auto_new;
        manu_auto_new = manu_auto_new(x_new_idx);
        % pull current layer
        cur_layer = undo.cur_layer;
        %========================================================================
        % Uptdate layer data, layer quality, and layer manual/auto
        layer_data_y = cur_echowin.eg.layer.y{cur_layer};
        layer_data_x = cur_echowin.eg.layer.x{cur_layer};
        layerQual = cur_echowin.eg.layer.qual{cur_layer};
        layerM_A = cur_echowin.eg.layer.type{cur_layer};
        % remove old pnts
        if ~isempty(x_old_gps)
          for idx = 1:length(x_old_gps)
            x_idx = logical(x_old_gps(idx) - layer_data_x);
            % but don't remove manual pts if operation is auto!
%             if any(manu_auto_new==2) % right now tools set all pts to auto in an auto operation so this should be ok (i.e. no mixed cmds)
%               if any(layerM_A(~x_idx)==1)
%                 keep_mask = logical(layerM_A(~x_idx)-2);
%                 x_idx(~x_idx) = keep_mask;
%               end
%             end
            layer_data_x = layer_data_x(x_idx);
            layer_data_y = layer_data_y(x_idx);
            layerQual = layerQual(x_idx);
            layerM_A = layerM_A(x_idx);
          end
        end
        % Add new pnts
        if ~isempty(x_new_gps)
          [layer_data_x,data_idx] = sort([layer_data_x,x_new_gps]);
          tmp = [layer_data_y,y_new_twtt];
          layer_data_y = tmp(data_idx);
          tmp = [layerQual,quality_new];
          layerQual = tmp(data_idx);
          tmp = [layerM_A,manu_auto_new];
          layerM_A = tmp(data_idx);
        end
        %=========================================================================
        % Convert layer_data to current unit
        % convert layer_data_x
        layer_x_curUnit = interp1(cur_echowin.eg.image_gps_time,...
          cur_echowin.eg.image_xaxis,...
          layer_data_x,'linear');
        layer_x_curUnit_manual = layer_x_curUnit(layerM_A == 1);
        layer_x_curUnit_auto = layer_x_curUnit(layerM_A == 2);
        % convert layer_data_y
        yaxis_choice = get(cur_echowin.left_panel.yaxisPM,'Value');
        layer_y_curUnit = zeros(1,length(layer_data_y));
        if yaxis_choice == 1 % TWTT
          layer_y_curUnit = layer_data_y*1e6;
          layer_y_curUnit_manual = layer_y_curUnit(layerM_A == 1);
          layer_y_curUnit_auto = layer_y_curUnit(layerM_A == 2);
        elseif yaxis_choice == 2 % WGS_84 Elevation
          elevation = interp1(cur_echowin.eg.gps_time,...
            cur_echowin.eg.elevation,...
            layer_data_x,'linear');
          surface = interp1(cur_echowin.eg.gps_time,...
            cur_echowin.eg.surface,...
            layer_data_x,'linear');
          physical_constants;
          for idx = 1: length(layer_data_y)
            if layer_data_y(idx) <= surface(idx)
              layer_y_curUnit(idx) = elevation(idx) - layer_data_y(idx)*c/2;
            else
              layer_y_curUnit(idx) = elevation(idx) - surface(idx)*c/2 - ...
                (layer_data_y(idx) - surface(idx))*c/(sqrt(er_ice)*2);
            end
          end
          layer_y_curUnit_manual = layer_y_curUnit(layerM_A == 1);
          layer_y_curUnit_auto = layer_y_curUnit(layerM_A == 2);
        elseif yaxis_choice == 3 % Depth/Range
          surface = interp1(cur_echowin.eg.gps_time,...
            cur_echowin.eg.surface,...
            layer_data_x,'linear');
          physical_constants;
          for idx = 1: length(layer_data_y)
            if layer_data_y(idx) <= surface(idx)
              layer_y_curUnit(idx) = layer_data_y(idx)*c/2;
            else
              layer_y_curUnit(idx) = surface(idx)*c/2 + (layer_data_y(idx) - surface(idx))*c/(sqrt(er_ice)*2);
            end
          end
          layer_y_curUnit_manual = layer_y_curUnit(layerM_A == 1);
          layer_y_curUnit_auto = layer_y_curUnit(layerM_A == 2);
        elseif yaxis_choice == 4 % Range bin
          layer_y_curUnit = interp1(cur_echowin.eg.time,cur_echowin.eg.image_yaxis,...
            layer_data_y,'linear');
          layer_y_curUnit_manual = layer_y_curUnit(layerM_A == 1);
          layer_y_curUnit_auto = layer_y_curUnit(layerM_A == 2);
        end
        %=======================================================================
        % Update layer handle
        % manual pts handle
        set(cur_echowin.layer_h(2*cur_layer-1),{'XData','YData'},{layer_x_curUnit_manual,layer_y_curUnit_manual});
        % auto pts handle
        set(cur_echowin.layer_h(2*cur_layer),{'XData','YData'},{layer_x_curUnit_auto,layer_y_curUnit_auto});
        % Update quality handles
        good_idx_manual = layerM_A == 1 & ~logical(layerQual-1);
        good_idx_auto = layerM_A == 2 & ~logical(layerQual-1);
        good_x_manual = layer_x_curUnit(good_idx_manual);
        good_x_auto = layer_x_curUnit(good_idx_auto);
        good_y_manual = layer_y_curUnit(good_idx_manual);
        good_y_auto = layer_y_curUnit(good_idx_auto);
        moderate_idx_manual = layerM_A == 1 & ~logical(layerQual-2);
        moderate_idx_auto = layerM_A == 2 & ~logical(layerQual-2);
        moderate_x_manual = layer_x_curUnit(moderate_idx_manual);
        moderate_x_auto = layer_x_curUnit(moderate_idx_auto);
        moderate_y_manual = layer_y_curUnit(moderate_idx_manual);
        moderate_y_auto = layer_y_curUnit(moderate_idx_auto);
        bad_idx_manual = layerM_A == 1 & ~logical(layerQual-3);
        bad_idx_auto = layerM_A == 2 & ~logical(layerQual-3);
        bad_x_manual = layer_x_curUnit(bad_idx_manual);
        bad_x_auto = layer_x_curUnit(bad_idx_auto);
        bad_y_manual = layer_y_curUnit(bad_idx_manual);
        bad_y_auto = layer_y_curUnit(bad_idx_auto);
        set(cur_echowin.quality_h(6*(cur_layer-1)+1),{'XData','YData'},{good_x_manual,good_y_manual});
        set(cur_echowin.quality_h(6*(cur_layer-1)+2),{'XData','YData'},{good_x_auto,good_y_auto});
        set(cur_echowin.quality_h(6*(cur_layer-1)+3),{'XData','YData'},{moderate_x_manual,moderate_y_manual});
        set(cur_echowin.quality_h(6*(cur_layer-1)+4),{'XData','YData'},{moderate_x_auto,moderate_y_auto});
        set(cur_echowin.quality_h(6*(cur_layer-1)+5),{'XData','YData'},{bad_x_manual,bad_y_manual});
        set(cur_echowin.quality_h(6*(cur_layer-1)+6),{'XData','YData'},{bad_x_auto,bad_y_auto});
        %=========================================================================
        % Update layer data vector
        cur_echowin.eg.layer.y{cur_layer} = layer_data_y;
        cur_echowin.eg.layer.x{cur_layer} = layer_data_x;
        cur_echowin.eg.layer.qual{cur_layer} = layerQual;
        cur_echowin.eg.layer.type{cur_layer} = layerM_A;
        cur_echowin.eg.layer.x_curUnit{cur_layer} = layer_x_curUnit;
        cur_echowin.eg.layer.y_curUnit{cur_layer} = layer_y_curUnit;
      end
    end
    cur_echowin.undo_stack_ptr = obj.pointer;
  end
  %%=========================================================================
  % indicate whether each echowin has been modified
  if obj.echowin_list(echowin_idx).isvalid && ~obj.ismodified()
    set(obj.echowin_list(echowin_idx).left_panel.savePB,'String','Save Layer');
  elseif obj.echowin_list(echowin_idx).isvalid && obj.ismodified()
    set(obj.echowin_list(echowin_idx).left_panel.savePB,'String','Save Layer*');
  end
end
return
