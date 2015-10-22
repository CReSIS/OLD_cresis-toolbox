function plot_layers(obj)
% echowin.plot_layers(obj)
%
% Update layer plot handles

physical_constants;
% ======================================================================
%% Convert the data along the y-axis according to the units
% perform y-axis conversion (from twtt)
yaxis_choice = get(obj.left_panel.yaxisPM,'Value');
if yaxis_choice == 1 % TWTT
  % convert layerPnts_y from TWTT to current unit
  layer_y_curUnit = cell(1,length(obj.eg.layer.y));
  for idx = 1:length(obj.eg.layer.y)
    layer_y_curUnit{idx} = obj.eg.layer.y{idx}*1e6;
  end
  
elseif yaxis_choice == 2 % WGS_84 Elevation
  % Convert layerPnts_y from TWTT to WGS_84 Elevation
  layer_y_curUnit = obj.eg.layer.y;
  elevation = interp1(obj.eg.gps_time,...
    obj.eg.elevation,...
    obj.eg.layer.x{1},'linear','extrap');
  surface = interp1(obj.eg.gps_time,...
    obj.eg.surface,...
    obj.eg.layer.x{1},'linear','extrap');
  for idx = 1: length(layer_y_curUnit)
    for pnt_idx = 1:length(layer_y_curUnit{idx})
      range = min(layer_y_curUnit{idx}(pnt_idx),surface(pnt_idx))*c/2 ...
        +max(0,layer_y_curUnit{idx}(pnt_idx)-surface(pnt_idx)) * c/(sqrt(er_ice)*2);
      layer_y_curUnit{idx}(pnt_idx) = elevation(pnt_idx) - range;
    end
    layer_y_curUnit{idx}(isnan(obj.eg.layer.y{idx})) = NaN;
  end
  
elseif yaxis_choice == 3 % Depth/Range
  % convert layerPnts_y from TWTT to depth
  layer_y_curUnit = obj.eg.layer.y;
  surface = interp1(obj.eg.gps_time,...
    obj.eg.surface,...
    obj.eg.layer.x{1},'linear','extrap');
  for idx = 1:length(layer_y_curUnit)
    for pnt_idx = 1:length(layer_y_curUnit{idx})
      layer_y_curUnit{idx}(pnt_idx) = min(layer_y_curUnit{idx}(pnt_idx),surface(pnt_idx))*c/2 ...
        +max(0,layer_y_curUnit{idx}(pnt_idx)-surface(pnt_idx)) * c/(sqrt(er_ice)*2);
    end
    layer_y_curUnit{idx}(isnan(obj.eg.layer.y{idx})) = NaN;
  end
  
elseif yaxis_choice == 4 % Range bin
  % convert gCtrl layerPnts_y from TWTT to range bin
  layer_y_curUnit = obj.eg.layer.y;
  for idx = 1:length(layer_y_curUnit)
    layer_y_curUnit{idx} = interp1(obj.eg.time,...
      1:length(obj.eg.time),...
      layer_y_curUnit{idx},'linear','extrap');
  end
end

%% Plot layers
%% WARNING: DO NOT IMPLEMENT WITH SCATTER... TOO SLOW RENDERING
layer_data_x = obj.eg.layer.x;
for idx = 1:length(layer_data_x)
  % Convert x-axis units
  layer_x_curUnit = interp1(obj.eg.image_gps_time,...
    obj.eg.image_xaxis,...
    layer_data_x{idx},'linear','extrap');
  
  % get manual/auto pts (use them for layer handles)
  layer_manual = obj.eg.layer.type{idx} == 1;
  
  % Manual points (plot this way to handle empty XData or YData
  set(obj.layer_h(2*(idx-1)+1),{'XData','YData'}, ...
    {layer_x_curUnit(layer_manual),layer_y_curUnit{idx}(layer_manual)});
  % Auto and manual points
  set(obj.layer_h(2*(idx-1)+2),{'XData','YData'}, ...
    {layer_x_curUnit,layer_y_curUnit{idx}});
  
  layer_good = obj.eg.layer.qual{idx} == 1;
  layer_moderate = obj.eg.layer.qual{idx} == 2;
  layer_derived = obj.eg.layer.qual{idx} == 3;
  layer_y_curUnit_good = layer_y_curUnit{idx};
  layer_y_curUnit_good(~layer_good) = NaN;
  layer_y_curUnit_moderate= layer_y_curUnit{idx};
  layer_y_curUnit_moderate(~layer_moderate) = NaN;
  layer_y_curUnit_derived= layer_y_curUnit{idx};
  layer_y_curUnit_derived(~layer_derived) = NaN;

  % Good manual points (plot this way to handle empty XData or YData
  set(obj.quality_h(6*(idx-1)+1),{'XData','YData'}, ...
    {layer_x_curUnit(layer_good&layer_manual),layer_y_curUnit{idx}(layer_good&layer_manual)});
  
  % Good auto points (plot this way to handle empty XData or YData
  set(obj.quality_h(6*(idx-1)+2),{'XData','YData'}, ...
    {layer_x_curUnit,layer_y_curUnit_good});

  % Moderate manual points (plot this way to handle empty XData or YData
  set(obj.quality_h(6*(idx-1)+3),{'XData','YData'}, ...
    {layer_x_curUnit(layer_moderate&layer_manual),layer_y_curUnit{idx}(layer_moderate&layer_manual)});
  
  % Moderate auto points (plot this way to handle empty XData or YData
  set(obj.quality_h(6*(idx-1)+4),{'XData','YData'}, ...
    {layer_x_curUnit,layer_y_curUnit_moderate});

  % Derived manual points (plot this way to handle empty XData or YData
  set(obj.quality_h(6*(idx-1)+5),{'XData','YData'}, ...
    {layer_x_curUnit(layer_derived&layer_manual),layer_y_curUnit{idx}(layer_derived&layer_manual)});
  
  % Derived auto points (plot this way to handle empty XData or YData
  set(obj.quality_h(6*(idx-1)+6),{'XData','YData'}, ...
    {layer_x_curUnit,layer_y_curUnit_derived});

  %% Update obj.eg layers with current units
  obj.eg.layer.x_curUnit{idx} = layer_x_curUnit;
  obj.eg.layer.y_curUnit{idx} = layer_y_curUnit{idx};
end

end
