function load_layers_init(obj)
% echowin.load_layers_init(obj)
%
% Load layer information from database and update layer plot handles
 %% OPS: Preallocating layer arrays
if strcmpi(obj.eg.layer_source,'OPS')
  for idx = 1:length(obj.eg.layer_id)
    obj.eg.layer.x{idx} = double(obj.eg.map_gps_time); % gps-time
    obj.eg.layer.y{idx} = NaN*zeros(size(obj.eg.map_id)); % twtt
    obj.eg.layer.qual{idx} = NaN*zeros(size(obj.eg.map_id)); % integer 1-3
    obj.eg.layer.type{idx} = NaN*zeros(size(obj.eg.map_id)); % this is either 1 (manual) or 2 (auto)
  end
  
  %% LayerData: Preallocating layer arrays
else
  for idx = 1:length(obj.eg.layer_id)
      obj.eg.layer.x{idx} = []; %gps time
      obj.eg.layer.y{idx} = []; % twtt
      obj.eg.layer.qual{idx} = []; % integer 1-3
      obj.eg.layer.type{idx} = []; % this is either 1 (manual) or 2 (auto)
  end
end
%% Plot layers
delete(obj.layer_h);
delete(obj.quality_h);
%% WARNING: DO NOT IMPLEMENT WITH SCATTER... TOO SLOW RENDERING
layer_data_x = obj.eg.layer.x;
for idx = 1:length(layer_data_x)
  % Manual points (plot this way to handle empty XData or YData
  obj.layer_h(2*(idx-1)+1) = plot(obj.right_panel.axes.handle,NaN,NaN,'bx');
  % Auto and manual points
  obj.layer_h(2*(idx-1)+2) = plot(obj.right_panel.axes.handle, ...
    NaN,NaN,'b--');

  % Good manual points (plot this way to handle empty XData or YData
  obj.quality_h(6*(idx-1)+1) = plot(obj.right_panel.axes.handle,1,1,'gx');
  
  % Good auto points (plot this way to handle empty XData or YData
  obj.quality_h(6*(idx-1)+2) = plot(obj.right_panel.axes.handle,1,1,'g--');

  % Moderate manual points (plot this way to handle empty XData or YData
  obj.quality_h(6*(idx-1)+3) = plot(obj.right_panel.axes.handle,1,1,'yx');
  
  % Moderate auto points (plot this way to handle empty XData or YData
  obj.quality_h(6*(idx-1)+4) = plot(obj.right_panel.axes.handle,1,1,'y--');

  % Derived manual points (plot this way to handle empty XData or YData
  obj.quality_h(6*(idx-1)+5) = plot(obj.right_panel.axes.handle,1,1,'rx');
  
  % Derived auto points (plot this way to handle empty XData or YData
  obj.quality_h(6*(idx-1)+6) = plot(obj.right_panel.axes.handle,1,1,'r--');
end

end
