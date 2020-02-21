function layerLB_callback(obj,source,event)
% layerLB_callback(obj,source,event)

% Ensure focus stays on figure to prevent hotkeys registering with this
% uicontrol.
% uicontrol(obj.right_panel.status_panel.statusText);

if strcmpi(get(obj.h_fig,'SelectionType'),'open')
  val = get(source,'Value');
  if ~isempty(val)
    val = val(1);
    ylims_new = [min(obj.eg.layers.y_curUnit{val}) max(obj.eg.layers.y_curUnit{val})];
    dy = obj.eg.image_yaxis(2)-obj.eg.image_yaxis(1);
    if diff(ylims_new)/dy < 20
      ylims = diff(ylim(obj.h_axes))*[-0.5 0.5] + mean(ylims_new);
    else
      ylims = diff(ylims_new)*[-2 2] + mean(ylims_new);
    end
    
    % Convert x_min, x_max to GPS time
    xlims = xlim(obj.h_axes);
    xlims = interp1(obj.eg.image_xaxis,obj.eg.image_gps_time,xlims,'linear','extrap');
    
    % Redraw echogram imagesc with new axis
    obj.redraw(xlims(1),xlims(2),ylims(1),ylims(2),struct('clipped',1,'ylim_force',true));
  end
  
else
  val = get(source,'Value');
  
  obj.eg.layers.selected_layers(:)=false;
  obj.eg.layers.selected_layers(val)=true;
  
  % Update plot based on selection
  obj.set_visibility();
end