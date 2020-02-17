function xaxisPM_callback(obj,hObj,event)
% echowin.xaxisPM_callback(obj,hObj,event)
%
% Callback for x-axis popup menu

physical_constants

% Ensure focus stays on figure to prevent hotkeys registering with this
% uicontrol.
uicontrol(obj.right_panel.status_panel.statusText);

xlims = xlim(obj.h_axes);
% Convert all x-limits to twtt
if obj.cur_xaxis == 2
  xlims(xlims<0) = xlims(xlims<0)*2/c*1e6;
  xlims(xlims>0) = xlims(xlims>0)*2*sqrt(er_ice)/c*1e6;
  obj.xlims(obj.xlims<0) = obj.xlims(obj.xlims<0)*2/c*1e6;
  obj.xlims(obj.xlims>0) = obj.xlims(obj.xlims>0)*2*sqrt(er_ice)/c*1e6;
end

val = get(obj.left_panel.xaxisPM,'Value');

% Update xlabel
if val == 1
  xlabel(obj.h_axes,'TWTT (us)');
elseif val == 2
  xlabel(obj.h_axes,'Depth (m)');
end

for idx = 1:length(obj.ascope.echowin)
  % Update label
  if val == 1
    xlabel(obj.h_axes,'TWTT (us)');
    set(obj.h_ascope(idx),'XData',(obj.ascope.twtt{idx} - obj.ascope.surf_twtt(idx))*1e6);
  elseif val == 2
    xlabel(obj.h_axes,'Depth (m)');
    depth = zeros(size(obj.ascope.twtt{idx}));
    above_surf = obj.ascope.twtt{idx} < obj.ascope.surf_twtt(idx);
    depth(above_surf) = (obj.ascope.twtt{idx}(above_surf) - obj.ascope.surf_twtt(idx)) * c/2;
    depth(~above_surf) = (obj.ascope.twtt{idx}(~above_surf) - obj.ascope.surf_twtt(idx)) * c/2/sqrt(er_ice);
    set(obj.h_ascope(idx),'XData',depth);
  end
end

% Convert all x-limits to depth if needed
if val == 2
  xlims(xlims<0) = xlims(xlims<0)/1e6 * c/2;
  xlims(xlims>0) = xlims(xlims>0)/1e6 * c/2/sqrt(er_ice);
  obj.xlims(obj.xlims<0) = obj.xlims(obj.xlims<0)/1e6 * c/2;
  obj.xlims(obj.xlims>0) = obj.xlims(obj.xlims>0)/1e6 * c/2/sqrt(er_ice);
end
xlim(obj.h_axes,obj.xlims);
obj.cur_xaxis = val;