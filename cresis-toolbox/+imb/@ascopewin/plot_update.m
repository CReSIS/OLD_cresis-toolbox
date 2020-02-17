function plot_update(obj)
% plot_update(obj)

% Update plot based on selection
set(obj.h_ascope(obj.ascope.selected),'Color','red');
set(obj.h_ascope(~obj.ascope.selected),'Color','black');

% Update plot based on visibility
set(obj.h_ascope(obj.ascope.visible),'Visible','on');
set(obj.h_ascope(~obj.ascope.visible),'Visible','off');

% Update list box entries
LB_strings = cell(1,length(obj.ascope.sys));
obj.xlims = [inf -inf];
obj.ylims = [inf -inf];
for idx = 1:length(obj.ascope.sys)
  if obj.ascope.visible(idx)
    LB_strings{idx} = sprintf('%s %s', ...
      obj.ascope.frm_str{idx},datestr(epoch_to_datenum(obj.ascope.gps_time(idx)),'HH:mm:SS'));
  else
    LB_strings{idx} = sprintf('<HTML><FONT color="red">%s %s</FONT></HTML>', ...
      obj.ascope.frm_str{idx},datestr(epoch_to_datenum(obj.ascope.gps_time(idx)),'HH:mm:SS'));
  end
  uistack(obj.h_ascope(idx),'bottom');
  if obj.ascope.xlims{idx}(1) < obj.xlims(1)
    obj.xlims(1) = obj.ascope.xlims{idx}(1);
  end
  if obj.ascope.xlims{idx}(2) > obj.xlims(2)
    obj.xlims(2) = obj.ascope.xlims{idx}(2);
  end
  if obj.ascope.ylims{idx}(1) < obj.ylims(1)
    obj.ylims(1) = obj.ascope.ylims{idx}(1);
  end
  if obj.ascope.ylims{idx}(2) > obj.ylims(2)
    obj.ylims(2) = obj.ascope.ylims{idx}(2);
  end
end
set(obj.left_panel.ascopeLB,'String',LB_strings);

notify(obj,'StateChange');

end
