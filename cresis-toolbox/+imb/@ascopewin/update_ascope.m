function update_ascope(obj,ascope)
% update_ascope(obj,ascope)
%
% ascope.echowin
% ascope.sys
% ascope.season_name
% ascope.frm_str
% ascope.gps_time
% ascope.twtt
% ascope.data
% ascope.surf_elev
% ascope.lat
% ascope.lon
% ascope.elev

new_idx = [];
if ~isnan(ascope.echowin)
  for idx = 1:length(obj.ascope.echowin)
    if ascope.echowin == obj.ascope.echowin(idx)
      new_idx = idx;
      break;
    end
  end
end
if isempty(new_idx)
  new_idx = length(obj.ascope.echowin)+1;
  obj.h_ascope(new_idx) = plot(NaN,NaN,'r-','parent',obj.h_axes);
end

obj.ascope.echowin(new_idx) = ascope.echowin;
obj.ascope.sys{new_idx} = ascope.sys;
obj.ascope.season_name{new_idx} = ascope.season_name;
obj.ascope.frm_str{new_idx} = ascope.frm_str;
obj.ascope.gps_time(new_idx) = ascope.gps_time;
obj.ascope.twtt{new_idx} = ascope.twtt;
obj.ascope.data{new_idx} = ascope.data;
obj.ascope.surf_twtt(new_idx) = ascope.surf_twtt;
obj.ascope.lat(new_idx) = ascope.lat;
obj.ascope.lon(new_idx) = ascope.lon;
obj.ascope.target_twtt(new_idx) = ascope.target_twtt;
obj.ascope.xlims{new_idx} = [ascope.twtt([1 end]) - ascope.surf_twtt]*1e6;
obj.ascope.ylims{new_idx} = [min(10*log10(ascope.data)) max(10*log10(ascope.data))];
obj.ascope.visible(new_idx) = true;
obj.ascope.selected(new_idx) = true;
set(obj.left_panel.ascopeLB,'Value',find(obj.ascope.selected));

val = get(obj.left_panel.xaxisPM,'Value');
if val == 2
  physical_constants;
  depth = zeros(size(obj.ascope.twtt{new_idx}));
  above_surf = obj.ascope.twtt{new_idx} < obj.ascope.surf_twtt(new_idx);
  depth(above_surf) = (obj.ascope.twtt{new_idx}(above_surf) - obj.ascope.surf_twtt(new_idx)) * c/2;
  depth(~above_surf) = (obj.ascope.twtt{new_idx}(~above_surf) - obj.ascope.surf_twtt(new_idx)) * c/2/sqrt(er_ice);
  set(obj.h_ascope(new_idx),{'XData','YData','Visible','Color'}, ...
    {depth, 10*log10(obj.ascope.data{new_idx}),'on','red'});
else
  set(obj.h_ascope(new_idx),{'XData','YData','Visible','Color'}, ...
    {(obj.ascope.twtt{new_idx} - obj.ascope.surf_twtt(new_idx))*1e6, ...
    10*log10(obj.ascope.data{new_idx}), ...
    'on','red'});
end

obj.plot_update();
