function memory(obj,vals)
% memory(obj,vals)

if strcmp(get(obj.h_fig,'visible'),'off')
  set(obj.h_fig,'visible','on');
  figure(obj.h_fig);
  
else
  if nargin == 1
    vals = find(~isnan(obj.ascope.echowin));
  end
  for new_idx = vals
    if ~isnan(obj.ascope.echowin(new_idx))
      ascope.echowin = NaN;
      ascope.sys = obj.ascope.sys{new_idx};
      ascope.season_name = obj.ascope.season_name{new_idx};
      ascope.frm_str = obj.ascope.frm_str{new_idx};
      ascope.gps_time = obj.ascope.gps_time(new_idx);
      ascope.twtt = obj.ascope.twtt{new_idx};
      ascope.data = obj.ascope.data{new_idx};
      ascope.surf_twtt = obj.ascope.surf_twtt(new_idx);
      ascope.lat = obj.ascope.lat(new_idx);
      ascope.lon = obj.ascope.lon(new_idx);
      ascope.target_twtt = obj.ascope.target_twtt(new_idx);
      ascope.xlims = obj.ascope.xlims{new_idx};
      ascope.ylims = obj.ascope.ylims{new_idx};
      
      obj.update_ascope(ascope);
    end
  end
end
