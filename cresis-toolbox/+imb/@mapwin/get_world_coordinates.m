% Loading world coordinates from season layerdata files for the selected
% seasons

function [wc_xs, wc_ys] = get_world_coordinates(obj)
  wc_xs = [];
  wc_ys = [];
  
  % Looping through the seasons
  for season_idx = 1:length(obj.cur_map_pref_settings.seasons)
    
    %Loading the season layerdata files
    load(char(strcat(obj.cur_map_pref_settings.system,'_param_', obj.cur_map_pref_settings.seasons{season_idx},'.mat')));
    wc_xs = [wc_xs wc_x];
    wc_ys = [wc_ys wc_y];
  end
end