% Load world coordinates from season layerdata files for the selected
% seasons

function [wc_xs, wc_ys, frms] = get_world_coordinates(obj)
  wc_xs = [];
  wc_ys = [];
  frms = [];
  
  % Looping through the seasons
  for season_idx = 1:length(obj.cur_map_pref_settings.seasons)
    
    %Loading the season layerdata files
    S = load(char(strcat('X:csarp_support\season_layerdata_files\',obj.cur_map_pref_settings.system, '_param_',obj.cur_map_pref_settings.seasons{season_idx},'_layerdata.mat')));
    [wc_x, wc_y] = imb.latlon_to_world(S.lat, S.lon);
    wc_xs = [wc_xs wc_x];
    wc_ys = [wc_ys wc_y];
    frms = [frms S.frm];
  end
end