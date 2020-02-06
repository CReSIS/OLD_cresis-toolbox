function load_crossovers(obj,source,event)
% echowin.load_crossovers(obj,source,event)
%
% Load crossover information from database and update crossover plot handles

if obj.crossovers_en
  fprintf(' Loading crossovers from database (%s)\n', datestr(now,'HH:MM:SS'));
  %% Crossover management
  ops_param = struct('properties',[]);
  ops_param.properties.location = obj.eg.cur_sel.location;
  ops_param.properties.lyr_name = obj.eg.layers.layer_names;
  ops_param.properties.frame = obj.eg.frame_names(obj.eg.frame_idxs);
  ops_param.properties.segment_id = ones(size(ops_param.properties.frame)) ...
    *double(obj.eg.cur_sel.segment_id);
  
  [status,data] = opsGetCrossovers(obj.eg.system,ops_param);
  
  % For each cross over, find its gps_time
  data.properties.gps_time = [];
  for idx = 1:length(data.properties.source_point_path_id)
    match_idx = find(data.properties.source_point_path_id(idx) == obj.eg.map_id);
    if isempty(match_idx)
      warning('DEBUG CODE: Should never happen. Crossover returned that is not part of the loaded frames, find closest point');
      [~,match_idx] = min(double(obj.eg.map_id) - double(data.properties.source_point_path_id(idx)));
      match_idx
    end
    data.properties.gps_time(idx,1) = obj.eg.map_gps_time(match_idx);
  end
  fprintf('  Done (%s)\n', datestr(now,'HH:MM:SS'));
  
else
  data.properties.source_point_path_id = [];
  data.properties.cross_point_path_id = [];
  data.properties.source_elev = [];
  data.properties.cross_elev = [];
  data.properties.layer_id = [];
  data.properties.frame_name = [];
  data.properties.twtt = [];
  data.properties.angle = [];
  data.properties.abs_error = [];
  data.properties.gps_time = [];
  
end

% Sort results (group cross overs together)
obj.eg.crossovers.gui.set_crossovers(data.properties);

% Merge structs to add these fields to crossover:
% source_point_path_id = integer array
% cross_point_path_id = integer array
% source_elev = double array
% cross_elev = double array
% layer_id = integer array
% frame_name = integer array
% properties.twtt = double array
% properties.angle = double array
% properties.abs_error = double array
obj.eg.crossovers = merge_structs(obj.eg.crossovers,data.properties);

obj.plot_crossovers();

end

