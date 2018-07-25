function [ OPS_Surface, OPS_Bottom, OPS_data, OPS_crossover_data, total_bins] = setup_OPS(params, param_override)
%PRE_LOAD_OPS 
% Because OPS server does not support full parallel file access, we make
% this helper function to pre-load all the files needed (e.g. crossovers)
% for the Viterbi 2D ice-bottom detection.

% We do this also because these files are going to be used repeatedly from
% each different tasks run in the cluster

OPS_Surface = containers.Map('KeyType','int32', 'ValueType', 'any');
OPS_Bottom = containers.Map('KeyType','int32', 'ValueType', 'any');
OPS_frame = containers.Map('KeyType','int32', 'ValueType', 'any');
OPS_data = containers.Map('KeyType','int32', 'ValueType', 'any');
OPS_crossover_data = containers.Map('KeyType','int32', 'ValueType', 'any');

for param_idx = 1:length(params)
  param = params(param_idx);  
  % comment this out to run entire season
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
%     fprintf('passing %d \n', param_idx);
    continue;
  end
  % end of comment this out
  
  param = merge_structs(param,param_override);
  layer_params.name   = 'surface';
  layer_params.source = 'ops';

  try
    OPS_Surface(param_idx) = opsLoadLayers(param,layer_params);    
  catch ME
%       warning(ME.getReport);
    fprintf('!! data not loaded with param %d\n', param_idx);
    continue;
  end
  fprintf('Surface loaded with param %d\n', param_idx);
  
  layer_params.name   = 'bottom';
  OPS_Bottom(param_idx) = opsLoadLayers(param,layer_params);
  fprintf('Reference loaded with param %d\n', param_idx);   
  
  opsAuthenticate(param,false);
  layer_name                   = 'bottom';
  sys                          = ct_output_dir(param.radar_name);
  ops_param                    = struct('properties',[]);
  ops_param.properties.season  = param.season_name;
  ops_param.properties.segment = param.day_seg;
  [status,OPS_frame(param_idx)]          = opsGetSegmentInfo(sys,ops_param);

  ops_param = struct('properties',[]);
  ops_param.properties.location = param.post.ops.location;
  ops_param.properties.season = param.season_name;
  ops_param.properties.start_gps_time = OPS_frame(param_idx).properties.start_gps_time(1);
  ops_param.properties.stop_gps_time = OPS_frame(param_idx).properties.stop_gps_time(end);
  ops_param.properties.nativeGeom = true;
  [~,OPS_data(param_idx)] = opsGetPath(sys,ops_param);

  query = sprintf('SELECT rds_segments.id FROM rds_seasons,rds_segments where rds_seasons.name=''%s'' and rds_seasons.id=rds_segments.season_id and rds_segments.name=''%s''',param.season_name,param.day_seg);
  [status,tables] = opsQuery(query);
  segment_id = tables{1};
  ops_param                       = struct('properties',[]);
  ops_param.properties.location   = param.post.ops.location;
  ops_param.properties.lyr_name   = layer_name;
  ops_param.properties.frame      = OPS_frame(param_idx).properties.frame;
  ops_param.properties.segment_id = ones(size(ops_param.properties.frame)) ...
    *double(segment_id);  
  [status,OPS_crossover_data(param_idx)] = opsGetCrossovers(sys,ops_param);  
end

total_bins = 0;
for idx = cell2mat(keys(OPS_Surface))
  total_bins = total_bins + size(OPS_Surface(idx).gps_time,2);
end

end

