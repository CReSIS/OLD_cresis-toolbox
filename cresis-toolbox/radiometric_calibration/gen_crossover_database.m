%% OM

% use ops matlab/crossovers
% start with a small block in a grid. 
% (for+for for traverseing through entire greenland later)

% focus on crossovers for now 
%(natural targets takes time, do it post-comprehensive)

% generate the complete db first, then categorize type of crossovers
% close flight paths, intersection angles>60 degrees ish or etc

%%

 sys = 'rds';

  % Get the segment ID
  ops_param = [];
  ops_param.properties.search_str = '2014';
  ops_param.properties.location = 'arctic';
  ops_param.properties.season = '2014_Greenland_P3';
%   [status,ops_data_segment] = opsGetFrameSearch(sys,ops_param)
  
%%
  ops_param = [];
  ops_param.properties.location = 'arctic';
  ops_param.properties.lyr_name = 'surface';
  ops_param.properties.frame = {};
  ops_param.properties.segment_id = {'2014'};
  ops_param.properties.segment_id = {};
  
%   for frm = frm_list
%     ops_param.properties.frame{end+1} = [param.day_seg '_' sprintf('%03d',frm)];
%     ops_param.properties.segment_id{end+1} = ops_data_segment.properties.segment_id;
%   end
  [status,ops_crossovers] = opsGetCrossovers(sys,ops_param);