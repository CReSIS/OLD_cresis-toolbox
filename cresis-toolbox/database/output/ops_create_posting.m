%
% OpenPolarServer Database -> layerData Posting Tool
%
% Posts data from the database to new layerData files.
%
% layerData format
%
%   struct with fields:
%     .GPS_time: [1xn double]
%     .Latitude: [1xn double]
%     .Longitude: [1xn double]
%     .Elevation: [1xn double]
%     .layerData: {1xm cell of structs}
%           .layer_points_id: [1xpts double]
%           .longitude: [1xp double]
%           .latitude: [1xp double]
%           .elevation: [1xp double]
%           .gps_time: [1xp double]
%           .twtt: [1xp double]
%           .type: [1xp double]
%           .quality: [1xp double]
%           .name: string (name of layer)
%           .description: string (description of layer)
%
% Input:
%   sys: radar system (rds,accum,snow,...)
%   param_fn: param spreadsheet path and filename
%   location: location of data ('arctic','antarctic',...)
%
% Output:
%   layerData files are written in new format with data from the database.
%
% Author: Kyle W. Purdon

%##########################################################################
% USER INPUT SECTION
%##########################################################################

sys = 'snow';
param_fn = 'C:\Users\kpurdon\Documents\scripts\params-cr1\snow_param_2011_Antarctica_DC8.xls';
location = 'antarctic';

% SET the layerdata directory:
%   'layerData' for normal location
%   'qlook' for normal location using radar echograms
%   fullfile('CSARP_post','layerData') for the posting directory
%   fullfile('CSARP_post','qlook') for the qlook directory in posting

% layer_post_directory = 'layerData';
% layer_post_directory = fullfile('CSARP_post','layerData');
layer_post_directory = fullfile('CSARP_post','qlook');


% temporarily hardcode this directory
output_base_dir = 'C:\Users\kpurdon\Documents\Projects\Other\multi-layer-post\';

%##########################################################################
% AUTOMATED SECTION
%##########################################################################

% confirm gRadar is a variable, else throw startup error.
if ~isvarname('gRadar')
  error('gRadar IS NOT A VARIABLE. RUN STARTUP.M');
end

% load the param sheet (just need good days)
params = read_param_xls(param_fn);
for param_idx = 1:length(params)
  if ~isempty(regexpi(params(param_idx).cmd.notes,'do not process')) ...
      || ~params(param_idx).cmd.generic
    continue;
  end
  
  param = params(param_idx);
  season_name = params(1).season_name;
  segment = param.day_seg;
  
  fprintf('Creating layerData for segment %s (%s) ... \n',segment,datestr(now,'HH:MM:SS'));
  cur_out_dir = ct_filename_out(param,output_base_dir,'CSARP_layerData');
  
  % get a cell of layerData filenames for the current segment
  layerdata_fns = get_filenames(ct_filename_out(param,layer_post_directory), ...
    'Data_','','.mat');
  
  % make the output directory if it does not exist
  if ~exist(cur_out_dir,'dir')
    mkdir(cur_out_dir);
  end
  
  for frame_idx = 1:length(layerdata_fns)
    
    % get information about the frame (start/stop gps, segment_id, name, season, x,y)
    [~,layer_fn_name] = fileparts(layerdata_fns{frame_idx});
    param.properties.search_str = layer_fn_name(end-14:end);
    param.properties.location = location;
    [~,info] = ops_search_frames(sys,param);
    
    % print the information to the screen
    fprintf('\tCreating layerData for frame %s (%s) ... \n',info.properties.frame,datestr(now,'HH:MM:SS'));
    
    % create the layerData output filename (temporary, just use layerData_fns)
    % fn_out = layerData_fns{frame_idx};
    fn_out = fullfile(cur_out_dir,strcat('Data_',info.properties.frame,'.mat'));
    
    % get the data for the frame (gps_time,twtt,type,quality,lyr_id)
    fprintf('\t\tQuery for frame data ...\n');
    param2.properties.location = location;
    param2.properties.season = info.properties.season;
    param2.properties.segment_id = info.properties.segment_id;
    param2.properties.start_gps_time = info.properties.start_gps_time;
    param2.properties.stop_gps_time = info.properties.stop_gps_time;
    param2.properties.lyr_name = 'all';
    [~,data] = ops_get_layer_points(sys,param2);
    
    % get the primary key and geometry for the frame (lat,lon,elev)
    fprintf('\t\tQuery for frame geometry ...\n');
    query = sprintf('SELECT layer_points_id,ST_X(layer_point),ST_Y(layer_point),ST_Z(layer_point),gps_time from %s_layer_points where segment_id=%d and gps_time BETWEEN %0.7f AND %0.7f ORDER BY gps_time;',...
      sys,param2.properties.segment_id,info.properties.start_gps_time,info.properties.stop_gps_time);
    [~,geom] = ops_query(query);
    
    % confirm that the data and geom return contain the same number of points
    if length(geom) ~= length(data.properties.gps_time)
      fprintf('GET LAYER POINTS AND GEOM QUERY LENGTH DO NOT MATCH UP. LAYERDATA WILL NOT WRITE FOR FRAME %s.',info.properties.frame);
      warning('DATA AND GEOM DIMENTION MISMATCH');
      return
    end
    
    % load the layerdata file
    fprintf('\t\tLoad stored layerData ...\n');
    layerdata = load(layerdata_fns{frame_idx});
    
    % get all of the unique layer id's in the database return
    layer_keys = unique(data.properties.lyr_id);
    
    % build the output structure
    out_layerdata = layerdata;
    out_layerdata.layerData = {};
    
    % add each layer to the output structure
    for layer_idx = 1:length(layer_keys)
      
      % get the layer information (name,description)
      layer_key = layer_keys(layer_idx);
      [~,lyr_info] = ops_query(sprintf('SELECT layer_name,description FROM %s_layers WHERE layer_id = %d AND status=''normal'';',sys,layer_key));
      
      % pre-allocate the inner geinetry structure
      layer_struct.layer_points_id = zeros(1,sum(data.properties.lyr_id==layer_key));
      layer_struct.longitude = zeros(1,sum(data.properties.lyr_id==layer_key));
      layer_struct.latitude = zeros(1,sum(data.properties.lyr_id==layer_key));
      layer_struct.elevation = zeros(1,sum(data.properties.lyr_id==layer_key));
      
      % save lat,lon,elevation from the geom return
      geom_idx = 1;
      for point_idx = 1:length(data.properties.gps_time)
          if data.properties.lyr_id(point_idx) == layer_key
              layer_struct.layer_points_id(1,geom_idx) = geom{1,point_idx};
              layer_struct.longitude(1,geom_idx) = geom{2,point_idx};
              layer_struct.latitude(1,geom_idx) = geom{3,point_idx};
              layer_struct.elevation(1,geom_idx) = geom{4,point_idx};
              geom_idx = geom_idx+1;
          end
      end
      
      % add the rest of the variables
      layer_struct.gps_time = data.properties.gps_time(data.properties.lyr_id==layer_key);
      layer_struct.twtt = data.properties.twtt(data.properties.lyr_id==layer_key);
      layer_struct.type = data.properties.type(data.properties.lyr_id==layer_key);
      layer_struct.quality = data.properties.quality(data.properties.lyr_id==layer_key);
      layer_struct.name = lyr_info{1};
      layer_struct.description = lyr_info{2};
      
      % add the output structure to the layerData cell in the root structure
      out_layerdata.layerData{end+1} = layer_struct;
      
    end
    
    % save the output file
    fprintf('\t\tSaving the output layerData ...\n');
    save(fn_out,'-struct','out_layerdata');
    
  end
end