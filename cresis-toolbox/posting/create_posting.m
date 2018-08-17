function create_posting(param, param_override)
% create_posting(param, param_override)
%
% Generalized function for posting data. Should be called from
% run_create_posting.
%
% Author: Shashanka Jagarlapudi, John Paden, Logan Smith, Theresa Stumpf
%
% See also make_layer_files, run_make_layer_files, run_picker, picker

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Input arguments check and setup
% =========================================================================

[output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name);

if ~isfield(param.post,'ops')|| isempty(param.post.ops) ...
    || isempty(param.post.ops.en)
  param.post.ops.en = 0;
end

if ~isfield(param.post,'frm_types') || isempty(param.post.frm_types)
  param.post.frm_types = {-1,0,-1,-1,-1};
end

if ~isfield(param.post,'img_type') || isempty(param.post.img_type)
  param.post.img_type = 'jpg';
end

if ~isfield(param.post,'img_dpi') || isempty(param.post.img_dpi)
  param.post.img_dpi = 150;
end

if ~isfield(param.post,'img') || isempty(param.post.img)
  % Data_YYYYMMDD_SS_FFF is image 0
  % Data_img_II_YYYYMMDD_SS_FFF is image II (where II is 1 or more)
  param.post.img = 0;
end

if ~isfield(param.post,'echo') || isempty(param.post.echo)
  param.post.echo = [];
end

if ~isfield(param.post.echo,'plot_params') || isempty(param.post.echo.plot_params)
  param.post.echo.plot_params = {'PaperPosition',[0.25 2.5 8 6]};
end

% er_ice, c = speed of light
physical_constants;

% layer_path_in: where the layer information will come from
if ~isempty(param.post.layer_dir) && param.post.ops.en == 0
  layer_path_in = ct_filename_out(param, ...
    fullfile(param.post.in_path,param.post.layer_dir), param.day_seg);
  use_data_files_for_layer = false;
else
  % An empty layer directory makes the program assume that there are no
  % layer files and it uses the data files instead
  layer_path_in = ct_filename_out(param, ...
    fullfile(param.post.in_path,param.post.data_dirs{1}), param.day_seg);
  use_data_files_for_layer = true;

  % We also need to authenticate the OPS user
  if param.post.ops.en
    opsAuthenticate([]);
  end
end

% post_path: the directory where outputs will be created
post_path = ct_filename_out(param,param.post.out_path,'CSARP_post',1);

% Do some conversions on the boolean post fields: if the field is empty,
% assign it to false
if isempty(param.post.maps_en)
  param.post.maps_en = 0;
end
if isempty(param.post.echo_en)
  param.post.echo_en = 0;
end
if isempty(param.post.layers_en)
  param.post.layers_en = 0;
end
if isempty(param.post.data_en)
  param.post.data_en = 0;
end
if isempty(param.post.csv_en)
  param.post.csv_en = 0;
end
if isempty(param.post.concat_en)
  param.post.concat_en = 0;
end
if isempty(param.post.pdf_en)
  param.post.pdf_en = 0;
end
if ~isfield(param.post,'ops') || isempty(param.post.ops.en)
  param.post.ops.en = 0;
end

if strcmp(param.post.img_type,'jpg')
  print_device = '-djpeg';
elseif strcmp(param.post.img_type,'png')
  print_device = '-dpng';
else
  error('Unsupported image type %s', param.post.img_type);
end
print_dpi = sprintf('-r%d', param.post.img_dpi);

%% Catalog layer and data files
% =========================================================================
fprintf('Catalog layer and data files %s (%s)\n', param.day_seg, datestr(now));

% =========================================================================

% Load the frames file
frames_fn = ct_filename_support(param, '', 'frames');
load(frames_fn);

% % Get a list of all the layer files
% layer_fns = get_filenames(layer_path_in,'Data_','','.mat',struct('recursive',1));
% % Layer files are named:
% %   Data_YYYYMMDD_SS_FFF.mat
% % Extract out:
% %  1. the layer file name
% %  2. frame ID
% %  3. day_seg
% % Also eliminate frames that are not in param.cmd.frms or that have
% % proc_mode set to not post.
% frm_idx = 0;
% frms = {};
% day_segs = {};
% for fn_idx = 1:length(layer_fns)
%   layer_fn = layer_fns{fn_idx};
%   [tmp layer_name] = fileparts(layer_fn);
%   if layer_name(6) == 'i'
%     % We don't process files with "_img_II" in their name
%     continue;
%   else
%     frm_id = layer_name(6:end);
%   end
%   frm = str2double(frm_id(end-2:end));
%   
%   if ~isempty(param.cmd.frms)
%     % Do just frames specified in the command worksheet (if the field
%     % is left empty we do them all)
%     if all(frm ~= param.cmd.frms)
%       continue;
%     end
%   end
%   if mod(floor(frames.proc_mode(frm)/10),10) ~= 0
%     % UUUUUUUURRRR
%     %           ^
%     %    This is the digit we want for controlling posting
%     % R must be zero to post.
%     continue;
%   end
%     
%   % Add frame to list
%   frm_idx = frm_idx + 1;
%   frms{frm_idx}.layer_fn = layer_fn;
%   frms{frm_idx}.layer_name = layer_name;
%   day_segs{frm_idx} = param.day_seg;
%   frms{frm_idx}.frm_id = frm_id;
%   frms{frm_idx}.data_fns = {};
%   frms{frm_idx}.data_dirs = {};
% end
% clear('layer_fns');

if isempty(param.cmd.frms)
  param.cmd.frms = 1:length(frames.frame_idxs);
end
% Remove frames that do not exist from param.cmd.frms list
[valid_frms,keep_idxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
if length(valid_frms) ~= length(param.cmd.frms)
  bad_mask = ones(size(param.cmd.frms));
  bad_mask(keep_idxs) = 0;
  warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
    param.cmd.frms(find(bad_mask,1)));
  param.cmd.frms = valid_frms;
end

% Build a list of frame information and eliminate frames that have
% proc_mode set to not post.
frm_idx = 0;
frms = {};
day_segs = {};
for frm = param.cmd.frms
  % Check to make sure this frame should be processed (frames file field
  % proc_mode and param.post.frm_types controls this)
  if ~ct_proc_frame(frames.proc_mode(frm),param.post.frm_types)
    continue;
  end
  
  % frm_id = string containing frame ID YYYYMMDD_SS_FFF
  frm_id = sprintf('%s_%03d', param.day_seg, frm);
  % layer_name: data filename string Data_YYYYMMDD_SS_FFF.mat or 
  %             Data_img_II_YYYYMMDD_SS_FFF.mat
  if ~use_data_files_for_layer || param.post.img == 0
    layer_name = sprintf('Data_%s', frm_id);
  else
    layer_name = sprintf('Data_img_%02d_%s', param.post.img, frm_id);
  end
  % layer_fn = full path to the layer filename
  layer_fn = fullfile(layer_path_in,[layer_name '.mat']);
  
  if ~exist(layer_fn,'file')
    error('Missing layer file %s', layer_fn);
  end
  
  % Add frame to list
  frm_idx = frm_idx + 1;
  frms{frm_idx}.layer_fn = layer_fn;
  frms{frm_idx}.layer_name = layer_name;
  day_segs{frm_idx} = param.day_seg;
  frms{frm_idx}.frm_id = frm_id;
  frms{frm_idx}.data_fns = {};
  frms{frm_idx}.data_dirs = {};
end
clear('layer_fns');

if isempty(frms)
  error('No layer files match criteria in %s', layer_path_in);
end

%%  Associate data files with layer files
% =========================================================================

for data_dir_idx = 1:length(param.post.data_dirs)
  % Get a data directory
  data_path_in = ct_filename_out(param, ...
    fullfile(param.post.in_path,param.post.data_dirs{data_dir_idx}), param.day_seg);
  if ~exist(data_path_in,'dir')
    continue;
  end
  
  % Get all the data files in the data directory and then sort out
  % just the filename and day_seg
  data_files = get_filenames(data_path_in,'Data_','','.mat',struct('recursive',1));
  clear data_files_name data_file_day_seg;
  data_files_name = {};
  data_file_frm_id = {};
  for file_idx = 1:length(data_files)
    [tmp data_files_name{file_idx}] = fileparts(data_files{file_idx});
    if data_files_name{file_idx}(6) == 'i'
      data_file_frm_id{file_idx} = data_files_name{file_idx}(13:end);
    else
      data_file_frm_id{file_idx} = data_files_name{file_idx}(6:end);
    end
  end
  
  % Find the data files in this directory associated with each layer file
  % Layer files are named:
  %   Data_YYYYMMDD_SS_FFF.mat
  % Data files are named:
  %   Data_YYYYMMDD_SS_FFF.mat
  %   Data_img_II_YYYYMMDD_SS_FFF.mat
  % The master data file (mdata) is the one which images will be made of.
  % If param.post.img == 0, then the file should be called:
  %   Data_YYYYMMDD_SS_FFF.mat
  % If no such file exists, it will look for param.post.img == 1. If 
  % param.post.img ~= 0, then the echogram name should be:
  %   Data_img_II_YYYYMMDD_SS_FFF.mat
  % where II == param.post.img.
  % If no master file exists, an exception is thrown.
  for frm_idx = 1:length(frms)
    match_idxs = strmatch(frms{frm_idx}.frm_id,data_file_frm_id);
    for match_idx = 1:length(match_idxs)
      frms{frm_idx}.data_fns{end+1} = data_files{match_idxs(match_idx)};
      frms{frm_idx}.data_dirs{end+1} = param.post.data_dirs{data_dir_idx};
    end
    % Get the master data file
    if data_dir_idx == 1
      if ~use_data_files_for_layer
        % Using layerData for layers: need to explicitly specify data file names
        if param.post.img == 0
          img_name = [frms{frm_idx}.layer_name(1:5) frms{frm_idx}.layer_name(6:end)];
        else
          img_name = [frms{frm_idx}.layer_name(1:5) sprintf('img_%02d_',param.post.img) frms{frm_idx}.layer_name(6:end)];
        end
        match_idx = strmatch(img_name,data_files_name);
        if isempty(match_idx) && param.post.img == 0
          % Legacy Support: This only happens when img==0 data file does not exist, so we
          % use img==1 in this case.
          img_name = [frms{frm_idx}.layer_name(1:5) sprintf('img_%02d_',1) frms{frm_idx}.layer_name(6:end)];
          match_idx = strmatch(img_name,data_files_name);
        end
      else
        % Using echogram data for layers: file must already be in the list
        match_idx = strmatch(frms{frm_idx}.layer_name,data_files_name);
      end
      if isempty(match_idx)
        error('Layer file without data file %s\n', frms{frm_idx}.layer_fn);
      end
      frms{frm_idx}.master_fn = data_files{match_idx};
    end
  end
end

%% Create map template
if param.post.maps_en
  fprintf(' Creating map template (%s)\n', datestr(now));

  map_param.type = param.post.map.type;
  map_param.location = param.post.map.location;
  if isfield(param.post.map,'geotiff')
    map_param.geotiff = param.post.map.geotiff;
  end
  if exist('map_info','var')
    map_param.fig_hand = map_info.fig_hand;
  else
    map_param.fig_hand = [];
  end
  map_param.resample = true;
  map_param.map_title = '';
  map_param.decimate_seg = true;
  map_param.year = str2double(param.day_seg(1:4));
  map_param.month = str2double(param.day_seg(5:6));
  map_param.day = str2double(param.day_seg(7:8));
  
  map_info = publish_map('setup',map_param);
  map_param.fig_hand = map_info.fig_hand;

  % =========================================================================
  % For each day segment, find the [x,y] extent by looking at all the frames
  % in this segment
  % =========================================================================
  [day_seg_list tmp day_segs_idx] = unique(day_segs);
  for seg_idx = 1:length(day_seg_list)
    min_x(seg_idx) = inf;
    max_x(seg_idx) = -inf;
    min_y(seg_idx) = inf;
    max_y(seg_idx) = -inf;
    frm_idxs = strmatch(day_seg_list{seg_idx},day_segs);
    day_seg_x{seg_idx} = [];
    day_seg_y{seg_idx} = [];
    for idx = 1:length(frm_idxs)
      lay = load(frms{frm_idxs(idx)}.layer_fn,'Latitude','Longitude');
      [x,y] = projfwd(map_info.proj,lay.Latitude,lay.Longitude);
      x = x/1e3;
      y = y/1e3;
      min_x(seg_idx) = min([x min_x(seg_idx)]);
      max_x(seg_idx) = max([x max_x(seg_idx)]);
      min_y(seg_idx) = min([y min_y(seg_idx)]);
      max_y(seg_idx) = max([y max_y(seg_idx)]);
      day_seg_x{seg_idx} = [x day_seg_x{seg_idx}];
      day_seg_y{seg_idx} = [y day_seg_y{seg_idx}];
    end
  end
end

%% OPS Setup
if param.post.ops.en
  %% HACK
  if length(param.post.ops.layers) == 2
    param.post.ops.layers = {'surface' 'bottom'};
  elseif length(param.post.ops.layers) == 1
    param.post.ops.layers = {'surface'};
  else
    error('length(param.post.ops.layers) must be 1 or 2');
  end
  
  %% Get the layer data for this segment
  ops_sys = output_dir;
  ops_param = [];
  ops_param.properties.location = param.post.ops.location;
  ops_param.properties.season = param.season_name;
  ops_param.properties.segment = param.day_seg;
  ops_param.properties.return_geom = 'geog';
  ops_layer = {};
  for layer_idx = 1:length(param.post.ops.layers)
    ops_param.properties.lyr_name = param.post.ops.layers{layer_idx};
    [~,ops_layer{layer_idx}] = opsGetLayerPoints(ops_sys,ops_param);
    ops_layer{layer_idx} = ops_layer{layer_idx}.properties;
  end
end

%% Main post loop
% =========================================================================

% For each frame, do the following:
%  1. produce the map
%  2. produce the echogram
%  3. produce the echogram w/ layer
%  4. copy the data files
%  5. make the csv files

for frm_idx = 1:length(frms)
  fprintf('  Posting frame %s, %d of %d (%s)\n', ...
    frms{frm_idx}.frm_id, frm_idx, length(frms), datestr(now));
  
  if param.post.ops.en
    %% Interpolate each layer for this segment onto the master layer's gps time
    
    lay = load(frms{frm_idx}.layer_fn,'GPS_time','Latitude','Longitude','Elevation');
    
    lay = opsInterpLayersToMasterGPSTime(lay,ops_layer,param.post.ops.gaps_dist);
    
    lay.Surface = lay.layerData{1}.value{2}.data;
    if length(param.post.ops.layers) == 2
      lay.Bottom = lay.layerData{2}.value{2}.data;
    elseif length(param.post.ops.layers) == 1
      lay.Bottom = NaN*zeros(size(lay.Surface));
    end
    
  else
    % Load GPS_time, Latitude, Longitude, and layerData from layer file
    if ~isempty(param.post.layer_dir)
      lay = load(frms{frm_idx}.layer_fn,'GPS_time','Latitude','Longitude','Elevation','layerData');
    else
      % An empty layer directory makes the program assume that there are no
      % layer files and it uses the data files instead
      fprintf('    %s\n', frms{frm_idx}.layer_fn);
      warning off;
      lay = load(frms{frm_idx}.layer_fn,'GPS_time','Latitude','Longitude','Elevation','Surface','Bottom');
      warning on;
    end
  end
  
  time_stamp_str = datestr(epoch_to_datenum(lay.GPS_time(1)),'HHMMSS');
  
  if param.post.maps_en || param.post.echo_en
    image_dir = fullfile(post_path,'images',day_segs{frm_idx});
    if ~exist(image_dir,'dir')
      mkdir(image_dir)
    end
  end
  if param.post.pdf_en
    pdf_dir = fullfile(post_path,'pdf',day_segs{frm_idx});
    if ~exist(pdf_dir,'dir')
      mkdir(pdf_dir)
    end
  end
  
  % =======================================================================
  % Plot Flightlines over Geotiff
  % =======================================================================
  if param.post.maps_en
    % Plot the map for each frame
    seg_idx = day_segs_idx(frm_idx);
    map_param.day_seg_x = day_seg_x{seg_idx}*1000;
    map_param.day_seg_y = day_seg_y{seg_idx}*1000;
    [frame_X,frame_Y] = projfwd(map_info.proj,lay.Latitude,lay.Longitude);
    map_param.frame_X = {frame_X};
    map_param.frame_Y = {frame_Y};
    map_info = publish_map('delete',map_param,map_info);
    map_param.map_title = sprintf('%s\n%s', param.cmd.mission_names, frms{frm_idx}.frm_id);
    map_info = publish_map('plot',map_param,map_info);
    set(map_info.h_title,'fontSize',8,'fontweight','normal');
    
    % Remove current file
    map_fn = sprintf('%s*_0maps.%s',frms{frm_idx}.frm_id,param.post.img_type);
    map_fn = fullfile(image_dir,map_fn);
    delete(map_fn);
    
    % Save file
    map_fn = sprintf('%s_0maps.%s',frms{frm_idx}.frm_id,param.post.img_type);
    map_fn = fullfile(image_dir,map_fn);
    set(map_info.fig_hand(1),'PaperUnits','inches');
    set(map_info.fig_hand(1),'PaperPosition',[0.5 0.5 10 7.5]);
    set(map_info.fig_hand(1),'PaperOrientation','Portrait');
    print(map_info.fig_hand(1),print_device,print_dpi,map_fn);
    
    if param.post.pdf_en
      % Remove current file
      map_fn = sprintf('%s*_0maps.pdf',frms{frm_idx}.frm_id);
      map_fn = fullfile(pdf_dir,map_fn);
      delete(map_fn);
      
      map_fn = sprintf('%s_%s_0maps.pdf',frms{frm_idx}.frm_id,time_stamp_str);
      map_fn = fullfile(pdf_dir,map_fn);
      set(map_info.fig_hand(1),'PaperOrientation','Landscape');
      
      set(map_info.fig_hand(1),'PaperUnits','inches');
      set(map_info.fig_hand(1),'PaperPosition',[0.5 0.5 10 7.5]);
      saveas(map_info.fig_hand(1),map_fn);
    end
  end
  
  % =======================================================================
  % Create Surface and Bottom Variables
  % =======================================================================
  if ~param.post.ops.en
    if ~isempty(param.post.layer_dir)
      lay.Bottom  = lay.layerData{2}.value{2}.data;
      lay.Surface = lay.layerData{1}.value{2}.data;
    else
      % An empty layer directory makes the program assume that there are no
      % layer files and it uses the data files instead
      lay.Bottom = NaN*zeros(size(lay.Surface));
      lay.layerData{1}.value{2}.data = lay.Surface;
      lay.layerData{2}.value{2}.data = lay.Bottom;
      lay.layerData{1}.quality = zeros(size(lay.Surface));
      lay.layerData{2}.quality = zeros(size(lay.Surface));
    end
    lay.Thickness = lay.Bottom-lay.Surface;
    neg_idxs = find(lay.Thickness < 0 & isfinite(lay.Thickness));
    if ~isempty(neg_idxs)
      sprintf('  Negative thickness detected');
      lay.Bottom(neg_idxs) = lay.Surface(neg_idxs);
    end
    lay.Thickness = lay.Bottom-lay.Surface;
  end
  
  % =======================================================================
  % Plot Echogram and Echogram w/ Layers
  % =======================================================================
  if param.post.echo_en
    % Load master file (mdata)
    fprintf('    Loading %s\n', frms{frm_idx}.master_fn);
    mdata = load(frms{frm_idx}.master_fn);
    mdata = uncompress_echogram(mdata); % Uncompress if necessary
    
    if size(mdata.Data,2) > 5000
      % Decimate data to less than 5000 along-track samples to prevent
      % failure in pdf generation
      old_rlines = 1:size(mdata.Data,2);
      new_rlines = linspace(1,size(mdata.Data,2),5000);
      mdata.Data = interp1(old_rlines,mdata.Data.',new_rlines).';
      mdata.Latitude = interp1(old_rlines,mdata.Latitude,new_rlines);
      mdata.Longitude = interp1(old_rlines,mdata.Longitude,new_rlines);
      mdata.Elevation = interp1(old_rlines,mdata.Elevation,new_rlines);
      mdata.GPS_time = interp1(old_rlines,mdata.GPS_time,new_rlines);
    end
    
    echo_param = param.post.echo;
    if exist('echo_info','var')
      echo_param.fig_hand = echo_info.fig_hand;
    else
      echo_param.fig_hand = [];
    end
    echo_param.num_x_tics = 6;

    echo_param.frm_id = frms{frm_idx}.frm_id;
    
    echo_param.ops.en = true;
    
    echo_info = publish_echogram(echo_param,mdata,lay);
    set(echo_info.h_surf,'Visible','off')
    set(echo_info.h_bot,'Visible','off')
    
    % Update title (handle segment wrapping from one day to the next)
    seg_year = str2double(param.day_seg(1:4));
    seg_month = str2double(param.day_seg(5:6));
    seg_day = str2double(param.day_seg(7:8));
    
    [year,month,day,hour,minute,sec] = datevec(epoch_to_datenum(lay.GPS_time(1)));
    days_offset = datenum(year,month,day) - datenum(seg_year,seg_month,seg_day);
    hour = hour + 24*days_offset;
    start_time_str = sprintf('%02d:%02d:%04.1f', hour, minute, sec);
    
    [year,month,day,hour,minute,sec] = datevec(epoch_to_datenum(lay.GPS_time(end)));
    days_offset = datenum(year,month,day) - datenum(seg_year,seg_month,seg_day);
    hour = hour + 24*days_offset;
    stop_time_str = sprintf('%02d:%02d:%04.1f', hour, minute, sec);
    
    echo_info.h_title = title(echo_info.ah_echo,sprintf('%s %s: "%s"  %s: %s to %s GPS', output_dir, param.season_name, ...
      param.cmd.mission_names, frms{frm_idx}.frm_id, start_time_str, ...
      stop_time_str),'Interpreter','none','FontWeight','normal','FontSize',9);
    
    % Remove current file
    echo_fn = sprintf('%s*_1echo.%s',frms{frm_idx}.frm_id,param.post.img_type);
    echo_fn = fullfile(image_dir,echo_fn);
    delete(echo_fn);
    
    echo_fn = sprintf('%s_1echo.%s',frms{frm_idx}.frm_id,param.post.img_type);
    echo_fn = fullfile(image_dir,echo_fn);
    set(echo_info.fig_hand(1),'PaperUnits','inches');
    set(echo_info.fig_hand(1),param.post.echo.plot_params{1},param.post.echo.plot_params{2});
    set(echo_info.fig_hand(1),'PaperOrientation','Portrait');
    fprintf('    Saving output %s\n', echo_fn);
    print(echo_info.fig_hand(1),print_device,print_dpi,echo_fn);
    
    if param.post.pdf_en
      % Remove current file
      echo_fn = sprintf('%s*_1echo.pdf',frms{frm_idx}.frm_id);
      echo_fn = fullfile(pdf_dir,echo_fn);
      delete(echo_fn);
    
      echo_fn = sprintf('%s_%s_1echo.pdf',frms{frm_idx}.frm_id,time_stamp_str);
      echo_fn = fullfile(pdf_dir,echo_fn);
      set(echo_info.fig_hand(1),'PaperOrientation','Landscape');
      set(echo_info.fig_hand(1),'PaperUnits','inches');
      set(echo_info.fig_hand(1),'PaperPosition',[0.5 0.5 10 7.5]);
      %print(echo_info.fig_hand(1),'-depsc','-r72',echo_fn);
      saveas(echo_info.fig_hand(1),echo_fn);
    end
    % =====================================================================
    % Plot echogram w/ layers
    if param.post.layers_en
      set(echo_info.h_surf,'Visible','on')
      set(echo_info.h_bot,'Visible','on')
      
      % Remove current file
      echo_fn = sprintf('%s*_2echo_picks.%s',frms{frm_idx}.frm_id,param.post.img_type);
      echo_fn = fullfile(image_dir,echo_fn);
      delete(echo_fn);
      
      echo_fn = sprintf('%s_2echo_picks.%s',frms{frm_idx}.frm_id,param.post.img_type);
      echo_fn = fullfile(image_dir,echo_fn);
      set(echo_info.fig_hand(1),param.post.echo.plot_params{1},param.post.echo.plot_params{2});
      set(echo_info.fig_hand(1),'PaperOrientation','Portrait');
      print(echo_info.fig_hand(1),print_device,print_dpi,echo_fn);
      
      if param.post.pdf_en
        % Remove current file
        echo_fn = sprintf('%s*_2echo_picks.pdf',frms{frm_idx}.frm_id);
        echo_fn = fullfile(pdf_dir,echo_fn);
        delete(echo_fn);

        echo_fn = sprintf('%s_%s_2echo_picks.pdf',frms{frm_idx}.frm_id,time_stamp_str);
        echo_fn = fullfile(pdf_dir,echo_fn);
        set(echo_info.fig_hand(1),'PaperOrientation','Landscape');
        set(echo_info.fig_hand(1),'PaperUnits','inches');
        set(echo_info.fig_hand(1),'PaperPosition',[0.5 0.5 10 7.5]);
        %print(echo_info.fig_hand(1),'-depsc','-r72',echo_fn);
        saveas(echo_info.fig_hand(1),echo_fn);
      end
    end
  end
  
  % Create .csv files
  if param.post.csv_en
    csv_dir = fullfile(post_path,'csv',day_segs{frm_idx});
    if ~exist(csv_dir,'dir')
      mkdir(csv_dir)
    end
    % Remove current file
    csv_fn = sprintf('Data_%s*.csv',frms{frm_idx}.frm_id);
    csv_fn = fullfile(csv_dir,csv_fn);
    delete(csv_fn);
    % Create new filename
    csv_fn = sprintf('Data_%s_%s.csv',frms{frm_idx}.frm_id,time_stamp_str);
    csv_fn = fullfile(csv_dir,csv_fn);
    
    csv_good_dir = fullfile(post_path,'csv_good',day_segs{frm_idx});
    if ~exist(csv_good_dir,'dir')
      mkdir(csv_good_dir)
    end
    % Remove current file
    csv_good_fn = sprintf('Data_%s*.csv',frms{frm_idx}.frm_id);
    csv_good_fn = fullfile(csv_good_dir,csv_good_fn);
    delete(csv_good_fn);
    % Create new filename
    csv_good_fn = sprintf('Data_%s_%s.csv',frms{frm_idx}.frm_id,time_stamp_str);
    csv_good_fn = fullfile(csv_good_dir,csv_good_fn);
    
    PThickness = (lay.Bottom-lay.Surface)*c/2/sqrt(param.post.echo.er_ice);
    PSurface = lay.Surface*c/2;
    PBottom = PThickness+PSurface;
    
    PThickness(~isfinite(PThickness)) = -9999;
    PBottom(~isfinite(PBottom)) = -9999;
    PSurface(~isfinite(PSurface)) = -9999;
    % Quality of thickness is the lowest/worst confidence level of the surface
    % and bottom picks which each have their own quality level
    % (higher quality numbers mean lower confidence, so we take the max)
    Quality = max(lay.layerData{1}.quality,lay.layerData{2}.quality(1));

    % Compute seconds of day relative to the data frame ID
    UTC_time = lay.GPS_time - utc_leap_seconds(lay.GPS_time(1));
    UTC_time_start = datenum_to_epoch(datenum(str2double(frms{frm_idx}.frm_id(1:4)), ...
      str2double(frms{frm_idx}.frm_id(5:6)), ...
      str2double(frms{frm_idx}.frm_id(7:8))));
    UTC_time = UTC_time - UTC_time_start;
    
    bad_idxs = isnan(lay.Latitude) | isnan(lay.Longitude) | isnan(lay.Elevation) | isnan(UTC_time);
    if any(bad_idxs)
      warning('There are NaN in the data: suggest correcting and rerunning');
      keyboard;
      lay.Latitude = lay.Latitude(~bad_idxs);
      lay.Longitude = lay.Longitude(~bad_idxs);
      lay.Elevation = lay.Elevation(~bad_idxs);
      PThickness = PThickness(~bad_idxs);
      PBottom = PBottom(~bad_idxs);
      PSurface = PSurface(~bad_idxs);
      UTC_time = UTC_time(~bad_idxs);
      Quality = Quality(~bad_idxs);
    end
    
    fid_csv = fopen(csv_fn,'w');
    fid_csv_good = fopen(csv_good_fn,'w');
    fprintf(fid_csv,'%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
      'LAT','LON','UTCTIMESOD','THICK','ELEVATION','FRAME','SURFACE','BOTTOM','QUALITY');
    fprintf(fid_csv_good,'%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
      'LAT','LON','UTCTIMESOD','THICK','ELEVATION','FRAME','SURFACE','BOTTOM','QUALITY');
    for txt_idx = 1:length(lay.Latitude)
      frm_id_num = [frms{frm_idx}.frm_id(1:8) frms{frm_idx}.frm_id(10:11) frms{frm_idx}.frm_id(13:15)];
      fprintf(fid_csv,'%2.6f,%2.6f,%5.4f,%6.2f,%4.4f,%s,%6.2f,%6.2f,%01d\n',...
        lay.Latitude(txt_idx),lay.Longitude(txt_idx),...
        UTC_time(txt_idx),PThickness(txt_idx),...
        lay.Elevation(txt_idx),frm_id_num,PSurface(txt_idx),PBottom(txt_idx),Quality(txt_idx));
      if PThickness(txt_idx) > -9999
          if isnan(lay.Latitude(txt_idx)) || isnan(lay.Longitude(txt_idx))
              fprintf('NaN values found in latitdue or logitude at frame: %3d,continue after debug the error\n',frm_id_num);
              keyboard
          end
          fprintf(fid_csv_good,'%2.6f,%2.6f,%5.4f,%6.2f,%4.4f,%s,%6.2f,%6.2f,%01d\n',...
          lay.Latitude(txt_idx),lay.Longitude(txt_idx),...
          UTC_time(txt_idx),PThickness(txt_idx),...
          lay.Elevation(txt_idx),frm_id_num,PSurface(txt_idx),PBottom(txt_idx),Quality(txt_idx));
      end
    end
    fclose(fid_csv);
    fclose(fid_csv_good);
  end
  
  % Create layer files (or copy if not using OPS)
  if param.post.layers_en
    if param.post.ops.en
      %% Create layer data file
      % Copy layer file
      layer_fn = frms{frm_idx}.layer_fn;
      [tmp layer_name] = fileparts(layer_fn);
      %     layer_dir = sprintf('CSARP_%s', param.post.layer_dir);
      layer_out_dir = fullfile(post_path,sprintf('CSARP_layerData'),day_segs{frm_idx});
      if ~exist(layer_out_dir,'dir')
        mkdir(layer_out_dir)
      end
      layer_out_fn = fullfile(layer_out_dir,[layer_name '.mat']);
      save(layer_out_fn,'-struct','lay','GPS_time','Latitude','Longitude','Elevation','layerData');
    else
      % Copy layer file
      layer_fn = frms{frm_idx}.layer_fn;
      [tmp layer_name] = fileparts(layer_fn);
      %     layer_dir = sprintf('CSARP_%s', param.post.layer_dir);
      layer_out_dir = fullfile(post_path,sprintf('CSARP_layerData'),day_segs{frm_idx});
      if ~exist(layer_out_dir,'dir')
        mkdir(layer_out_dir)
      end
      layer_out_fn = fullfile(layer_out_dir,[layer_name '.mat']);
      copyfile(layer_fn, layer_out_fn);
    end
  end
    
  % Copy data files
  if param.post.data_en
    % Copy data files
    for data_idx = 1:length(frms{frm_idx}.data_fns)
      data_fn = frms{frm_idx}.data_fns{data_idx};
      [tmp data_name] = fileparts(data_fn);
      data_dir = sprintf('CSARP_%s', frms{frm_idx}.data_dirs{data_idx});
      data_out_dir = fullfile(post_path,data_dir,day_segs{frm_idx});
      if ~exist(data_out_dir,'dir')
        mkdir(data_out_dir)
      end
      data_out_fn = fullfile(data_out_dir,[data_name '.mat']);
      
      copyfile(data_fn, data_out_fn);
      
      if 0
        % DATA FILES SHOULD NOT BE UPDATED! They represent the surface
        % and bottom used during processing.
        
        % Interpolate onto the data files GPS time
        data_file = load(data_out_fn,'Time','GPS_time','Elevation');
        
        % Create fast-time correction vector
        %   The layer file may contain a different elevation profile
        %   than the data file. The surface and bottom need to be adjusted
        %   to account for these differences.
        elev_interp    = interp1(lay.GPS_time,lay.Elevation,data_file.GPS_time,'linear','extrap');
        fast_time_correction = (data_file.Elevation - elev_interp)/(c/2);
        
        warning('off','MATLAB:interp1:NaNinY');
        Surface = interp1(lay.GPS_time,lay.Surface,data_file.GPS_time,'linear','extrap') ...
          + fast_time_correction;
        Bottom = interp1(lay.GPS_time,lay.Bottom,data_file.GPS_time,'linear','extrap') ...
          + fast_time_correction;
        warning('on','MATLAB:interp1:NaNinY');
        save(data_out_fn,'Bottom','Surface','-APPEND');
      end
    end
  end
end


% =======================================================================
% Create by-segment and by-season concatenated and browse files
%   CSV and KML formats
% =======================================================================

if param.post.concat_en
  fprintf(' Creating csv and kml files (%s)\n', datestr(now));
  
  csv_dir = fullfile(post_path,'csv',param.day_seg);
  kml_base_dir = fullfile(post_path,'kml');
  
  [csv_dir_path csv_dir_name] = fileparts(csv_dir);
  out_fn = fullfile(csv_dir_path,sprintf('Data_%s.csv',csv_dir_name));
  concatenate_thickness_files(csv_dir,'*.csv',out_fn,',');
  
  % Create KML browse files for each segment
  % Extract day_seg from filename
  in_fn = out_fn;
  [in_fn_dir in_fn_name] = fileparts(in_fn);
  kml_out_fn = fullfile(kml_base_dir, ['Browse_' in_fn_name '.kml']);
  day_seg = in_fn_name(6:end);
  kml_write_cresis(in_fn, kml_out_fn, day_seg,'segment',[inf 40]);
  
  % Repeat for csv_good and kml_good
  csv_dir = fullfile(post_path,'csv_good',param.day_seg);
  kml_base_dir = fullfile(post_path,'kml_good');
  
  [csv_dir_path csv_dir_name] = fileparts(csv_dir);
  out_fn = fullfile(csv_dir_path,sprintf('Data_%s.csv',csv_dir_name));
  concatenate_thickness_files(csv_dir,'*.csv',out_fn,',');
  
  % Create KML browse files for each segment
  % Extract day_seg from filename
  in_fn = out_fn;
  [in_fn_dir in_fn_name] = fileparts(in_fn);
  kml_out_fn = fullfile(kml_base_dir, ['Browse_' in_fn_name '.kml']);
  day_seg = in_fn_name(6:end);
  kml_write_cresis(in_fn, kml_out_fn, day_seg,'segment',[inf 40]);
end

% =======================================================================
% Create pdf files
% =======================================================================
if param.post.pdf_en
  fprintf(' Creating combined pdf files (%s)\n', datestr(now));
  gs_path = 'C:\Progra~1\gs\gs9.16\bin\gswin64.exe';
  if ispc && ~exist(gs_path,'file')
    warning('  Can not do this final step on a PC without ghostscript (gs version 9.02) commands. Rerun posting with just pdf enabled on the post worksheet on a linux machine and it will skip to this combine step so you don''t have to recreate all the temporary frame pdf files again.');
  else
    pdf_base_dir = fullfile(post_path,'pdf');
    % The PDF creation used to do all available directories:
    %pdf_dirs = get_filenames(pdf_base_dir,'[0-9]','','',struct('type','d'));
    % Now the PDF creation just does the current segment:
    pdf_dirs = {fullfile(pdf_base_dir,param.day_seg)};

    for dir_idx = 1:length(pdf_dirs)
      in_search_str = fullfile(pdf_dirs{dir_idx},'*.pdf');
      [pdf_dir_path pdf_dir_name] = fileparts(pdf_dirs{dir_idx});
      out_fn = fullfile(pdf_dir_path,sprintf('%s.pdf',pdf_dir_name));
      fprintf('  Creating pdf %s\n', out_fn);
      if ispc
        pdf_fns = get_filenames(pdf_dirs{dir_idx},'','','.pdf');
        pdf_fns = sprintf('%s ', pdf_fns{:});
        sys_cmd = sprintf('%s -dNOPAUSE -dBATCH -dSAFER -sOutputFile=%s -sDEVICE=pdfwrite -f %s', ...
          gs_path, out_fn, pdf_fns);
      else
        sys_cmd = sprintf('gs -dNOPAUSE -dBATCH -dSAFER -sOutputFile=%s -sDEVICE=pdfwrite -f %s </dev/null', ...
          out_fn, in_search_str);
      end
      [status,result] = system(sys_cmd);
      if status > 1
        warning('pdf creation may have failed');
        fprintf('Not deleting temporary files so that pdf can be run again\n');
        fprintf('Ghostscript version (gs --version):\n');
        system('gs --version');
        fprintf('Try using gs version 9.02: you may need to download this\n');
      else
        % Remove temporary files
        delete(in_search_str);
        rmdir(pdf_dirs{dir_idx});
      end
    end
  end
end

fprintf('Done %s (%s)\n', param.day_seg, datestr(now));

try
  delete(map_info.fig_hand)
end
try
  delete(echo_info.fig_hand)
end

return;
