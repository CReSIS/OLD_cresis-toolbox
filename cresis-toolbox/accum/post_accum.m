function post_accum(param, param_override)
% post_accum(param_fn, day_seg_filter)
%
% Create posting function called from master (for accumulation radar
% and FMCW radars).
%
% Author: John Paden
%
% See also: master

% =====================================================================
% General Setup
% =====================================================================
% clear; % Useful when running as script
% close all; % Optional
physical_constants;
post_accum_tstart = tic;
fprintf('\n\n==============================================\n\n');

% =====================================================================
% User Settings
% =====================================================================
%param = []; % Uncomment if running as a script
if ~exist('param','var') || isempty(param)
  param = read_param_xls('/cresis/projects/dev/csarp_support/params/snow_param_2011_Greenland_P3.xls','20110325_11');
  
  clear('param_override');
  %param_override.sched.type = 'no scheduler';
  
  % Input checking
  if ~exist('param','var')
    error('A struct array of parameters must be passed in\n');
  end
  global gRadar;
  if exist('param_override','var')
    param_override = merge_structs(gRadar,param_override);
  else
    param_override = gRadar;
  end
  
elseif ~isstruct(param)
  % Functional form
  param();
end
param = merge_structs(param, param_override);

% =======================================================================
% Loop through and process each file, store processed output
%  (one output file per input file, one echogram jpg per input file, one
%   map jpg per input file, create pdfs and kml browse files)
%   auto-track surface for optimal echogram jpg
% =======================================================================
fprintf('\n=============================================================\n');
fprintf('Accum/FMCW post processing for %s (%.1f sec)\n', param.day_seg, toc(post_accum_tstart));
fprintf('=============================================================\n');

use_vectors_flag = true;
if use_vectors_flag
  % Use vectors file to create list of files that can be posted
  
  vectors_fn = ct_filename_support(param,param.post.map.vectors_fn,'vectors');
  load(vectors_fn);
  
  fprintf('Getting files for %s (%.1f sec)\n', vectors_fn, toc(post_accum_tstart));
  file_idxs = 1:length(vectors.file);
  fprintf('  Found %d files in %s\n', length(file_idxs), vectors_fn);
  if isempty(param.cmd.frms)
    for file_num = file_idxs
      param.cmd.frms(end+1) = vectors.file(file_num).idx;
    end
  end
else
  % Use data files to create list of files that can be posted
  
  base_dir = fullfile(param.vectors.file.base_dir,param.vectors.file.adc_folder_name);

  fprintf('Getting files for %s (%.1f sec)\n', base_dir, toc(post_accum_tstart));
  fns = get_filenames(base_dir,param.vectors.file.file_prefix,'','.dat');
  fprintf('  Found %d files in %s\n', length(fns), base_dir);

  if param.qlook.file.stop_idx > length(fns)
    stop_idx = length(fns);
  else
    stop_idx = param.qlook.file.stop_idx;
  end
  file_idxs = param.qlook.file.start_idx:stop_idx;

  if isempty(param.cmd.frms)
    % If frames to process is empty, process all frames
    for file_num = 1:length(file_idxs)
      fn = fns{file_idxs(file_num)};
      finfo = fname_info_accum(fn);
      param.cmd.frms(end+1) = finfo.file_idx;
    end  
  end
end

fprintf('Loading map image template (%.1f sec)\n', toc(post_accum_tstart));

if ~strcmpi(param.post.map.type,'none') && param.post.map.en
  % =========================================================================
  % Create map template which will be used for plotting the map of
  % each segment
  % =========================================================================
  fprintf('Creating map template (%.1f sec)\n', toc(post_accum_tstart));
  
  map_param.type = param.post.map.type;
  map_param.location = param.post.map.location;
  map_param.fig_hand = 1;
  map_param.map_title = '';
  map_param.decimate_seg = false;
  map_param.resample = false;
  map_param.year = str2double(param.day_seg(1:4));
  map_param.month = str2double(param.day_seg(5:6));
  map_param.day = str2double(param.day_seg(7:8));

  map_info = publish_map('setup',map_param);

  % =======================================================================
  % For each day segment, find the [x,y] extent by looking at all the frames
  % in this segment
  % =======================================================================
  vectors_fn = ct_filename_support(param,param.post.map.vectors_fn,'vectors');
  load(vectors_fn);
  
  [x,y] = projfwd(map_info.proj,vectors.lat,vectors.lon);
  x = x/1e3;
  y = y/1e3;
  min_x = min(x);
  max_x = max(x);
  min_y = min(y);
  max_y = max(y);
  frm_idx = 1;
  day_segs_idx(frm_idx) = 1;
  day_seg_x{1} = x;
  day_seg_y{1} = y;
end

% =====================================================================
% Determine which frames will actually be posted
% =====================================================================

% Get small records file info for start/stop idxs
frames_fn = ct_filename_support(param,'','frames');
if exist(frames_fn, 'file')
  frames = load(frames_fn);
else
  frames.proc_mode = [];
end

file_nums = 1:length(file_idxs);
file_num_mask = logical(zeros(size(file_nums)));
for file_num = file_nums
  if use_vectors_flag
    frm = vectors.file(file_num).idx;
  else
    fn = fns{file_idxs(file_num)};
    finfo = fname_info_accum(fn);
    frm = finfo.file_idx;
  end
  
  if ~isempty(find(param.cmd.frms == frm, 1)) ...
    && (isempty(frames.proc_mode) || mod(floor(frames.proc_mode(frm+1)/10),10) == 0)
    file_num_mask(file_num) = 1;
  end
end
file_nums = file_nums(file_num_mask);

% =====================================================================
% Post frames
% =====================================================================
vectors_idx = 0;
% Initialize list of numeric frames IDs that will be concatenated into 
% a single posted echogram to []
frm_concat_list = [];
if strcmp(param.post.img_type,'jpg')
  print_device = '-djpeg';
elseif strcmp(param.post.img_type,'png')
  print_device = '-dpng';
else
  error('Unsupported image type %s', param.post.img_type);
end
print_dpi = sprintf('-r%d', param.post.img_dpi);

for file_num_idx = 1:length(file_nums)
  file_num = file_nums(file_num_idx);
  if use_vectors_flag
    frm = vectors.file(file_num).idx;
  else
    fn = fns{file_idxs(file_num)};
    finfo = fname_info_accum(fn);
    frm = finfo.file_idx;
  end
  
  % =====================================================================
  % Read in frame
  % =====================================================================
  frm_id = sprintf('%s_%03.0f', param.day_seg, frm);
  in_name = sprintf('Data_%s.mat', frm_id);
  in_dir = ct_filename_out(param, ...
    param.post.in_dir, 'CSARP_qlook');
  in_fn = fullfile(in_dir,in_name);
   
  fprintf('  Loading file %s (%.1f sec)\n', in_fn, toc(post_accum_tstart));
  if length(frm_concat_list) == 0
    % If this is the first frame to be loaded to the list:
    mdata = load(in_fn);
  else
    % If there are already frames that have been loaded in the list, we
    % concatenate:
    tmp = load(in_fn);
    mdata.Data = cat(2,mdata.Data,tmp.Data);
    mdata.Surface = cat(2,mdata.Surface,tmp.Surface);
    mdata.Latitude = cat(2,mdata.Latitude,tmp.Latitude);
    mdata.Longitude = cat(2,mdata.Longitude,tmp.Longitude);
    mdata.Elevation = cat(2,mdata.Elevation,tmp.Elevation);
    mdata.GPS_time = cat(2,mdata.GPS_time,tmp.GPS_time);
  end

  % =====================================================================
  % Check to see if we should post the data we have accumulated
  %   1. Post if we have accumulated the number that the user has asked
  %      us to combine together
  %   2. Post if we are on the last file to be posted
  %   3. Post if there is a data gap before the next good file
  % =====================================================================
  frm_concat_list(end+1) = frm;
  if length(frm_concat_list) >= param.post.num_frm_combine ...
      || file_num == file_nums(end) ...
      || (file_num_idx < length(file_nums) && file_nums(file_num_idx+1) > file_num + 1)
    % We have accumulated enough files, create posting output

    % Create frame id string
    if length(frm_concat_list) > 1
      frm_id = sprintf('%s_%03.0f_%03.0f', param.day_seg, frm_concat_list(1), frm_concat_list(end));
    else
      frm_id = sprintf('%s_%03.0f', param.day_seg, frm_concat_list);
    end

    % Create output directory
    image_dir = fullfile(ct_filename_out(param, ...
      param.post.out_dir, 'CSARP_post', true),'images',param.day_seg);
    if ~exist(image_dir,'dir')
      mkdir(image_dir)
    end

    % Create time stamp to be used in output filenames
    time_stamp_str = datestr(epoch_to_datenum(mdata.GPS_time(1)),'HHMMSS');

    if ~strcmpi(param.post.map.type,'none') && param.post.map.en
      % ===================================================================
      % Plot Flightlines over Geotiff
      % ===================================================================
      lay.Latitude = [];
      lay.Longitude = [];
      for file_idx = 1:length(vectors.file)
        if ~isempty(find(vectors.file(file_idx).idx == frm_concat_list))
          lay.Latitude(end+1) = vectors.lat(file_idx);
          lay.Longitude(end+1) = vectors.lon(file_idx);
        end
      end

      % Plot the map for each frame
      map_param.day_seg_x = day_seg_x{1}*1000;
      map_param.day_seg_y = day_seg_y{1}*1000;
      [frame_X,frame_Y] = projfwd(map_info.proj,mdata.Latitude,mdata.Longitude);
      map_param.frame_X = {frame_X};
      map_param.frame_Y = {frame_Y};
      map_info = publish_map('delete',map_param,map_info);
      map_info = publish_map('plot',map_param,map_info);

      % Remove current file
      map_fn = sprintf('%s*_0maps.*',frm_id);
      map_fn = fullfile(image_dir,map_fn);
      delete(map_fn);

      % Save file
      map_fn = sprintf('%s_%s_0maps.%s',frm_id,time_stamp_str,param.post.img_type);
      map_fn = fullfile(image_dir,map_fn);
      fprintf('    Saving output %s\n', map_fn);
      print(map_param.fig_hand,print_device,print_dpi,map_fn);

    end

    % ===================================================================
    % Create echogram plot
    % ===================================================================
    echo_param.fig_hand = 2;
    echo_param.num_x_tics = 6;
    echo_param.frm_id = frm_id;
    echo_param.depth = param.post.depth_rng;
    echo_param.elev_comp = true;
    echo_param.er_ice = param.post.er_ice;

    lay.GPS_time = mdata.GPS_time;
    lay.Elevation = mdata.Elevation;
    lay.layerData{1}.value{1}.data = NaN*zeros(size(mdata.Surface));
    lay.layerData{1}.value{2}.data = mdata.Surface;
    lay.layerData{2}.value{1}.data = NaN*zeros(size(mdata.Surface));
    lay.layerData{2}.value{2}.data = NaN*zeros(size(mdata.Surface));
    echo_info = publish_echogram(echo_param,mdata,lay);
    set(echo_info.h_surf,'Visible','off');
    set(echo_info.h_bot,'Visible','off');
    if length(frm_concat_list) > 1
      title(sprintf('"%s" %s Frame IDs: %s', param.radar_name, param.post.mission_name, frm_id),'Interpreter','none');
    else
      title(sprintf('"%s" %s Frame ID: %s', param.radar_name, param.post.mission_name, frm_id),'Interpreter','none');
    end
    
    % Remove current file
    echo_fn = sprintf('%s*_1echo.%s',frm_id,param.post.img_type);
    echo_fn = fullfile(image_dir,echo_fn);
    delete(echo_fn);
    
    echo_fn = sprintf('%s_%s_1echo.%s',frm_id,time_stamp_str,param.post.img_type);
    echo_fn = fullfile(image_dir,echo_fn);
    set(echo_param.fig_hand,'PaperOrientation','Portrait');
    set(echo_param.fig_hand,param.post.plot_params{1},param.post.plot_params{2});
    fprintf('    Saving output %s\n', echo_fn);
    print(echo_param.fig_hand,print_device,print_dpi,echo_fn);
    
    % ===================================================================
    % Reset variables to start concatenating a new output
    % ===================================================================
    frm_concat_list = [];
  end
end

return;



