% function status = block_data(param,param_override)
% status = block_data(param,param_override)
%
% param: Parameter structure from read_param_xls parameter spreadsheet
%
% param.block_data: Structure which controls the size of each block
%  .block_size: number of columns in each block
%  .block_overlap: the percentage of overlap between each block
%  .top_gap: number of rows before the first layer
%  .bottom_pad : number of rows after the deepest layer
%  .surface_flat_en:	Enable/Disable surface flattening
%  .surface_rel_layers_flat_en:	Optional feature when surface filtering is enabled. Enable this feature to flatten the layers relative to the filtered surface.
%  .surface_filter_len:	Specifies the length of the filter for filtering the surface
%  .pre_detrend_filter_en:	Enable/Disable filtering before detrending
%  .post_detrend_filter_en:	Enable/Disable filtering after detrending (before normalization)
%  .uncompress_en:	Depending on the echogram data product used (e.g qlook, post), the echogram may be compressed. This flag when true uncompresses the compressed data using uncompress_echogram function prior to any processing.
%  .early_trunc:	Truncate data immediately after surface flattening (before detrending and normalizing)
%  .late_trunc:	Truncate data after all data manipulation( i.e detrending and normalizing ) is done.
%  .debug_plot:	Set to true for debug plots.
%  .detrend_debug:	Set to true for detrend debug plots.
%  .echo_path:	Path to echogram data, typically an argument of ct_filename_out function e.g 'CSARP\standard' => ct_filename_out(param,'CSARP\standard').
%  .out_fn:	Path where output blocks and files are saved. Currently, this is passed as an argument to ct_filename_tmp to save the outputs in KU user's scratch
%  .layers_source:	This specifies where the layer data is loaded from(e.g layerdata, records, lidar, etc). This forms a field of the layer_params struct passed into opsLoadLayers. See runOpsLoadLayers.m
%  .layerdata_source:	When layers_source is layerdata, this string specifies the layerdata (e.g layer_koenig, layer, post) to be loaded. This field is also one of the fields of the layer_params struct passed into opsLoadLayers. See runOpsLoadLayers.m
%  .regexp:	When layers_source is layerdata, all the layers with layer names that match this regular expression pattern are loaded. This field is also one of the fields of the layer_params struct passed into opsLoadLayers. See runOpsLoadLayers.m
%
% Authors: Ibikunle ( Adapted from John Paden's koenig_mat_loader )
%
% See also: run_block_data.m, block_data.m, run_unblock_data.m,
% unblock_data.m

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

[output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name);

physical_constants;

%% Input Checks: cmd
% =========================================================================

% Remove frames that do not exist from param.cmd.frms list
frames = frames_load(param);
param.cmd.frms = frames_param_cmd_frms(param,frames);

records = records_load(param);
ref = records_reference_trajectory_load(param,records);
along_track = geodetic_to_along_track(ref.lat,ref.lon,ref.elev);

%% Input Checks: block_data
% =========================================================================

% block_data.block_along_track: scalar double, along-track length of each
% block in meters, default is 5000 m
if ~isfield(param.block_data,'block_along_track') || isempty(param.block_data.block_along_track)
  param.block_data.block_along_track = 5e3; % Use default block size
end

% block_data.block_Nx: scalar integer, number of columns in each output
% block, default is 256
if ~isfield(param.block_data,'block_Nx') || isempty(param.block_data.block_Nx)
  param.block_data.block_Nx = 256; % Use default block size
end

% block_data.block_overlap: scalar double, percentage overlap for each
% block, default is 0.5, valid range is 0 to 1
if ~isfield(param.block_data,'block_overlap') || isempty(param.block_data.block_overlap)
  param.block_data.block_overlap = 0.5;
end
param.block_data.block_overlap = min(1,param.block_data.block_overlap);
param.block_data.block_overlap = max(0,param.block_data.block_overlap);

% block_data.echo_img: scalar integer, image number, default is 0 (the
% combined image Data_YYYYMMDD_SS_FFF.mat).
if ~isfield(param.block_data,'echo_img') || isempty(param.block_data.echo_img)
  param.block_data.echo_img = 0; % Use default echo_img
end

% block_data.file: structure controlling file saving and which output files
% will be generated
if ~isfield(param.block_data,'file') || isempty(param.block_data.file)
  param.block_data.file = [];
end
if ~isfield(param.block_data.file,'img_en') || isempty(param.block_data.file.img_en)
  param.block_data.file.img_en = true;
end
if ~isfield(param.block_data.file,'layer_bin_en') || isempty(param.block_data.file.layer_bin_en)
  param.block_data.file.layer_bin_en = true;
end
if ~isfield(param.block_data.file,'layer_mult_en') || isempty(param.block_data.file.layer_mult_en)
  param.block_data.file.layer_mult_en = true;
end
if ~isfield(param.block_data.file,'layer_seg_en') || isempty(param.block_data.file.layer_seg_en)
  param.block_data.file.layer_seg_en = true;
end
if ~isfield(param.block_data.file,'mat_en') || isempty(param.block_data.file.mat_en)
  param.block_data.file.mat_en = true;
end

% block_data.out_path: string specifying which output directory to put the
% block images in
if ~isfield(param.block_data,'out_path') || isempty(param.block_data.out_path)
  param.block_data.out_path = 'block_data';
end
out_fn_dir = ct_filename_out(param,param.block_data.out_path);

% block_data.rows: structure controlling truncation in row dimension
if ~isfield(param.block_data,'rows') || isempty(param.block_data.rows)
  param.block_data.rows = [];
end

% block_data.rows.t0_pad: integer controlling how many rows of the image
% will be preserved beyond the top layer. Default is inf which will
% preserve all available rows on the top of the image.
if ~isfield(param.block_data.rows,'t0_pad') || isempty(param.block_data.rows.t0_pad)
  param.block_data.rows.t0_pad = inf;
end

% block_data.rows.t1_pad: integer controlling how many rows of the image
% will be preserved below the bottom layer. Default is inf which will
% preserve all available rows on the bottom of the image.
if ~isfield(param.block_data.rows,'t1_pad') || isempty(param.block_data.rows.t1_pad)
  param.block_data.rows.t1_pad = inf;
end

%% Create blocks
% =========================================================================
dx = param.block_data.block_along_track*param.block_data.block_overlap;
X = param.block_data.block_along_track;
x0 = 0:dx:along_track(end)-X;
x1 = x0 + X;

%% Determine which blocks to create
% =========================================================================
block_mask = false(size(x0));
for frm = param.cmd.frms
  start_x = along_track(frames.frame_idxs(frm));
  if frm == length(frames.frame_idxs)
    stop_x = along_track(end);
  else
    stop_x = along_track(frames.frame_idxs(frm+1)-1);
  end
  for block_idx = 1:length(x0)
    if x0(block_idx) < stop_x && x1(block_idx) > start_x
      block_mask(block_idx) = true;
    end
  end
end

%% Load layers
% =========================================================================
[layers,layer_params] = opsLoadLayers(param,param.block_data.layer_params);

%% Block Loop
% =========================================================================
frm_mask = false(size(frames.frame_idxs));
echogram_fn_dir = ct_filename_out(param,param.block_data.echo_path);
mdata = {};
for block_idx = find(block_mask)
  
  %% Block: Frames to load
  % =========================================================================
  cat_data = [];
  for frm = 1:length(frames.frame_idxs)
    start_x = along_track(frames.frame_idxs(frm));
    if frm == length(frames.frame_idxs)
      stop_x = along_track(end);
    else
      stop_x = along_track(frames.frame_idxs(frm+1)-1);
    end
    if start_x < x1(block_idx) && stop_x > x0(block_idx)
      if ~frm_mask(frm)
        frm_mask(frm) = true;
        % Load frame into memory
        if param.block_data.echo_img == 0
          echogram_fn_name = sprintf('Data_%s_%03d.mat', param.day_seg, frm);
        else
          echogram_fn_name = sprintf('Data_img_%02d_%s_%03d.mat', param.block_data.echo_img, param.day_seg, frm);
        end
        echogram_fn = fullfile(echogram_fn_dir,echogram_fn_name);
        mdata{frm} = load_L1B(echogram_fn);
      end
      
      % Concatenate data
      cat_data = echo_concatenate(cat_data,mdata{frm});
        
    else
      if frm_mask(block_idx)
        frm_mask(frm) = false;
        % Remove frame from memory
        mdata{frm} = [];
      end
    end
  end
  
  %% Block: Condition data
  % =======================================================================

  %% Block: Extract block
  % =======================================================================
  start_rec = find(along_track >= x0(block_idx),1);
  stop_rec = find(along_track <= x1(block_idx),1,'last');
  start_idx = find(cat_data.GPS_time >= records.gps_time(start_rec),1);
  stop_idx = find(cat_data.GPS_time <= records.gps_time(stop_rec),1,'last');
  
  Nx = length(cat_data.GPS_time);
  new_axis = linspace(start_idx,stop_idx,param.block_data.block_Nx);
  cat_data.Data = interp1((1:Nx).',cat_data.Data.',new_axis.').';
  cat_data.Elevation = interp1(1:Nx,cat_data.Elevation,new_axis);
  cat_data.GPS_time = interp1(1:Nx,cat_data.GPS_time,new_axis);
  cat_data.Heading = interp1(1:Nx,cat_data.Heading,new_axis);
  cat_data.Latitude = interp1(1:Nx,cat_data.Latitude,new_axis);
  cat_data.Longitude = interp1(1:Nx,cat_data.Longitude,new_axis);
  cat_data.Roll = interp1(1:Nx,cat_data.Roll,new_axis);
  cat_data.Pitch = interp1(1:Nx,cat_data.Pitch,new_axis);
  cat_data.Surface = interp1(1:Nx,cat_data.Surface,new_axis);
  
  %% Block: Extract layers
  % =======================================================================
  Nt = length(cat_data.Time);
  Nx = length(cat_data.GPS_time);
  layers_bin_bitmap = zeros(size(cat_data.Data),'uint8');
  layers_bitmap = zeros(size(cat_data.Data),'uint16');
  layers_segment_bitmap = zeros(size(cat_data.Data),'uint16');
  layers_vector = nan(length(layers),Nx);
  for idx = 1:length(layers)
    twtt = interp1(layers(idx).gps_time,layers(idx).twtt,cat_data.GPS_time);
    twtt_bin = round(interp1(cat_data.Time, 1:length(cat_data.Time), twtt));
    good_idxs = find(isfinite(twtt_bin) & twtt_bin >= 1 & twtt_bin <= Nt);
    layers_bin_bitmap(twtt_bin(good_idxs) + Nt*(good_idxs-1)) = 1;
    layers_bitmap(twtt_bin(good_idxs) + Nt*(good_idxs-1)) = idx;
    for col = 1:Nx
      if isfinite(twtt_bin(col)) && twtt_bin(col) >= 1 && twtt_bin(col) <= Nt
        layers_segment_bitmap(twtt_bin(col):end,col) = idx;
      end
    end
    layers_vector(idx,good_idxs) = twtt_bin(good_idxs);
  end
  
  min_layer = interp_finite(min(layers_vector),0) - param.block_data.rows.t0_pad;
  max_layer = interp_finite(max(layers_vector),0) + param.block_data.rows.t1_pad;
  min_layer = min(min_layer);
  max_layer = max(max_layer);
  min_layer = min(Nt,max(1,min_layer));
  max_layer = min(Nt,max(1,max_layer));
  
  layers_bin_bitmap = layers_bin_bitmap(min_layer:max_layer,:);
  layers_bitmap = layers_bitmap(min_layer:max_layer,:);
  layers_segment_bitmap = layers_segment_bitmap(min_layer:max_layer,:);
  layers_vector = layers_vector - min_layer + 1;
  cat_data.Data = cat_data.Data(min_layer:max_layer,:);
  cat_data.Time = cat_data.Time(min_layer:max_layer);
  
  %% Block: Save
  % =======================================================================
  
  if param.block_data.file.img_en
    out_fn_name = sprintf('img_%s_%04d.png',param.day_seg,block_idx);
    out_fn = fullfile(out_fn_dir,out_fn_name);
    fprintf('%s\t%s\n', out_fn, datestr(now,'yyyymmdd_HHMMSS'));
    min_val = min(cat_data.Data(:));
    max_val = max(cat_data.Data(:));
    img_tmp = uint8((cat_data.Data-min_val)/(max_val-min_val)*255);
    imwrite(uint8(img_tmp), out_fn);
    clear img_tmp;
  end
  
  if param.block_data.file.layer_bin_en
    out_fn_name = sprintf('layer_bin_%s_%04d.png',param.day_seg,block_idx);
    out_fn = fullfile(out_fn_dir,out_fn_name);
    fprintf('%s\t%s\n', out_fn, datestr(now,'yyyymmdd_HHMMSS'));
    imwrite(uint8(layers_bin_bitmap), out_fn);
  end
  
  if param.block_data.file.layer_mult_en
    out_fn_name = sprintf('layer_mult_%s_%04d.png',param.day_seg,block_idx);
    out_fn = fullfile(out_fn_dir,out_fn_name);
    fprintf('%s\t%s\n', out_fn, datestr(now,'yyyymmdd_HHMMSS'));
    imwrite(uint8(layers_bitmap), out_fn);
  end
  
  if param.block_data.file.layer_seg_en
    out_fn_name = sprintf('layer_seg_%s_%04d.png',param.day_seg,block_idx);
    out_fn = fullfile(out_fn_dir,out_fn_name);
    fprintf('%s\t%s\n', out_fn, datestr(now,'yyyymmdd_HHMMSS'));
    imwrite(uint8(layers_segment_bitmap), out_fn);
  end
  
  if ~param.block_data.file.mat_en
    cat_data = rmfield(cat_data,'Data');
  end
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  out_fn_name = sprintf('%s_%04d.mat',param.day_seg,block_idx);
  out_fn = fullfile(out_fn_dir,out_fn_name);
  cat_data.param_block_data = param;
  cat_data.layers_vector = layers_vector;
  fprintf('%s\t%s\n', out_fn, datestr(now,'yyyymmdd_HHMMSS'));
  ct_save(out_fn,'-struct','cat_data');
  
end
status = true;