
function status = block_data(params,block_data_param)
% status = block_data(params,block_data_param)

% params: Parameter structure from read_param_xls parameter spreadsheet

% block_data_param: Structure which controls the size of each block
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


% Authors: Ibikunle ( Adapted from John Paden's koenig_mat_loader )



physical_constants;

block_data = block_data_param.block_data;

% Check conditions
if ~isfield(block_data,'block_size') || isempty(block_data.block_size)
  block_size = 256; % Use default block size
  fprintf(' Invalid/ empty block size: default value %d is used',block_size);
else
  block_size = block_data.block_size;
end

if ~isfield(block_data,'block_overlap') || isempty(block_data.block_overlap) && block_data.block_overlap > 1
  error('Invalid/ empty block overlap')
else
  block_overlap = block_data.block_overlap;
end

if ~isfield(block_data,'top_gap') || isempty(block_data.top_gap)
  error('Please specify the value of top gap used')
else
  top_gap = block_data.top_gap;
end

if ~isfield(block_data,'bottom_pad') || isempty(block_data.bottom_pad)
  error('Please specify the value of bottom gap used')
else
  bottom_pad = block_data.bottom_pad;
end

echo_path = block_data.echo_path;

surface_flat_en = block_data.surface_flat_en;

detrend_en = block_data.detrend_en;
filter_echo_en = block_data.filter_echo_en;


fn = block_data.out_fn;

layers_name = block_data.layers_name;
layers_source = block_data.layers_source;
layerdata_source = block_data.layerdata_source;
regexp_string = block_data.regexp;


% Debug Yay or nay?? :)
debug_plot = block_data.debug_plot;
detrend_debug = block_data.debug_plot;






%% Automated Section

for param_idx = 1:length(params)
  param = params(param_idx);
  day_seg_path_echo = ct_filename_out(param,echo_path);
  echo_fns_all = get_filenames(day_seg_path_echo, 'Data_2012','0','.mat',struct('recursive',true));
  start = 1; stop =length(echo_fns_all);
  echo_fns = echo_fns_all( start:stop);
  
  %% Load layer for the entire day_seg
  % Check radar type and load layers
  
  if regexp(param.radar_name,'snow')
    idx = 1;
    layer_params(idx).name = 'surface';
    layer_params(idx).source = layers_source;
    layer_params(idx).layerdata_source = layerdata_source;
    %     idx = idx + 1;
    %     layer_params(idx).name = 'bottom';
    %     layer_params(idx).source = layers_source;
    %     layer_params(idx).layerdata_source = layerdata_source;
    idx = idx + 1;
    layer_params(idx).regexp = regexp_string;
    layer_params(idx).source = layers_source;
    layer_params(idx).layerdata_source = layerdata_source; % This might need to be updated
    
    layers_struct = opsLoadLayers(param,layer_params);
    
    %     layers_struct = load('layers_cell.mat'); % remove afterwards!
    %     layers_struct = layers_cell.layers_cell;
    
    
  elseif regexp(param.radar_name,'mcords')
    
    idx = 1;
    layer_params(idx).name = 'surface';
    layer_params(idx).source = layers_source;
    layer_params(idx).layerdata_source = layerdata_source;
    idx = idx + 1;
    layer_params(idx).name = 'bottom';
    layer_params(idx).source = layers_source;
    layer_params(idx).layerdata_source = layerdata_source;
    idx = idx + 1;
    layer_params(idx).regexp = regexp_string;
    layer_params(idx).source = layers_source;
    layer_params(idx).layerdata_source = 'layer_MacGregor'; % This might need to be updated
    
    layers_struct = opsLoadLayers(param,layer_params);
    
    
  else
    warning('Radar type not recognized');
    break;
    
  end

  % Surface
  surf = [];
  surf.twtt = layers_struct(1).twtt ;
  surf.gps_time = layers_struct(1).gps_time;
  surf.elev = layers_struct(1).elev ;
  
  loaded_layers = zeros( length(layers_struct), length(surf.twtt) ) ;
  
  for iter_idx = 1 : size(loaded_layers,1)
    loaded_layers(iter_idx,:) = layers_struct(iter_idx).twtt ;
  end
  
  % Put surface and bottom in the right positions
  %   loaded_layers(1,:) = surf.twtt;
  %   loaded_layers(end+1,:) = loaded_layers(2,:) ; % move bottom to the last entry
  %   loaded_layers(2:end-1,:) = loaded_layers(3:end,:) ;
  %   loaded_layers(end,:) =[]; % remove duplicate bottom
  
  % Filter Surface
  filter_len = block_data.surface_filter_len;
  if filter_len ~=1 && surface_flat_en == 1 % Only filter surface if surface flattening is enabled
    
    surf.twtt_filtered = surf.twtt - surf.elev/(c/2);
    surf.twtt_filtered = fir_dec(surf.twtt_filtered,ones(1,filter_len)/filter_len,1);
    surf.twtt_filtered = surf.twtt_filtered + surf.elev/(c/2);
    
    if block_data.surface_rel_layers_flat_en
      % To remove high frequency variations in internal layers,
      % flatten the internal layers relative to filtered surface
      
      layers = bsxfun(@minus,loaded_layers(2:end,:),surf.twtt);
      layers = bsxfun(@plus,layers, surf.twtt_filtered);
      layers = [surf.twtt; layers];    
        
      % Layers in a cell ( This is needed for opsInterpLayersToMasterGPSTime)
      layers_cell = [];
      for layer_idx = 1: length(layers_struct)
        layers_tmp = layers_struct(layer_idx);
        layers_tmp.twtt = layers(layer_idx,:);
        layers_cell{end+1} = layers_tmp;
      end
      
    else
      layers = loaded_layers;
    end
    
    
  else
    
    surf.twtt_filtered = surf.twtt; % No filtering
    layers = loaded_layers;      
  
    % Layers in a cell ( This is needed for opsInterpLayersToMasterGPSTime)
    layers_cell = [];
    for layer_idx = 1: length(layers_struct)
      layers_cell{end+1} = layers_struct(layer_idx);
    end

    
  end
  
  % ============== End load surface and layers ===========
  
  
  % Create output folder if it does not exist
  out_dir = ct_filename_tmp(param,'',fn,sprintf('frames_%03d_%03d',start,stop));
  if ~exist(out_dir,'dir')
    mkdir(out_dir)
  end
  
  % Load frames and records data
  
  load(ct_filename_support(param,'','frames')); % load frames
  records = load(ct_filename_support(param,'','records')); % load records
  
  
  % Initializations for left-over data
  left_over = [];
  processed_day_seg = {};
  processed_years = {};
  longest_col = 0;
  frame_overlap = 0; % Flag to indicate that a block is an overlap of frames
  rangeline_offset = 0;
  
  % Iterate over each frame in the day segment
  for fn_idx = 1 : length(echo_fns)
       
    echo_fn = echo_fns{fn_idx};
    [~,fn_name] = fileparts(echo_fn);
    fprintf('%d of %d (%s)\n', fn_idx, length(echo_fns), datestr(now));
    fprintf('  %s\n', echo_fn);
    
    % Load echogram data
    tmp1 = load(echo_fn);    
    
    if 0
      %  View loaded echogram
      figure(10);
      imagesc(lp(tmp1.Data)); colormap(1-gray(256))
      title(sprintf('Loaded data %s',fn_name),'Interpreter','none')
    end
    
    if block_data.uncompress_en
      tmp = uncompress_echogram(tmp1);
    else
      tmp = tmp1;
    end
    
    [Nt,Nx] = size(tmp.Data);
    
    
    % Pre detrend along track incoherent averaging
    if block_data.pre_detrend_filter_en
      tmp.Data = echo_filt(tmp.Data,filter_len);
    end
    
    
    % Check fields in the tmp struct
    
    if ~isfield(tmp,'Roll')
      tmp.Roll = interp1(records.gps_time,records.roll,tmp.GPS_time);
    end
    
    if ~isfield(tmp,'Heading')
      tmp.Heading = interp1(records.gps_time,records.heading,tmp.GPS_time);
    end
    
    if ~isfield(tmp,'Pitch')
      tmp.Pitch = interp1(records.gps_time,records.pitch,tmp.GPS_time);
    end
    
    %% Check if roll data is in radians and convert to degrees
    
    tmp.Roll = tmp.Roll * (180/pi); % Convert roll to degrees
    
    
    % Remove rangelines that overlap with previous frame using GPS_time.
    
    if (fn_idx > 1) && (fn_idx < length(echo_fns))
      
      keep_idx = find(tmp.GPS_time>=records.gps_time(frame_idxs(fn_idx)) & tmp.GPS_time<= records.gps_time(frame_idxs(fn_idx+1)-1));
      rangeline_offset = keep_idx(1) -1;
      
      
      tmp.Data = tmp.Data(:,keep_idx);
      tmp.GPS_time =  tmp.GPS_time(keep_idx);
      
      tmp.Latitude = tmp.Latitude(keep_idx);
      tmp.Longitude = tmp.Longitude(keep_idx);
      tmp.Elevation = tmp.Elevation(keep_idx);
      tmp.Surface = tmp.Surface(keep_idx);
      tmp.Pitch = tmp.Pitch(keep_idx);
      tmp.Heading = tmp.Heading(keep_idx);
      tmp.Roll = tmp.Roll(keep_idx);
      
    end
    
    mdata = tmp; % This will be modified when there's left_over data     
   
    
    % Interpolate current frame using opsInterpLayersToMasterGPSTime    
    tmp_layers = opsInterpLayersToMasterGPSTime(tmp,layers_cell,[300 60]);    
    
    % curr_layers is the layers twtt for the current frame
    curr_layers = zeros( size(layers,1),size(tmp.Data,2) );
    
    for idx = 1: size(layers,1)
      curr_layers(idx,:) = tmp_layers.layerData{1,idx}.value{1,2}.data;
    end
    
  % Interpolate surface data for the current frame from the entire day_seg
    mdata.curr_surf_twtt = curr_layers(1,:); % filtered surface
    mdata.unfiltered_surface = interp1(surf.gps_time,surf.twtt,tmp.GPS_time); % unfiltered surface
    mdata.offset_rounding = mdata.curr_surf_twtt - mdata.unfiltered_surface; % filtering offset ( This is zero when there's no filtering )
    
    % Check surface    
    if 0
      figure;clf;
      imagesc([],tmp.Time,lp(tmp.Data)); colormap(1-gray(256))
      hold on; plot( mdata.curr_surf_twtt )
      title(sprintf('Plot of surface on Echogram data %s',fn_name),'Interpreter','none');
    end   
    
    
    if surface_flat_en == 0
      
      layer_rangebin1 = round(interp1(tmp.Time,1:length(tmp.Time),curr_layers));
      surf_index = layer_rangebin1(1,:);
      mdata.layer_rangebin = layer_rangebin1 - top_gap* ones(size(layer_rangebin1));
      
      %       Check top_gap value
      if any(mdata.layer_rangebin(1,:)  < 0 )
        error('top_gap value is too large');
        break;
      end
      
      longest_col = max( max(mdata.layer_rangebin(:)),longest_col);
      flattened_Data = tmp.Data;
      
      
    else
      % Surface flattening is enabled (e.g for snow data)
      
      surf_index = round(interp1(tmp.Time,1:length(tmp.Time),mdata.curr_surf_twtt)); % filtered surface index
      if any(~isfinite(surf_index))
        
        if all(~isfinite(surf_index))
          warning (sprintf('No surface found for %s: skipping to the next frame',fn_name));
          continue;
        else
          def_val = surf_index(find(~isnan(surf_index),1)); % should I use mean instead?
          surf_index = interp_finite(surf_index,def_val);
        end
        
      end
      
      % Convert current layers (twtt) to rangebin index
      layer_rangebin1 = zeros(size(curr_layers));
      layer_rangebin1(1,:) = surf_index; % surface index using filtered surface
      
      for layers_idx = 2:size(layers,1)
        layer_rangebin1(layers_idx,:) = round(interp1(tmp.Time,1:length(tmp.Time),curr_layers(layers_idx,:)));
      end
      
      
      %% Create a new (flattened) data matrix
      shift  = surf_index - top_gap;
      flattened_Data = [tmp.Data ; nan(Nt,Nx)];
      
      for rline = 1:Nx
        flattened_Data(:,rline) = circshift(flattened_Data(:,rline), -shift(rline));
      end
      
      % Adjust internal layers after shifting
      layer_rangebin = layer_rangebin1 - repmat(layer_rangebin1(1,:)-top_gap,size(layers,1),1) ; % change layers index to match flattened surface
      
      % I need to confirm if this wouldn't lead to any issues
      flattened_Data(isnan(flattened_Data)) = 0;
      
      % Handle situation where internal layers may go above surface
      % after shifting particularly for last(bad) frames
      layer_rangebin(layer_rangebin < min(layer_rangebin(1,:) ) ) = min(layer_rangebin(1,:) );
      mdata.layer_rangebin = layer_rangebin;
      
      stop_col = max(mdata.layer_rangebin(:)) + bottom_pad ;% Stop col using flattened surface
      longest_col = max(longest_col,stop_col);
      
      if block_data.early_trunc
        flattened_Data = flattened_Data(1:longest_col,:);
        mdata.Time = tmp.Time(1:longest_col);
      end
      
      
    end
    % =========== End of layer adjustments and surface flattening =========== %
    
    if debug_plot
      
      % Check Loaded data
      figure(100); imagesc(lp(tmp.Data));colormap(1-gray);
      a = title(sprintf(' Frame %s ',fn_name));
      ylabel(' Row index')
      set(a,'Interpreter','none' )
      for idx = 1:size(layer_rangebin1,1)
        hold on;
        plot(layer_rangebin1(idx,:))
      end
      keyboard
      
      % Check flattened_data
      figure(101); imagesc(lp(flattened_Data));colormap(1-gray);
      a = title(sprintf('Flattened Frame %s ',fn_name));
      ylabel(' Row index')
      set(a,'Interpreter','none' )
      for idx = 1:size(layer_rangebin,1)
        hold on;
        plot(mdata.layer_rangebin(idx,:))
      end
      
    end
    
    % Is Pre_detrend filter enabled?
    if filter_echo_en % If filtering is enabled
      flattened_Data = echo_filt(flattened_Data,11);
    end
    
    frame_rline = (1:length(surf_index)) + rangeline_offset; % range line index
    mdata.surf_index = surf_index;
    
    
    %% Create layer raster from layer_rangebin
    Layer = zeros(size(flattened_Data));
    for layer_idx = 1:size(layer_rangebin,1)
      for layer_col = 1:size(layer_rangebin,2)
        temp = round(layer_rangebin(layer_idx,layer_col)); % This idx represents the column of the layer
        %         temp = temp + top_gap; % Add top_gap offset
        if ~isnan(temp)
          if Layer(temp,layer_col) == 0
            Layer(temp,layer_col) = layer_idx;
          end
        else
          continue;
        end
      end
    end
    
    % ================  End: Create layer  =========================%
    
    
    
    %%     Merge data and meta-data when there's left-over from previous frame
    
    if ~isempty(left_over) && fn_idx ~= length(echo_fns)
      threshold = (size(flattened_Data,1) - size(left_over.data,1)); % Check which has more rows ("range bins")
      
      if threshold < 0 % left_over has more, zero pad flattened_Data
        flattened_Data = [flattened_Data ; zeros(abs(threshold),size(flattened_Data,2)) ];
        Layer = [Layer ; zeros(abs(threshold),size(flattened_Data,2)) ];
      end
      
      appended_data = zeros(size(flattened_Data,1),size(left_over.data,2));
      appended_data(1:size(left_over.data,1),:) = left_over.data;
      appended_layer = zeros(size(Layer,1),size(left_over.layer,2));
      appended_layer(1:size(left_over.layer,1),:) = left_over.layer;
      
      flattened_Data = [appended_data  flattened_Data ];
      Layer = [appended_layer Layer];
      
      % Modified meta data
      mdata.GPS_time =  [left_over.GPS_time_left tmp.GPS_time ];
      mdata.Latitude = [left_over.Latitude_left mdata.Latitude];
      mdata.Longitude = [left_over.Longitude_left mdata.Longitude];
      mdata.Elevation = [left_over.Elevation_left mdata.Elevation];
      mdata.Roll = [left_over.Roll_left mdata.Roll];
      mdata.Pitch = [left_over.Pitch_left mdata.Pitch];
      mdata.Heading = [left_over.Heading_left mdata.Heading];
      
      mdata.surf_index = [left_over.surf_index_left mdata.surf_index];
      mdata.curr_surf_twtt = [left_over.Surface_left mdata.curr_surf_twtt];
      
      mdata.unfiltered_surface = [left_over.unfiltered_surface_left mdata.unfiltered_surface];
      
      mdata.offset_rounding = [left_over.offset_rounding_left mdata.offset_rounding];
      mdata.layer_rangebin = [left_over.layer_rangebin_left mdata.layer_rangebin];
      frame_rline = [left_over.frame_rline_left frame_rline];
      frame_overlap = 1;
      
    end
    
    
    
    %% Search for valid data and break data into blocks
    good_rline = ~any(isnan(flattened_Data));
    rline = find(good_rline,1);
    
    while ~isempty(rline)
      rline_end = find(~good_rline(rline:end),1)-2 + rline;
      if isempty(rline_end)
        rline_end = size(flattened_Data,2);
      end
      if rline_end >= rline + block_size
        rline_end = rline + block_size-1;
      end
      
      if rline_end-rline+1 >= block_size || fn_idx == length(echo_fns)
        
        %% Create directories to save data
        
        day_seg = fn_name(6:16);
        
        if isempty(processed_day_seg) || ~ismember(day_seg,processed_day_seg)
          processed_day_seg{end+1} = day_seg;
        end
        
        year = fn_name(6:9);        
        if isempty(processed_years) || ~ismember(year,processed_years)
          processed_years{end+1} = year;
          
          block = 1;
          
          % Create output directories for the day_seg
          base_layer_bin = fullfile(out_dir,'layer_bin');
          base_layer = fullfile(out_dir,'layer');
          base_image = fullfile(out_dir,'image');
          base_figures = fullfile(out_dir,'figures');
          
          
          if ~exist(base_layer_bin, 'dir')
            mkdir(base_layer_bin);
          end
          
          if ~exist(base_layer, 'dir')
            mkdir(base_layer);
          end
          
          if ~exist(base_image, 'dir')
            mkdir(base_image);
          end
          
          if ~exist(base_figures, 'dir')
            mkdir(base_figures);
          end
          
        else
          block= block+1;
          
        end
        
        %% Break (data and layer) into Nt by block_size images (slow time)
        
        data = flattened_Data(:,rline:rline_end);
        layer = Layer(:,rline:rline_end);
        raster = uint8(layer);
        
        if detrend_debug
          figure(400);clf;
          imagesc(lp(data)); colormap(1-gray(256));
          title(sprintf('Data before detrending:  block%d, %s',block,fn_name),'Interpreter','None' )
          
          figure(500);clf;
          title(sprintf('A-scope before detrending:  block%d, %s',block,fn_name),'Interpreter','None' )
          for idx = 1:size(data,2)
            hold on;
            plot(data(:,idx));
          end
        end
        
        
        % For snow data; bottom should be set to nan so echo_detrend uses data
        % from surface to Nt
        if regexp(param.radar_name,'snow')
          dt_bottom = nan(1,size(data,2));
        else
          dt_bottom = interp1(1:length(mdata.Time),mdata.Time,mdata.layer_rangebin(end,rline:rline_end));
        end
        
        if detrend_en
          detrend_param = block_data.norm_detrend_params;
          detrend_param.layer_top = interp1(1:length(mdata.Time),mdata.Time,mdata.layer_rangebin(1,rline:rline_end));
          detrend_param.layer_bottom= dt_bottom;
          detrend_param.roll= mdata.Roll(rline:rline_end);
          detrend_param.block = block;
          detrend_param.fn_name = fn_name;
          detrend_param.out_dir = out_dir;
          
          % Let's see...what's the min and max roll for this block?
          sprintf('Max %2.3f, Min %2.3f', max(detrend_param.roll), min(detrend_param.roll))
          
          % It's time to detrend!!! Detrend to correct large dynamic range in data
          dt_Data = echo_detrend(struct('Time',mdata.Time,'Data',lp(data)),detrend_param);
          dt_mask = isfinite(dt_Data);
          
          if detrend_debug
            figure(401);clf;
            imagesc(dt_Data); colormap(1-gray(256));
            title(sprintf('Data after detrending:  block%d, %s',block,fn_name),'Interpreter','None' )
            
            figure(501);clf;
            title(sprintf('A-scope after detrending:  %s block%d,',fn_name,block),'Interpreter','None' )
            
            for idx = 1:size(dt_Data,2)
              hold on;
              plot(dt_Data(:,idx))
            end
          end
          
        else
          
          dt_Data = 10*log10(data); % Convert to log; no detrending
          
        end
        
        
        % Post detrend along track filtering: reduce 'em noise...
        if block_data.post_detrend_filter_en
          filter_len = block_data.norm_detrend_params.filter_len;
          
          % Convert -inf to nan for nan_fir_dec
          dt_Data(isinf(dt_Data)) = nan;
          filt_Data = 10*log10( nan_fir_dec( 10.^(dt_Data/10),ones(1,filter_len)/filter_len,1) );
          
        else
          filt_Data = dt_Data;
        end
        
        if detrend_debug
          figure(402);clf;
          imagesc(filt_Data); colormap(1-gray(256));
          title(sprintf('Data after detrending and filtering:  block%d, %s',block,fn_name),'Interpreter','None' )
          
          figure(502);clf;
          title(sprintf('A-scope after detrending and filtering: log scale %s block%d,',fn_name,block),'Interpreter','None' )
          
          for idx = 1:size(filt_Data,2)
            hold on;
            plot(filt_Data(:,idx))
          end
          
        end
        
        
        scale_min = block_data.norm_detrend_params.scale_min;
        scale_max = block_data.norm_detrend_params.scale_max;
        
        % "window" specifies that rangebins that'd be used to estimate the noise floor
        
        if any(~isnan(dt_bottom))
          % Estimate noise floor
          %           nf = nanmean(echo_noise(struct('Time',mdata.Time,'Data',mfilt_Data),struct('window',[dt_bottom+5e-6; inf(size(dt_bottom))])));
          % Adaptively set scale_min and scale_max ??
          
          norm_Data = echo_norm(struct('Time',mdata.Time,'Data',filt_Data),struct('scale',[scale_min scale_max],'window',[dt_bottom+5e-6; inf(size(dt_bottom))]));
        else
          valid_max_range_dB = [0 inf];
          norm_window = block_data.norm_detrend_params.norm_window;
          
          % Estimate noise floor
          nf = echo_noise(struct('Time',mdata.Time,'Data',filt_Data),struct(...
            'window_units','%','window',norm_window,'valid_max_range_dB',valid_max_range_dB));
          
          nf = nanmean(nf(isfinite(nf)));
          
          surf_idx = min(mdata.layer_rangebin(1,:));
          sc_Data = filt_Data(surf_idx:end,:); % scaling data: data used to determine min/max scale
          sc_Data = sc_Data(dt_mask(surf_idx:end,:));
          sc_min = nanmin(sc_Data(:));
          sc_max = nanmax(sc_Data(:));
          
          scale_min = ( nf - sc_min  ) / ( sc_max - sc_min );
          scale_max = 1 - scale_min;
          
          norm_Data = echo_norm(struct('Time',mdata.Time,'Data',filt_Data),struct('scale',[scale_min scale_max],'window_units','%','window',norm_window,'valid_max_range_dB',valid_max_range_dB));
        end
        
        if block_data.late_trunc
          norm_Data = norm_Data ( top_gap : longest_col+bottom_pad,:); % Truncate after the longest data
        end
        
        if detrend_debug % Again :(
          figure(403);clf;
          imagesc(norm_Data);
          title(sprintf('Image after normalizing block%d, %s',block,fn_name),'Interpreter','None' )
          colormap(gray(256));
          caxis([0 1]);
          link_figures([400,401,403]);
          
          figure(503);clf;
          title(sprintf('A-scope after normalizing: log scale block%d, %s',block,fn_name),'Interpreter','None' )
          for idx = 1:size(norm_Data,2)
            hold on;
            plot(norm_Data(:,idx))
          end
        end
        
        % Check if there's a row of zeros and truncate
        % This may occur as a result of zero-padding for left over case
        if  any ( all(~isfinite(data),2) )
          zero_index = find ( all(~isfinite(data),2),1,'first');
          data = data(1:zero_index-1,:) ;
          layer = layer(1:zero_index-1,:) ;
        end
        
        
        data2 = single(norm_Data); % Just to make sure output data is "single"
        
        
        %% Mirror last block with insufficient rline
        if rline_end-rline+1 < block_size && fn_idx == length(echo_fns)
          % Mirror last data to complete 256 cols
          
          extra = block_size - size(data,2);
          data2 = [data2 data2(:,end-extra+1:end) ];
          layer = [layer layer(:,end-extra+1:end) ];
          raster = [raster raster(:,end-extra+1:end) ];
          
        end
        
        if debug_plot % Check flattened_data
          figure(102);
          
          subplot(211)
          imagesc(lp(flattened_Data));colormap(1-gray);
          a = title(sprintf(' Flattened Frame %s ',fn_name));
          ylabel(' Row index')
          set(a,'Interpreter','none' )
          for idx = 1:size(mdata.layer_rangebin,1)
            hold on;
            plot(layer_rangebin(idx,:))
          end
          
          subplot(212)
          imagesc((data));colormap(1-gray);
          
          a = title(sprintf(' Flattened block Frame %s_block%d ',fn_name,block));
          ylabel(' Row index')
          set(a,'Interpreter','none' )
          
          for idx = 1:size(layer_rangebin,1)
            hold on;
            plot(mdata.layer_rangebin(idx,rline:rline_end))
          end
          
          keyboard
        end
        
        
        %% Save meta-data
        
        tmp_block_data = [];
        tmp_block_data.block = block;
        
        tmp_block_data.data = data2;
        tmp_block_data.surf_index = mdata.surf_index(rline:rline_end);
        tmp_block_data.GPS_time = mdata.GPS_time (rline:rline_end);
        tmp_block_data.Latitude = mdata.Latitude (rline:rline_end);
        tmp_block_data.Longitude = mdata.Longitude (rline:rline_end);
        tmp_block_data.Elevation = mdata.Elevation (rline:rline_end);
        tmp_block_data.Roll = mdata.Roll (rline:rline_end);
        tmp_block_data.Pitch = mdata.Pitch (rline:rline_end);
        tmp_block_data.Heading = mdata.Heading (rline:rline_end);
        tmp_block_data.frame_overlap_flag = frame_overlap;
        
        tmp_block_data.norm_detrend_params = block_data.norm_detrend_params;
        
        tmp_block_data.fn_name = fn_name;
        tmp_block_data.day_seg = day_seg;
        tmp_block_data.filt_surface  = mdata.curr_surf_twtt(rline:rline_end); % filtered surface
        tmp_block_data.unfiltered_surface = mdata.unfiltered_surface(rline:rline_end);
        
        tmp_block_data.raster = raster;
        tmp_block_data.layer = mdata.layer_rangebin(:,rline:rline_end);
        tmp_block_data.rline = frame_rline(rline);
        tmp_block_data.rline_end = frame_rline(rline_end);
        
        tmp_block_data.offset_rounding = mdata.offset_rounding(:,rline:rline_end);
        
        tmp_block_data.block_size = block_data.block_size ;
        tmp_block_data.block_overlap = block_data.block_overlap ;
        tmp_block_data.top_gap = block_data.top_gap ;
        tmp_block_data.bottom_pad = block_data.bottom_pad;
        
        
        out_fn = fullfile(out_dir,'image',sprintf('image_%06d.mat',block));
        fprintf('    Save %s\n', out_fn);
        save(out_fn,'-struct','tmp_block_data');
        
        
        
        % Save image as .tiff
        
        out_fn = fullfile(out_dir,'image',sprintf('image_%06d.tiff',block));
        fprintf('    Save %s\n', out_fn);
        
        t = Tiff(out_fn, 'w');
        tagstruct.ImageLength = size(data2, 1);
        tagstruct.ImageWidth = size(data2, 2);
        tagstruct.Compression = Tiff.Compression.None;
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
        tagstruct.BitsPerSample = 32;
        tagstruct.SamplesPerPixel = size(data2,3);
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        t.setTag(tagstruct);
        t.write(data2);
        t.close();
        
        
        %% Save images as .png
        
        out_fn = fullfile(out_dir,'layer',sprintf('layer_%06d.png',block));
        fprintf('    Save %s\n', out_fn);
        imwrite(raster,out_fn);
        
        out_fn = fullfile(out_dir,'layer_bin',sprintf('layer_binary_%06d.png',block));
        fprintf('    Save %s\n', out_fn);
        imwrite(logical(raster),out_fn);
        
        out_fn = fullfile(out_dir,'figures',sprintf('image_fig_%06d.png',block));
        
        
        %% Create and save echogram figure with layers plotted
        
        hdl2 = figure(1000);
        curr_plot = imagesc(data2);colormap(1-gray);
        
        a = title(sprintf(' Flattened Frame %s_block%d ',fn_name,block));
        ylabel(' Row index')
        set(a,'Interpreter','none' )
        
        for idx = 1:size(mdata.layer_rangebin,1)
          hold on;
          plot(mdata.layer_rangebin(idx,rline:rline_end));
        end
        
        saveas(curr_plot,out_fn)
        close(hdl2)
        
        %          ============  End of save image of echo and layers   ==========
        
        
        % Re-initialize frame_overlap and rline for next loop
        frame_overlap = 0; % turn off frame overlap flag
        rline = find(good_rline(rline_end+1:end),1)+rline_end - block_overlap*block_size;
        
        left_over = []; % re-initialize
        
      else
        
        left_over.data = flattened_Data(:,rline:end);
        
        left_over.layer = Layer(:,rline:end);
        
        left_over.GPS_time_left =  mdata.GPS_time(:,rline:end);
        left_over.Latitude_left =  mdata.Latitude(:,rline:end);
        left_over.Longitude_left = mdata.Longitude(:,rline:end);
        left_over.Elevation_left = mdata.Elevation(:,rline:end);
        left_over.Roll_left =     mdata.Roll(:,rline:end);
        left_over.Pitch_left =    mdata.Pitch(:,rline:end);
        left_over.Heading_left = mdata.Heading(:,rline:end);
        
        left_over.surf_index_left = mdata.surf_index(:,rline:end);
        left_over.Surface_left   = mdata.curr_surf_twtt (:,rline:end);
        left_over.offset_rounding_left = mdata.offset_rounding(:,rline:end);
        left_over.unfiltered_surface_left = mdata.unfiltered_surface(:,rline:end);
        
        left_over.layer_rangebin_left = mdata.layer_rangebin(:, rline:end);
        left_over.frame_rline_left = frame_rline(rline:end);
        
        rline = find(good_rline(rline_end+1:end),1)+rline_end; % This forces rline to be empty
        break
        
      end
      
    end % End while
    
    
  end
  status = true ;
end