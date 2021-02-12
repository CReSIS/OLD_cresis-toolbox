% function combine_passes(param,param_override)
% combine_passes(param,param_override)
%
% Combines multiple passes (SAR focussed echograms from array processing or
% SAR focused single look complex images) into a single file for processing
% by multipass.multipass.
%
% param = struct with processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Example:
%  See multipass.run_combine_passes.m for examples of how to run this function.
%
% Authors: Cody Barnett, Bailey Miller, John Paden
%
% See also: multipass.combine_passes.m, multipass.run_combine_passes.m,
% multipass.multipass.m, multipass.run_multipass.m

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.combine_passes.pass_name, datestr(now));
fprintf('=====================================================================\n');

%% Input checking

physical_constants;

% dist_min: scalar containing the minimum distance in meters from the
% start/stop points that a trajectory must be to consider it to have
% entered or left the start/stop location. Default is 300 m.
if ~isfield(param.combine_passes,'dist_min') || isempty(param.combine_passes.dist_min)
  param.combine_passes.dist_min = 300;
end
dist_min = param.combine_passes.dist_min;

% master_pass_idx: index into passes array of which pass the param
% structure will be pulled from. This will control which output season the
% data are stored in. Default is 1 (the first pass).
if ~isfield(param.combine_passes,'master_pass_idx') || isempty(param.combine_passes.master_pass_idx)
  param.combine_passes.master_pass_idx = 1;
end

% out_path: string containing out_path to be passed to ct_filename_out
if ~isfield(param.combine_passes,'out_path') || isempty(param.combine_passes.out_path)
  param.combine_passes.out_path = 'multipass';
end

% start: structure with lat/lon field containing the start location of the
% passes (start.lat and start.lon)
start = [];
start.lat = param.combine_passes.start.lat;
start.lon = param.combine_passes.start.lon;

% stop: same as start only for the stop location. Technically the
% start/stop ordering is arbitrary since the direction of a pass does not
% matter. Therefore a pass can enter at the stop location and leave at the
% start location and vice versa.
stop = [];
stop.lat = param.combine_passes.stop.lat;
stop.lon = param.combine_passes.stop.lon;

% pass_name: string containing the pass name which is used to generate the
% filename in CSARP_multipass. Must be specified.
pass_name = param.combine_passes.pass_name;

% input_type: string containing "echo" for echogram combining and "sar" for
% SAR (coherent single look complex) image combining
if ~isfield(param.combine_passes,'input_type') || isempty(param.combine_passes.input_type)
  param.combine_passes.input_type = 'echo';
end

% surf_layer: surface layer structure that is used to combine images (only
% used for "sar" inputs)
if ~isfield(param.combine_passes,'surf_layer') || isempty(param.combine_passes.surf_layer)
  param.combine_passes.surf_layer = struct();
end
if ~isfield(param.combine_passes.surf_layer(1),'name') || isempty(param.combine_passes.surf_layer(1).name)
  param.combine_passes.surf_layer(1).name = 'surface';
end
if ~isfield(param.combine_passes.surf_layer(1),'source') || isempty(param.combine_passes.surf_layer(1).source)
  param.combine_passes.surf_layer(1).source = 'layerData';
end

% passes: struct array containing information about each pass to be
% combined. Must be specified.
passes = param.combine_passes.passes;
% Check input arguments of each pass
for passes_idx = 1:length(passes)
  % day_seg: string containing segent ID string YYYYMMDD_SS
  if ~isfield(passes(passes_idx),'day_seg') || isempty(passes(passes_idx).day_seg)
    error('passes(%d).day_seg must be defined.', passes_idx);
  end
  
  % frms: interger list of frames to load
  if ~isfield(passes(passes_idx),'frms') || isempty(passes(passes_idx).frms)
    error('passes(%d).frms must be defined.', passes_idx);
  end
  
  % param_fn: string to parameter spreadsheet.
  if ~isfield(passes(passes_idx),'param_fn') || isempty(passes(passes_idx).param_fn)
    error('passes(%d).param_fn must be defined.', passes_idx);
  end
  
  if ~isfield(passes(passes_idx),'in_path') || isempty(passes(passes_idx).in_path)
    if strcmp(param.combine_passes.input_type,'echo')
      % ECHOGRAM: Default is CSARP_standard
      passes(passes_idx).in_path = 'standard';
    else
      % SAR: Default is determined by load_sar_data.m
      passes(passes_idx).in_path = '';
    end
  end
  
  % wf_adc: Contains a cell array of images. Each image is a wf-adc pair list to load. Each image may
  % contain only one wf-adc pair. Therefore each entry in the cell array
  % must be 1x2. This field is only required for input_type == "sar".
  if strcmp(param.combine_passes.input_type,'sar')
    % Ensure that each wf_adc list is a cell array
    if ~iscell(passes(passes_idx).imgs)
      passes(passes_idx).imgs = {passes(passes_idx).imgs};
    end
    for img = 1:length(passes(passes_idx).imgs)
      if size(passes(passes_idx).imgs{img}, 1) ~= 1
        error('Each entry in the passes(%d).wf_adc list must contain only a single row/wf-adc pair, but detected %d rows for image %d.', ...
          passes_idx, size(passes(passes_idx).imgs{img},1), img);
      end
    end
  end
end
% Store any updates back into the parameter structure since this will be
% saved in the output file.
param.combine_passes.passes = passes;

%% Load echogram data for each pass
% =========================================================================
metadata = {};
data = {};
for passes_idx = 1: length(passes)
  %% Load: Load param xls
  param_fn = ct_filename_param(passes(passes_idx).param_fn);
  % Only load the parameter spreadsheet if it is not already in memory
  if passes_idx == 1 || ~strcmp(passes(passes_idx).day_seg,passes(passes_idx-1).day_seg)
    param_pass = read_param_xls(param_fn,passes(passes_idx).day_seg);
    param_pass = merge_structs(param_pass, param_override);
    
    if strcmp(param.combine_passes.input_type,'sar')
      %% Load: Load surface layer
      % Only required to combine SAR images together with img_combine
      surf_layer = opsLoadLayers(param_pass,param.combine_passes.surf_layer);
    end
  end
  
  % Use the master pass parameters to update the main "param" structure.
  % The master pass determines which season the output files will be stored
  % in.
  if passes_idx == param.combine_passes.master_pass_idx
    param = merge_structs(param_pass,param);
  end
  
  %% Load: Load Data
  if strcmp(param.combine_passes.input_type,'echo')
    %% Load: Load "echo" Data
    echo_fn_dir{passes_idx} = ct_filename_out(param_pass, passes(passes_idx).in_path); %creates directory for pass, from given data format
    % Loop to load each echogram frame
    for frm_idx = 1:length(passes(passes_idx).frms)
      echo_fn_name{passes_idx} = sprintf('Data_%s_%03d.mat',passes(passes_idx).day_seg,passes(passes_idx).frms(frm_idx));
      echo_fn{passes_idx} = fullfile(echo_fn_dir{passes_idx},echo_fn_name{passes_idx}); %develops path for rds data to then load
      fprintf('Loading %s (%s)\n', echo_fn{passes_idx}, datestr(now));
      tmp_data = load_L1B(echo_fn{passes_idx}); %loads data into tmp_data
      if frm_idx == 1
        metadata{passes_idx} = struct('day_seg',passes(passes_idx).day_seg,...
          'frms',passes(passes_idx).frms,'param_records', tmp_data.param_records,...
          'param_sar',tmp_data.param_sar,'param_array',tmp_data.param_array,...
          'time',tmp_data.Time,'param_pass',param_pass); %creates tmp struct with given fields
        
        data{passes_idx} = [];
        metadata{passes_idx}.gps_time = [];
        metadata{passes_idx}.lat = [];
        metadata{passes_idx}.lon = [];
        metadata{passes_idx}.elev = [];
        metadata{passes_idx}.roll = [];
        metadata{passes_idx}.pitch = [];
        metadata{passes_idx}.heading = [];
        metadata{passes_idx}.surface = [];
        metadata{passes_idx}.bottom = [];
        metadata{passes_idx}.fcs.origin = [];
        metadata{passes_idx}.fcs.x = [];
        metadata{passes_idx}.fcs.y = [];
        metadata{passes_idx}.fcs.z = [];
        metadata{passes_idx}.fcs.pos = [];
      end
      metadata{passes_idx}.gps_time = [metadata{passes_idx}.gps_time ,tmp_data.GPS_time];
      metadata{passes_idx}.lat = [metadata{passes_idx}.lat ,tmp_data.Latitude];
      metadata{passes_idx}.lon = [metadata{passes_idx}.lon ,tmp_data.Longitude];
      metadata{passes_idx}.elev = [metadata{passes_idx}.elev ,tmp_data.Elevation];
      metadata{passes_idx}.roll = [metadata{passes_idx}.roll ,tmp_data.Roll];
      metadata{passes_idx}.pitch = [metadata{passes_idx}.pitch ,tmp_data.Pitch];
      metadata{passes_idx}.heading = [metadata{passes_idx}.heading ,tmp_data.Heading];
      metadata{passes_idx}.surface = [metadata{passes_idx}.surface ,tmp_data.Surface];
      metadata{passes_idx}.bottom = [metadata{passes_idx}.bottom ,tmp_data.Bottom];
      % Handle old format, but print error message
      if isfield(tmp_data.param_array,'array_param')
        warning('OLD SAR DATA FORMAT. IF MULTIPLE WF-ADC PAIRS WERE USED IN THE DATA PRODUCT, THE FLIGHT COORDINATE SYSTEM WILL BE INCORRECT AND THE DATA SHOULD BE REPROCESSED.');
        metadata{passes_idx}.fcs.origin = [metadata{passes_idx}.fcs.origin ,tmp_data.param_array.array_param.fcs{1}{1}.origin];
        metadata{passes_idx}.fcs.x = [metadata{passes_idx}.fcs.x ,tmp_data.param_array.array_param.fcs{1}{1}.x];
        metadata{passes_idx}.fcs.y = [metadata{passes_idx}.fcs.y ,tmp_data.param_array.array_param.fcs{1}{1}.y];
        metadata{passes_idx}.fcs.z = [metadata{passes_idx}.fcs.z ,tmp_data.param_array.array_param.fcs{1}{1}.z];
        metadata{passes_idx}.fcs.pos = [metadata{passes_idx}.fcs.pos ,tmp_data.param_array.array_param.fcs{1}{1}.pos];
      elseif iscell(tmp_data.param_array.array_proc.fcs)
        warning('OLD SAR DATA FORMAT. IF MULTIPLE WF-ADC PAIRS WERE USED IN THE DATA PRODUCT, THE FLIGHT COORDINATE SYSTEM WILL BE INCORRECT AND THE DATA SHOULD BE REPROCESSED.');
        metadata{passes_idx}.fcs.origin = [metadata{passes_idx}.fcs.origin ,tmp_data.param_array.array_proc.fcs{1}{1}.origin];
        metadata{passes_idx}.fcs.x = [metadata{passes_idx}.fcs.x ,tmp_data.param_array.array_proc.fcs{1}{1}.x];
        metadata{passes_idx}.fcs.y = [metadata{passes_idx}.fcs.y ,tmp_data.param_array.array_proc.fcs{1}{1}.y];
        metadata{passes_idx}.fcs.z = [metadata{passes_idx}.fcs.z ,tmp_data.param_array.array_proc.fcs{1}{1}.z];
        metadata{passes_idx}.fcs.pos = [metadata{passes_idx}.fcs.pos ,tmp_data.param_array.array_proc.fcs{1}{1}.pos];
      else
        metadata{passes_idx}.fcs.origin = [metadata{passes_idx}.fcs.origin ,tmp_data.param_array.array_proc.fcs.origin];
        metadata{passes_idx}.fcs.x = [metadata{passes_idx}.fcs.x ,tmp_data.param_array.array_proc.fcs.x];
        metadata{passes_idx}.fcs.y = [metadata{passes_idx}.fcs.y ,tmp_data.param_array.array_proc.fcs.y];
        metadata{passes_idx}.fcs.z = [metadata{passes_idx}.fcs.z ,tmp_data.param_array.array_proc.fcs.z];
        metadata{passes_idx}.fcs.pos = [metadata{passes_idx}.fcs.pos ,tmp_data.param_array.array_proc.fcs.pos];
      end
      data{passes_idx} = [data{passes_idx} ,tmp_data.Data];
    end
  else
    %% Load: Load "sar" Data
    
    % Setup param_sar structure used to load SAR data
    param_sar = param_pass;
    % (wf,adc) pairs to load
    param_sar.load_sar_data.imgs = passes(passes_idx).imgs;
    % Debug level (1 = default)
    param_sar.load_sar_data.debug_level = 2;
    
    for frm_idx = 1:length(passes(passes_idx).frms)
      % Load SAR data for this frame
      param_sar.load_sar_data.frm = passes(passes_idx).frms(frm_idx);
      % load_sar_data loads all images at once
      [data{end+1},metadata{end+1}] = load_sar_data(param_sar);
      
      metadata{end}.frms = passes(passes_idx).frms(frm_idx);
      metadata{end}.param_pass = param_pass;
      
      %% Do image combining
      % Combines low-gain and high-gain images into a single image
      param_mode = 'array';
      param_img_combine = param_pass;
      % array.img_comb should have 3*(length(imgs)-1) fields, if it has more
      % than required, then truncate to 3*(length(imgs)-1)
      param_img_combine.array.img_comb = param_pass.array.img_comb(1:...
        min([length(param_pass.array.img_comb), 3*(length(passes(passes_idx).imgs)-1)]));
      param_img_combine.array.imgs = passes(passes_idx).imgs;
      if length(param_img_combine.array.imgs) > 1 && length(param_img_combine.array.img_comb) == 2*(length(param_img_combine.array.imgs)-1)
        error('Spreadsheet has the wrong number of entries for param.array.img_comb. It has the right number for the old combine method. Usually this is fixed by inserting "-inf" for the second coefficient which controls how the receiver blanking is used. For example [3e-6 1e-6; 10e-6 3e-6] becomes [3e-6 -inf 1e-6; 10e-6 -inf 3e-6].');
      end
      param_img_combine.load.frm = metadata{end}.frms;
      param_img_combine.day_seg = passes(passes_idx).day_seg;
      
      data_in = [];
      data_in.Data =data{end};
      data_in.Time = {};
      for img = 1:length(passes(passes_idx).imgs)
        wf = passes(passes_idx).imgs{img}(1);
        data_in.Time{img}= metadata{end}.wfs(wf).time;
      end
      data_in.GPS_time = metadata{end}.fcs{1}{1}.gps_time;
      if length(data{end})>1 %Combine if using multiple waveforms
        [data{end}, metadata{end}.time] = img_combine(param_img_combine,param_mode,surf_layer,data_in);
      else
        data{end} = data{end}{1};
      end
      if 0
        %Check image combine output
        figure(1); clf; h_axes = [];
        for idx = 1:length(data_in.Data)
          h_axes(idx) = subplot(length(data_in.Data),1,idx);
          imagesc(data_in.GPS_time,data_in.Time{idx},lp(data_in.Data{idx}))
          xlabel('GPS Time')
          ylabel('TWTT')
        end
        
        figure(2); clf;
        h_axes(end+1) = axes('parent',2);
        imagesc(data_in.GPS_time,metadata{end}.time,lp(data{end}))
        xlabel('GPS Time')
        ylabel('TWTT')
        linkaxes(h_axes)
        keyboard
      end
    end
  end
end
clear tmp_data tmp_fn frm_idx passes_idx

%% Extract passes from data
% =========================================================================
% NOTE: This code looks through each pass that is specified and extracts
% out the pass (or passes if there is more than one). This code handles the
% case where a single pass contains 0, 1, 2 or more valid passes. Each
% valid pass will be extracted out so the output "pass" structure array may
% contain more or less entries than the input data/metadata cell arrays. If
% each input cell contains a single pass, then "pass" will be the same
% length as data/metadata with a one-to-one relationship in the same order.

% Get start/stop locations in earth centered earth fixed coordinates
[start.x,start.y,start.z] = geodetic2ecef(start.lat/180*pi,start.lon/180*pi,0,WGS84.ellipsoid);
[stop.x,stop.y,stop.z] = geodetic2ecef(stop.lat/180*pi,stop.lon/180*pi,0,WGS84.ellipsoid);

% Initialize output "pass" variable
pass = [];

for data_idx = 1:length(data)
  frm_id = sprintf('%s_%03d', metadata{data_idx}.param_records.day_seg, metadata{data_idx}.frms(1));
  if strcmpi(param.combine_passes.input_type,'echo')
    fprintf('Input Pass: %s %d to %d\n', frm_id);
  else
    wf = passes(data_idx).imgs{1}(1,1);
    adc = passes(data_idx).imgs{1}(1,2);
    fprintf('Input Pass: %s wf %d adc %d\n', frm_id, wf, adc);
  end
  
  % Find the distance to the start
  start_ecef = [start.x;start.y;start.z];
  stop_ecef = [stop.x;stop.y;stop.z];
  radar_ecef = [];
  [radar_ecef.x,radar_ecef.y,radar_ecef.z] = geodetic2ecef(metadata{data_idx}.lat/180*pi, ...
    metadata{data_idx}.lon/180*pi,0*metadata{data_idx}.elev, ...
    WGS84.ellipsoid);
  radar_ecef = [radar_ecef.x; radar_ecef.y; radar_ecef.z];
  
  %% Extract: Collect the closest point every time the trajectory passes near (<dist_min) the start point
  dist = bsxfun(@minus, radar_ecef, start_ecef);
  dist = sqrt(sum(abs(dist).^2));
  
  start_idxs = [];
  start_points = dist < dist_min; % Find all radar points within dist_min from start
  start_idx = find(start_points,1); % Get the index of the first point on the trajectory that is within dist_min
  while ~isempty(start_idx)
    stop_idx = find(start_points(start_idx:end)==0,1); % Get the first point past the start point that is outside of dist_min
    if isempty(stop_idx)
      [~,new_idx] = min(dist(start_idx:end)); % Within the first section of the trajectory that is less than dist_min, find the index of the minimum point
      new_idx = new_idx + start_idx-1; % Convert it to absolute index
      start_idxs = [start_idxs new_idx]; % Add this index to the start_idxs array
      start_idx = []; % If there is no point past the outside, then terminate
    else
      [~,new_idx] = min(dist(start_idx+(0:stop_idx-1))); % Within the first section of the trajectory that is less than dist_min, find the index of the minimum point
      new_idx = new_idx + start_idx-1; % Convert it to absolute index
      start_idxs = [start_idxs new_idx]; % Add this index to the start_idxs array
      new_start_idx = find(start_points(start_idx+stop_idx-1:end),1); % Find the next passby of the start point
      start_idx = new_start_idx + start_idx+stop_idx-1-1; % Convert it to absolute index
    end
  end
  
  %% Extract: Collect the closest point every time the trajectory passes near (<dist_min) the stop point
  stop_dist = bsxfun(@minus, radar_ecef, stop_ecef);
  stop_dist = sqrt(sum(abs(stop_dist).^2));
  
  stop_idxs = [];
  start_points = stop_dist < dist_min;
  start_idx = find(start_points,1);
  while ~isempty(start_idx) % This loop works in the same way as previous "start_idxs" loop
    stop_idx = find(start_points(start_idx:end)==0,1);
    if isempty(stop_idx)
      [~,new_idx] = min(stop_dist(start_idx:end));
      new_idx = new_idx + start_idx-1;
      stop_idxs = [stop_idxs new_idx];
      start_idx = [];
    else
      [~,new_idx] = min(stop_dist(start_idx+(0:stop_idx-1)));
      new_idx = new_idx + start_idx-1;
      stop_idxs = [stop_idxs new_idx];
      new_start_idx = find(start_points(start_idx+stop_idx-1:end),1);
      start_idx = new_start_idx + start_idx+stop_idx-1-1;
    end
  end
  
  if 0
    plot(dist,'b');
    hold on;
    plot(start_idxs, dist(start_idxs), 'ro');
    plot(stop_dist,'k');
    plot(stop_idxs, stop_dist(stop_idxs), 'ro');
    hold off;
    pause;
  end
  
  %% Extract: extract the data out of each pass in this frame sequence
  idxs = [start_idxs stop_idxs]; % Concatenate into one long 1 by N array
  [idxs,sort_idxs] = sort(idxs); % Sort the array
  start_mask = [ones(size(start_idxs)) zeros(size(stop_idxs))]; % Create another 1 by N array that indicates which indices are start_idxs
  start_mask = start_mask(sort_idxs);
  no_passes_flag = true;
  
  for pass_idx = 2:length(idxs)
    if start_mask(pass_idx) ~= start_mask(pass_idx-1) % If we have a start then stop or stop then start, we assume this is a SAR "pass"
      start_idx = idxs(pass_idx-1); % Get the first index of this pass
      stop_idx = idxs(pass_idx);% Get the last index of this pass
      no_passes_flag = false;
      
      fprintf('  Valid Pass: %s %d to %d\n', frm_id, start_idx, stop_idx);
      
      %% Extract: Extract the pass and save it
      if start_mask(pass_idx-1)
        rlines = start_idx:stop_idx;
        pass(end+1).direction = 1;
      else
        rlines = stop_idx:-1:start_idx;
        pass(end+1).direction = -1;
      end
      
      pass(end).wf = 1;
      if strcmpi(param.combine_passes.input_type,'echo')
        pass(end).data = data{data_idx}(:,rlines);
        pass(end).wfs.time = metadata{data_idx}.time;
        if isfield(metadata{data_idx}.param_array,'combine')
          % Old format
          if iscell(metadata{data_idx}.param_array.combine.imgs{1})
            wf = metadata{data_idx}.param_array.combine.imgs{1}{1}(1);
          else
            wf = metadata{data_idx}.param_array.combine.imgs{1}(1);
          end
        else
          if iscell(metadata{data_idx}.param_array.array.imgs{1})
            wf = metadata{data_idx}.param_array.array.imgs{1}{1}(1);
          else
              % Old format
            wf = metadata{data_idx}.param_array.array.imgs{1}(1);
          end
        end
        if ~isfield(metadata{data_idx}.param_array.radar.wfs,'fc')
          % Old format
          pass(end).wfs.fc = 0.5*(metadata{data_idx}.param_array.radar.wfs(wf).f0+metadata{data_idx}.param_array.radar.wfs(wf).f1);
        else
          pass(end).wfs.fc = metadata{data_idx}.param_array.radar.wfs(wf).fc;
        end
        pass(end).time = metadata{data_idx}.time;
        
        pass(end).gps_time = metadata{data_idx}.gps_time(rlines);
        pass(end).roll = metadata{data_idx}.roll(rlines);
        pass(end).pitch = metadata{data_idx}.pitch(rlines);
        pass(end).heading = metadata{data_idx}.heading(rlines);
        
        pass(end).x = metadata{data_idx}.fcs.x(:,rlines);
        pass(end).y = metadata{data_idx}.fcs.y(:,rlines);
        pass(end).z = metadata{data_idx}.fcs.z(:,rlines);
        pass(end).origin = metadata{data_idx}.fcs.origin(:,rlines);
        if size(metadata{data_idx}.fcs.pos,2) ~= size(metadata{data_idx}.fcs.origin,2)
          warning('Old file format with pos field error. No valid data in pos field. Setting to all zeros (i.e. the reference position of the array and not the actual position of the wf-adc pair).');
          pass(end).pos = zeros(3,length(rlines));
        else
          pass(end).pos = metadata{data_idx}.fcs.pos(:,rlines);
        end
        pass(end).surface = metadata{data_idx}.surface(:,rlines);
      else
        pass(end).data = data{data_idx}(:,rlines);
        pass(end).wfs = metadata{data_idx}.wfs;
        pass(end).time = metadata{data_idx}.time;
        
        pass(end).gps_time = metadata{data_idx}.fcs{1}{1}.gps_time(rlines);
        pass(end).roll = metadata{data_idx}.fcs{1}{1}.roll(rlines);
        pass(end).pitch = metadata{data_idx}.fcs{1}{1}.pitch(rlines);
        pass(end).heading = metadata{data_idx}.fcs{1}{1}.heading(rlines);
        
        pass(end).x = metadata{data_idx}.fcs{1}{1}.x(:,rlines);
        pass(end).y = metadata{data_idx}.fcs{1}{1}.y(:,rlines);
        pass(end).z = metadata{data_idx}.fcs{1}{1}.z(:,rlines);
        pass(end).origin = metadata{data_idx}.fcs{1}{1}.origin(:,rlines);
        pass(end).pos = metadata{data_idx}.fcs{1}{1}.pos(:,rlines);
        pass(end).surface = metadata{data_idx}.fcs{1}{1}.surface(:,rlines);
      end
      
      pass(end).lat = metadata{data_idx}.lat(rlines);
      pass(end).lon = metadata{data_idx}.lon(rlines);
      pass(end).elev = metadata{data_idx}.elev(rlines);
      
      pass(end).param_records = metadata{data_idx}.param_records;
      pass(end).param_sar = metadata{data_idx}.param_sar;
      pass(end).param_pass = metadata{data_idx}.param_pass;
      pass(end).param_pass.cmd.frms = metadata{data_idx}.frms;
      pass(end).input_type = param.combine_passes.input_type;
      
    end
  end
  
  if no_passes_flag
    warning('Frame %s_%03d has no passes. Closest distance from start %.0f m. Closest distance from stop %.0f m.', ...
      metadata{data_idx}.param_sar.day_seg, metadata{data_idx}.frms(1), min(dist), min(stop_dist));
  end
  
end

%% Save results
% =========================================================================
out_fn = fullfile(ct_filename_out(param,param.combine_passes.out_path,'',1),[pass_name '.mat']);
fprintf('Saving %s (%s)\n', out_fn, datestr(now));
out_fn_dir = fileparts(out_fn);
if ~exist(out_fn_dir,'dir')
  mkdir(out_fn_dir);
end
param_combine_passes = param;
save(out_fn,'-v7.3','pass','param_combine_passes');
