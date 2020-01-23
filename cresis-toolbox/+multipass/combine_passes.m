% function combine_passes(param,param_override)
% combine_passes(param,param_override)
%
% 

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.combine_passes.pass_name, datestr(now));
fprintf('=====================================================================\n');

%% Input checking

dist_min = param.combine_passes.dist_min;
passes = param.combine_passes.passes;
start = param.combine_passes.start;
stop = param.combine_passes.stop;

if ~isfield(param.combine_passes,'layer') || isempty(param.combine_passes.layer)
  param.combine_passes.layer = struct();
end
if length(param.combine_passes.layer) < 2
  param.combine_passes.layer(2) = struct();
end
if ~isfield(param.combine_passes.layer(1),'name') || isempty(param.combine_passes.layer(1).name)
  param.combine_passes.layer(1).name = 'surface';
end
if ~isfield(param.combine_passes.layer(1),'source') || isempty(param.combine_passes.layer(1).source)
  param.combine_passes.layer(1).source = 'layerData';
end
if ~isfield(param.combine_passes.layer(2),'name') || isempty(param.combine_passes.layer(2).name)
  param.combine_passes.layer(2).name = 'bottom';
end
if ~isfield(param.combine_passes.layer(1),'source') || isempty(param.combine_passes.layer(2).source)
  param.combine_passes.layer(2).source = 'layerData';
end

%% Load echogram data
%loads data for chosen seasons and frames listed in run_combine passes.
%Creates tmp struct with chosen parameter fields. 
%Adapated from combine_passes SAR processing script.

%Authors: John Paden, Cody Barnett, Bailey Miller  


%loop for loading data
layers = {}; Datacomb = {}; timecomb = {};
metadata = {};
data = {};
for passes_idx = 1: length(passes)
  %% Load param xls
   param_fn = ct_filename_param(passes(passes_idx).param_fn); %gets filename
   if passes_idx == 1 || ~strcmp(passes(passes_idx).day_seg,passes(passes_idx-1).day_seg)
     param_multipass = read_param_xls(param_fn,passes(passes_idx).day_seg); %reads parameter sheet for given pass
   end
   
   if passes_idx == master_pass_idx
     param = merge_structs(param_multipass,param);
   end
     %% Load layers
  % Load layers
  layers{passes_idx} = opsLoadLayers(param_multipass,param.combine_passes.layer);

%     % Interpolate all layers onto a common reference (ref)
%     for lay_idx = 1:length(pass(pass_idx).layers)
%       pass(pass_idx).layers(lay_idx).twtt ...
%         = interp_finite(interp1(pass(pass_idx).layers(lay_idx).gps_time, ...
%         pass(pass_idx).layers(lay_idx).twtt, ...
%         pass(pass_idx).gps_time,'linear'));
%     end
   %% Format Data
   param_multipass = merge_structs(param_multipass, param_override); %merges param_multipass and param_override into one struct
   
   if strcmp(param.combine_passes.echo_sar,'echo') || isempty(param.combine_passes.echo_sar)
     echo_fn_dir{passes_idx} = ct_filename_out(param_multipass, passes(passes_idx).in_path); %creates directory for pass, from given data format
     for frm_idx = 1:length(passes(passes_idx).frms) %loop for individual frame loading
        echo_fn_name{passes_idx} = sprintf('Data_%s_%03d.mat',passes(passes_idx).day_seg,passes(passes_idx).frms(frm_idx));
        echo_fn{passes_idx} = fullfile(echo_fn_dir{passes_idx},echo_fn_name{passes_idx}); %develops path for rds data to then load
        fprintf('Loading %s (%s)\n', echo_fn{passes_idx}, datestr(now));
        tmp_data = load_L1B(echo_fn{passes_idx}); %loads data into tmp_data file
        if frm_idx == 1
          metadata{passes_idx} = struct('day_seg',passes(passes_idx).day_seg,...
            'frms',passes(passes_idx).frms,'param_records', tmp_data.param_records,...
            'param_sar',tmp_data.param_sar,'param_array',tmp_data.param_array,...
            'time',tmp_data.Time,'param_multipass',param_multipass); %creates tmp struct with given fields

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
          metadata{passes_idx}.param_array =[];
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
        metadata{passes_idx}.fcs.origin = [metadata{passes_idx}.fcs.origin ,tmp_data.param_array.array_proc.fcs.origin];
        metadata{passes_idx}.fcs.x = [metadata{passes_idx}.fcs.x ,tmp_data.param_array.array_proc.fcs.x];
        metadata{passes_idx}.fcs.y = [metadata{passes_idx}.fcs.y ,tmp_data.param_array.array_proc.fcs.y];
        metadata{passes_idx}.fcs.z = [metadata{passes_idx}.fcs.z ,tmp_data.param_array.array_proc.fcs.z];      
        metadata{passes_idx}.fcs.pos = [metadata{passes_idx}.fcs.pos ,tmp_data.param_array.array_proc.fcs.pos];
        data{passes_idx} = [data{passes_idx} ,tmp_data.Data];
  %       if passes_idx==1 && frm_idx==1 %condition to make first field in struct, if developed moves to else statement and adds to end
  %         rds_data(passes_idx) = tmp_struct;
  %       else
  %         rds_data(end+1) = tmp_struct;
  %       end
      end
   elseif strcmp(param.combine_passes.echo_sar,'sar')
%       param_sar = [];
%       param_sar.day_seg = passes(passes_idx).day_seg;
%       param_sar = read_param_xls(ct_filename_param(passes(passes_idx).param_fn),param_sar.day_seg);
      param_sar = param_multipass;

      param_sar.load_sar_data.fn = ''; % Leave empty for default

      % Start and stop chunk to load (inf for second element loads to the end)
      param_sar.load_sar_data.chunk = [1 inf];

      param_sar.load_sar_data.sar_type = 'fk';
      
      param_sar.load_sar_data.subap = 1;

      % (wf,adc) pairs to load
      if ~iscell(passes(passes_idx).wf_adc)
        passes(passes_idx).wf_adc = {passes(passes_idx).wf_adc};
      end
      param_sar.load_sar_data.imgs = passes(passes_idx).wf_adc;

      % Combine waveforms parameters
      param_sar.load_sar_data.wf_comb = 10e-6;

      % Debug level (1 = default)
      param_sar.load_sar_data.debug_level = 2;

      % Combine receive channels
      param_sar.load_sar_data.combine_channels = 0;

      % Take abs()^2 of the data (only runs if combine_channels runs)
      param_sar.load_sar_data.incoherent = 0;

      % Combine waveforms (only runs if incoherent runs)
      param_sar.load_sar_data.combine_waveforms = 0;

      % Parameters for local_detrend (cmd == 5 disables, only runs if incoherent runs)
      param_sar.load_sar_data.detrend.cmd = 3;
      param_sar.load_sar_data.detrend.B_noise = [100 200];
      param_sar.load_sar_data.detrend.B_sig = [1 10];
      param_sar.load_sar_data.detrend.minVal = -inf;
      
      for frm_idx = 1:length(passes(passes_idx).frms)
        param_sar.load_sar_data.frame = passes(passes_idx).frms(frm_idx);

        [data{end+1},metadata{end+1}] = load_sar_data(param_sar);

        metadata{end}.frms = passes(passes_idx).frms(frm_idx);
        metadata{end}.param_multipass = param_multipass;
        
        %% Do waveform combination
        param_mode = 'array';
        param_img_combine = param_multipass;
        param_img_combine.array.out_path = 'CSARP_post/standard';
        param_img_combine.sar.out_path = 'sar';
        param_img_combine.sar.img_comb = param_img_combine.array.img_comb;
        param_img_combine.array.img_comb = param_multipass.array.img_comb(1:...
          min([length(param_multipass.array.img_comb), 3*(length(passes(passes_idx).wf_adc)-1)]));
        param_img_combine.array.imgs = passes(passes_idx).wf_adc;
        if length(param_img_combine.array.img_comb)<3
          param_img_combine.array.img_comb = [param_img_combine.array.img_comb(1) -inf param_img_combine.array.img_comb(2)];
        end
        param_img_combine.load.frm = metadata{end}.frms;
        param_img_combine.day_seg = passes(passes_idx).day_seg;
        
        timevecs ={};
        for wf_id = 1:length(passes(passes_idx).wf_adc)
          timevecs{wf_id}= metadata{end}.wfs(wf_id).time;
        end
        data_in.Data =data{end};
        data_in.Time = timevecs;
        data_in.GPS_time = {metadata{end}.fcs{1}{1}.gps_time, metadata{end}.fcs{2}{1}.gps_time};
        if length(data{end})>1 %Combine if using multiple waveforms
          [data{end}, timecomb{end+1}] = img_combine(param_img_combine,param_mode,layers{end}(1),data_in);
        else
          data{end} = data{end}{1};
        end
        if 0
          %Check image combine output
          figure(1)
          subplot(2,1,1)
          imagesc(data_in.GPS_time{1},data_in.Time{1},lp(data_in.Data{1}))
          xlabel('GPS Time')
          ylabel('TWTT')
          subplot(2,1,2)
          imagesc(data_in.GPS_time{1},data_in.Time{1},lp(data_in.Data{2}))
          xlabel('GPS Time')
          ylabel('TWTT')

          figure(2)
          imagesc(data_in.GPS_time{1},timecomb{end},lp(data{end}))
          xlabel('GPS Time')
          ylabel('TWTT')
          keyboard
        end
      end
   end
end 
clear tmp_data tmp_fn frm_idx passes_idx

%% Find start/stop points and extract radar passes
physical_constants;
[start.x,start.y,start.z] = geodetic2ecef(start.lat/180*pi,start.lon/180*pi,0,WGS84.ellipsoid);
[stop.x,stop.y,stop.z] = geodetic2ecef(stop.lat/180*pi,stop.lon/180*pi,0,WGS84.ellipsoid);

pass = [];

%% Go through each frame and extract the pass(es) from that frame
% NOTE: This code looks for every pass in the frame (i.e. a frame may
% contain multiple passes and this code should find each).
for passes_idx = 1:length(passes)
  % Find the distance to the start
  start_ecef = [start.x;start.y;start.z];
  stop_ecef = [stop.x;stop.y;stop.z];
  radar_ecef = [];
  [radar_ecef.x,radar_ecef.y,radar_ecef.z] = geodetic2ecef(metadata{passes_idx}.lat/180*pi, ...
    metadata{passes_idx}.lon/180*pi,0*metadata{passes_idx}.elev, ...
    WGS84.ellipsoid);
  radar_ecef = [radar_ecef.x; radar_ecef.y; radar_ecef.z];
  
  %% Collect the closest point every time the trajectory passes near (<dist_min) the start point
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
  
  %% Collect the closest point every time the trajectory passes near (<dist_min) the stop point
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
  
  %% Extract the data out of each pass in this frame
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
      
      frm_id = sprintf('%s_%03d', metadata{passes_idx}.param_records.day_seg, metadata{passes_idx}.frms(1));
      
      fprintf('New Segment: %s %d to %d\n', frm_id, start_idx, stop_idx);
  
      %% Extract the pass and save it
      if start_mask(pass_idx-1)
        rlines = start_idx:stop_idx;
        pass(end+1).direction = 1;
      else
        rlines = stop_idx:-1:start_idx;
        pass(end+1).direction = -1;
      end
      
      pass(end).wf = 1;
      if strcmp(param.combine_passes.echo_sar,'echo') || isempty(param.combine_passes.echo_sar)
        pass(end).data = data{passes_idx}(:,rlines);
        pass(end).wfs.time = metadata{passes_idx}.time;
        pass(end).time = metadata{passes_idx}.time;
        
        pass(end).gps_time = metadata{passes_idx}.gps_time(rlines);
        pass(end).roll = metadata{passes_idx}.roll(rlines);
        pass(end).pitch = metadata{passes_idx}.pitch(rlines);
        pass(end).heading = metadata{passes_idx}.heading(rlines);
        
        pass(end).x = metadata{passes_idx}.fcs.x(:,rlines);
        pass(end).y = metadata{passes_idx}.fcs.y(:,rlines);
        pass(end).z = metadata{passes_idx}.fcs.z(:,rlines);
        pass(end).origin = metadata{passes_idx}.fcs.origin(:,rlines);
        pass(end).pos = metadata{passes_idx}.fcs.pos(:,rlines);
        pass(end).surface = metadata{passes_idx}.surface(:,rlines);
      elseif strcmp(param.combine_passes.echo_sar,'sar')
        pass(end).data = data{passes_idx}(:,rlines);
        pass(end).wfs = metadata{passes_idx}.wfs;
        pass(end).time = timecomb{passes_idx};
        
        pass(end).gps_time = metadata{passes_idx}.fcs{1}{1}.gps_time(rlines);
        pass(end).roll = metadata{passes_idx}.fcs{1}{1}.roll(rlines);
        pass(end).pitch = metadata{passes_idx}.fcs{1}{1}.pitch(rlines);
        pass(end).heading = metadata{passes_idx}.fcs{1}{1}.heading(rlines);
        
        pass(end).x = metadata{passes_idx}.fcs{1}{1}.x(:,rlines);
        pass(end).y = metadata{passes_idx}.fcs{1}{1}.y(:,rlines);
        pass(end).z = metadata{passes_idx}.fcs{1}{1}.z(:,rlines);
        pass(end).origin = metadata{passes_idx}.fcs{1}{1}.origin(:,rlines);
        pass(end).pos = metadata{passes_idx}.fcs{1}{1}.pos(:,rlines);
        pass(end).surface = metadata{passes_idx}.fcs{1}{1}.surface(:,rlines);
      end
      
      pass(end).lat = metadata{passes_idx}.lat(rlines);
      pass(end).lon = metadata{passes_idx}.lon(rlines);
      pass(end).elev = metadata{passes_idx}.elev(rlines);
      
      pass(end).param_records = metadata{passes_idx}.param_records;
      pass(end).param_sar = metadata{passes_idx}.param_sar;
      pass(end).param_multipass = metadata{passes_idx}.param_multipass;
      pass(end).param_multipass.cmd.frms = passes(passes_idx).frms;
      pass(end).combine_passes = passes(passes_idx);
      pass(end).echo_sar = param.combine_passes.echo_sar;
      pass(end).layers = layers{passes_idx};
      
      
      
    end
  end
  if no_passes_flag
    warning('Frame %s_%03d has no passes. Closest distance from start %.0f m. Closest distance from stop %.0f m.', metadata{passes_idx}.param_sar.day_seg, metadata{passes_idx}.frms(1), min(dist), min(stop_dist));
  end
  
end

%% Save the results
out_fn = fullfile(ct_filename_out(param,'multipass','',1),[param.combine_passes.pass_name '.mat']);
fprintf('  Saving %s\n', out_fn);
out_fn_dir = fileparts(out_fn);
if ~exist(out_fn_dir,'dir')
  mkdir(out_fn_dir);
end
param_combine_passes = param;
save(out_fn,'-v7.3','pass','param_combine_passes');
