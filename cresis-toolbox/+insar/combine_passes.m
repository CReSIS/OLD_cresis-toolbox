%% User Settings
passes = struct('frm',{},'wf_adc',{});

if 1
  %% OIB P3 Greenland: North Line
  param_fn = 'rds_param_2014_Greenland_P3.xls';
  pass_name='rds_thule';
  % Found using check_region
  wf = 2;
%   passes(end+1) = struct('frm','20140429_01_005','wf_adc',[wf 2]);
  for adc = 2:16
    passes(end+1) = struct('frm','20140429_01_067','wf_adc',[wf adc]);
  end
%   passes(end+1) = struct('frm','20140514_01_004','wf_adc',[wf 2]);
%   passes(end+1) = struct('frm','20140514_01_066','wf_adc',[wf 2]);
%   passes(end+1) = struct('frm','20140515_02_004','wf_adc',[wf 2]);
%   passes(end+1) = struct('frm','20140515_02_069','wf_adc',[wf 2]);
%   passes(end+1) = struct('frm','20140521_02_033','wf_adc',[wf 2]);

%   frms{end+1} = '20140429_01_005';
%   frms{end+1} = '20140429_01_067';
%   %   frms{end+1} = '20140429_01_004'; %sar processed stop is far
%   %   frms{end+1} = '20140429_01_066'; %outside of window
%   %   frms{end+1} = '20140501_01_044'; %sar processed %outside of window
%   %   frms{end+1} = '20140501_01_045'; %sar processed %Large baseline
%   %   frms{end+1} = '20140502_01_060'; %High altitude flight
%   %   frms{end+1} = '20140507_01_068'; %outside of window
%   frms{end+1} = '20140514_01_004';
%   frms{end+1} = '20140514_01_066';
%   frms{end+1} = '20140515_02_004';
%   frms{end+1} = '20140515_02_069';
%   %   frms{end+1} = '20140521_01_003'; %Lots of random noise
%   frms{end+1} = '20140521_02_033';
  %Found using check_region
  start = struct('lat', 77.10,'lon', -62.3);
  stop = struct('lat', 77.13, 'lon', -61.9);
  %Set dist min
  dist_min = 300;
elseif 0
  %% TO DTU Greenland: North Line
  param_fn = 'rds_param_2016_Greenland_TOdtu.xls';
  frms{end+1} = '20161107_02_003';
  frms{end+1} = '20161107_02_006';
  frms{end+1} = '20161107_03_003';
  frms{end+1} = '20161107_03_006';
  frms{end+1} = '20161107_07_003';
  frms{end+1} = '20161108_01_003';
  frms{end+1} = '20161108_01_006';
  frms{end+1} = '20161108_01_009';
  frms{end+1} = '20161108_02_002';
  frms{end+1} = '20161110_01_003';
  frms{end+1} = '20161110_01_006';
  frms{end+1} = '20161110_01_009';
  frms{end+1} = '20161110_02_003';
  frms{end+1} = '20161110_03_001';
  pass_name = 'north';
  start.lat = 69.174061;
  start.lon = -49.668096;
  stop.lat = 69.196834;
  stop.lon = -48.906555;
  dist_min = 2000;
elseif 0
  %% TO DTU Greenland: Middle Line
  param_fn = 'rds_param_2016_Greenland_TOdtu.xls';
  % frms{end+1} = '20161107_02_001'; % 11.7 km short of stop point
  frms{end+1} = '20161107_02_004';
  frms{end+1} = '20161107_03_001';
  frms{end+1} = '20161107_03_004';
  frms{end+1} = '20161107_07_001';
  frms{end+1} = '20161108_01_001';
  frms{end+1} = '20161108_01_004';
  frms{end+1} = '20161108_01_007'; % Surface does not show up in interferogram
  frms{end+1} = '20161108_01_010'; % Surface does not show up in interferogram
  frms{end+1} = '20161110_01_001'; % Surface is faint
  frms{end+1} = '20161110_01_004';
  frms{end+1} = '20161110_01_007';
  frms{end+1} = '20161110_02_001'; % 2 km short of stop point
  frms{end+1} = '20161110_02_004';
  frms{end+1} = '20161110_03_007';
  frms{end+1} = '20161110_03_008';
  frms{end+1} = '20161111_01_001'; % 460 m short of start point
  frms{end+1} = '20161111_01_002';
  frms{end+1} = '20161111_01_003';
  frms{end+1} = '20161111_01_004';
  frms{end+1} = '20161111_01_005';
  frms{end+1} = '20161111_01_006';
  frms{end+1} = '20161111_01_007';
  frms{end+1} = '20161111_01_008';
  pass_name = 'middle';
  start.lat = 69.163;
  start.lon = -49.675;
  stop.lat = 69.176;
  stop.lon = -48.888;
  dist_min = 2500;
elseif 0
  %% TO DTU Greenland: South Line
  param_fn = 'rds_param_2016_Greenland_TOdtu.xls';
  frms{end+1} = '20161107_02_002';
  frms{end+1} = '20161107_02_005';
  frms{end+1} = '20161107_03_002';
  frms{end+1} = '20161107_03_005';
  frms{end+1} = '20161107_07_002';
  frms{end+1} = '20161108_01_002';
  frms{end+1} = '20161108_01_005';
  frms{end+1} = '20161108_01_008'; % Surface does not show up in interferogram
  frms{end+1} = '20161108_02_001'; % Surface does not show up in interferogram
  frms{end+1} = '20161110_01_002';
  frms{end+1} = '20161110_01_005';
  frms{end+1} = '20161110_01_008';
  frms{end+1} = '20161110_02_002';
  frms{end+1} = '20161110_02_005';
  pass_name = 'south';
  start.lat = 69.149;
  start.lon = -49.666;
  stop.lat = 69.162;
  stop.lon = -48.846;
  dist_min = 2000;
elseif 0
  %% TO DTU Iceland: South Line
  param_fn = 'rds_param_2016_Greenland_TOdtu.xls';
  frms{end+1} = '20161101_03_004';
  frms{end+1} = '20161101_01_002';
  frms{end+1} = '20161101_02_002';
  frms{end+1} = '20161101_03_002';
  pass_name = 'iceland_south';
  start.lat = 64.802;
  start.lon = -18.850;
  stop.lat = 64.807;
  stop.lon = -19.105;
  dist_min = 2000;
elseif 1
  %% TO DTU Iceland: North Line
  param_fn = 'rds_param_2016_Greenland_TOdtu.xls';
  frms{end+1} = '20161101_01_003';
  frms{end+1} = '20161101_02_001';
  frms{end+1} = '20161101_02_003';
  frms{end+1} = '20161101_03_003';
  frms{end+1} = '20161101_04_001';
  pass_name = 'iceland_north';
  start.lat = 64.843;
  start.lon = -19.203;
  stop.lat = 64.804;
  stop.lon = -18.909;
  dist_min = 2000;  
elseif 0
  %% G1XB Russell Glacier: Good Quality Line
  param_fn = 'rds_param_2016_Greenland_G1XB.xls';
  frms{end+1} = '20160413_01_001';
  frms{end+1} = '20160413_01_002';
  frms{end+1} = '20160413_02_001';
  frms{end+1} = '20160413_02_002';
  frms{end+1} = '20160413_02_003';
  frms{end+1} = '20160413_02_004';
  pass_name = 'good_line';
  start.lat = 67.092809;
  start.lon = -50.204091;
  stop.lat = 67.096958;
  stop.lon = -50.054023;
  dist_min = 100;
elseif 0
  %% G1XB Russell Glacier: Medium Quality Line
  param_fn = 'rds_param_2016_Greenland_G1XB.xls';
  frms{end+1} = '20160416_01_002';
  frms{end+1} = '20160416_01_003';
  frms{end+1} = '20160416_01_004';
  frms{end+1} = '20160416_01_005';
  % frms{end+1} = '20160417_04_002'; % LARGE BASELINE
  % frms{end+1} = '20160417_04_003'; % LARGE BASELINE
  frms{end+1} = '20160417_04_004';
  frms{end+1} = '20160417_04_005';
  pass_name = 'medium_line';
  start.lat = 67.097188;
  start.lon = -50.219048;
  stop.lat = 67.101868;
  stop.lon = -50.047311;
  dist_min = 100;
else
  %% G1XB Russell Glacier: Bad Quality Line
  param_fn = 'rds_param_2016_Greenland_G1XB.xls';
  frms{end+1} = '20160417_01_001';
  % frms{end+1} = '20160417_01_002'; GPS BAD
  frms{end+1} = '20160417_02_002';
  frms{end+1} = '20160417_02_003';
  frms{end+1} = '20160417_02_004';
  frms{end+1} = '20160417_02_005';
  frms{end+1} = '20160417_03_002';
  frms{end+1} = '20160417_03_003';
  frms{end+1} = '20160417_03_004';
  frms{end+1} = '20160417_03_005';
  pass_name = 'bad_line';
  start.lat = 67.102454;
  start.lon = -50.192878;
  stop.lat = 67.108618;
  stop.lon = -49.968087;
  dist_min = 100;
end

%% Automated
for passes_idx = 1:length(passes)
  if ~isfield(passes(passes_idx),'wf_adc') || isempty(passes(passes_idx).wf_adc)
    passes(passes_idx).wf_adc = [1 1];
  end
end

if 1
  % Load SAR data
  metadata = [];
  data = [];
  for passes_idx = 1:length(passes)
    param = [];
    param.day_seg = passes(passes_idx).frm(1:11);
    param = read_param_xls(ct_filename_param(param_fn),param.day_seg);
    
    param.load_sar_data.fn = ''; % Leave empty for default
    
    % Start and stop chunk to load (inf for second element loads to the end)
    param.load_sar_data.chunk = [1 inf];
    
    param.load_sar_data.sar_type = 'fk';
    
    frm = str2double(passes(passes_idx).frm(13:end));
    param.load_sar_data.frame = frm;
    
    param.load_sar_data.subap = 1;
    
    % (wf,adc) pairs to load
    param.load_sar_data.imgs = {passes(passes_idx).wf_adc};
    
    % Combine waveforms parameters
    param.load_sar_data.wf_comb = 10e-6;
    
    % Debug level (1 = default)
    param.load_sar_data.debug_level = 2;
    
    % Combine receive channels
    param.load_sar_data.combine_channels = 0;
    
    % Take abs()^2 of the data (only runs if combine_channels runs)
    param.load_sar_data.incoherent = 0;
    
    % Combine waveforms (only runs if incoherent runs)
    param.load_sar_data.combine_waveforms = 0;
    
    % Parameters for local_detrend (cmd == 5 disables, only runs if incoherent runs)
    param.load_sar_data.detrend.cmd = 3;
    param.load_sar_data.detrend.B_noise = [100 200];
    param.load_sar_data.detrend.B_sig = [1 10];
    param.load_sar_data.detrend.minVal = -inf;
    
    [data{passes_idx},metadata{passes_idx}] = load_sar_data(param);
    
    metadata{passes_idx}.frm = frm;
  end
end

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
      
      frm_id = sprintf('%s_%03d', metadata{passes_idx}.param_sar.day_seg, metadata{passes_idx}.frm);
      
      fprintf('New Segment: %s %d to %d\n', frm_id, start_idx, stop_idx);
  
      %% Extract the pass and save it
      if start_mask(pass_idx-1)
        rlines = start_idx:stop_idx;
        pass(end+1).direction = 1;
      else
        rlines = stop_idx:-1:start_idx;
        pass(end+1).direction = -1;
      end
      
      pass(end).frm = passes(passes_idx).frm;
      pass(end).wf = passes(passes_idx).wf_adc(1);
      pass(end).adc = passes(passes_idx).wf_adc(2);
      pass(end).data = data{passes_idx}{1}(:,rlines);
      
      pass(end).gps_time = metadata{passes_idx}.fcs{1}{1}.gps_time(rlines);
      pass(end).lat = metadata{passes_idx}.lat(rlines);
      pass(end).lon = metadata{passes_idx}.lon(rlines);
      pass(end).elev = metadata{passes_idx}.elev(rlines);
      pass(end).roll = metadata{passes_idx}.fcs{1}{1}.roll(rlines);
      pass(end).pitch = metadata{passes_idx}.fcs{1}{1}.pitch(rlines);
      pass(end).heading = metadata{passes_idx}.fcs{1}{1}.heading(rlines);
      
      pass(end).Lsar = metadata{passes_idx}.fcs{1}{1}.Lsar;
      pass(end).wfs = metadata{passes_idx}.wfs;
      pass(end).param_records = metadata{passes_idx}.param_records;
      pass(end).param_sar = metadata{passes_idx}.param_sar;
      pass(end).surface = metadata{passes_idx}.fcs{1}{1}.surface(:,rlines);
      
      pass(end).x = metadata{passes_idx}.fcs{1}{1}.x(:,rlines);
      pass(end).y = metadata{passes_idx}.fcs{1}{1}.y(:,rlines);
      pass(end).z = metadata{passes_idx}.fcs{1}{1}.z(:,rlines);
      pass(end).origin = metadata{passes_idx}.fcs{1}{1}.origin(:,rlines);
      pass(end).pos = metadata{passes_idx}.fcs{1}{1}.pos(:,rlines);
    end
  end
  if no_passes_flag
    warning('Frame %s_%03d has no passes.', metadata{passes_idx}.param_sar.day_seg, metadata{passes_idx}.frm);
  end
  
end

%% Save the results
out_fn = fullfile(ct_filename_out(param,'insar','',1),[pass_name '.mat']);
fprintf('  Saving %s\n', out_fn);
out_fn_dir = fileparts(out_fn);
if ~exist(out_fn_dir,'dir')
  mkdir(out_fn_dir);
end
param_insar = param;
save(out_fn,'-v7.3','pass','param_insar');
