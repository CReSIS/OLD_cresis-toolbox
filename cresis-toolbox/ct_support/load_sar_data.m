function [data,metadata] = load_sar_data(param)
% [data,metadata] = load_sar_data(param)
%
% Loads and concatenates sar data. Currently requires all data
% to be in the same directory. Returns data, position, and header
% information associated with the sar data.
% Known Issue: Does not remove overlapping parts of chunks
%
% Options for:
%   combine_channels (0 = returns 3D matrix)
%   incoherent = takes abs()^2
%   combine_waveforms = combines waveforms
%   local detrending
%   * no data decimation is done except when combining
%   * does not support multiple tx channels yet
%
% Author: John Paden, Logan Smith
%
% See also: run_master.m, master.m, run_array.m, array.m, load_sar_data.m,
% array_proc.m, array_task.m, array_combine_task.m
%
% Also used in: run_load_sar_data.m

% Path to the input data
year = str2double(param.day_seg(1:4));
month = str2double(param.day_seg(5:6));
day = str2double(param.day_seg(7:8));
seg = str2double(param.day_seg(10:11));

%% Check input arguments
% =====================================================================
if ~isfield(param.load_sar_data,'debug_level')
  param.load_sar_data.debug_level = 1;
end
if ~isfield(param.load_sar_data,'combine_channels')
  param.load_sar_data.combine_channels = 1;
end
if ~isfield(param.load_sar_data,'incoherent')
  param.load_sar_data.incoherent = 1;
end
if ~isfield(param.load_sar_data,'combine_waveforms')
  param.load_sar_data.combine_waveforms = 1;
end
if ~isfield(param.load_sar_data,'detrend') || ~isfield(param.load_sar_data.detrend,'cmd')
  param.load_sar_data.detrend.cmd = 3;
end
if ~isfield(param.load_sar_data.detrend,'B_noise')
  param.load_sar_data.detrend.B_noise = [100 200];
end
if ~isfield(param.load_sar_data.detrend,'B_sig')
  param.load_sar_data.detrend.B_sig = [10 20];
end
if ~isfield(param.load_sar_data,'fn') || isempty(param.load_sar_data.fn)
  param.load_sar_data.fn = 'sar';
end
% The base path for all the data
base_path = ct_filename_out(param,param.load_sar_data.fn,'');
if ~isfield(param.load_sar_data.detrend,'minVal')
  param.load_sar_data.detrend.minVal = -inf;
end
if ~isfield(param.load_sar_data,'wf_comb')
  param.load_sar_data.wf_comb           = 10e-6;
end

physical_constants;

%% Initialize memory for outputs
data     = cell(length(param.load_sar_data.imgs),1);
metadata = [];

for subap_idx = 1:length(param.load_sar_data.subap)
  subap = param.load_sar_data.subap(subap_idx);

  
  if param.load_sar_data.sar_type(1) == 'f'
    % F-k data
    in_path = fullfile(base_path, ...
      sprintf('%s_data_%03d_%02d_01/','fk',param.load_sar_data.frame, subap));
  else
    error('Not supported');
  end
  
  %% Get the list of files for this subaperture set to determine what chunks are available
  img = 1;
  wf_adc_list = param.load_sar_data.imgs{img};
  wf_adc_idx = 1;
  wf = wf_adc_list(wf_adc_idx,1);
  adc = wf_adc_list(wf_adc_idx,2);
  fns = get_filenames(in_path,sprintf('wf_%02.0f_adc_%02.0f_chk_',wf,adc),'','.mat');
  valid_chks = [];
  for fns_idx = 1:length(fns)
    [~,fn_name] = fileparts(fns{fns_idx});
    valid_chks(end+1) = str2double(fn_name(end-2:end));
  end
  param.load_sar_data.chunk(param.load_sar_data.chunk==inf) = max(valid_chks);
  chks_to_load = param.load_sar_data.chunk(1):param.load_sar_data.chunk(end);
  
  %% Remove chunks that do not exist from chunks_to_load list
  [valid_chks,keep_idxs] = intersect(chks_to_load, valid_chks);
  if length(valid_chks) ~= length(chks_to_load)
    bad_mask = ones(size(chks_to_load));
    bad_mask(keep_idxs) = 0;
    warning('Nonexistent chunks specified in chks_to_load (e.g. chunk "%g" is invalid), removing these', ...
      chks_to_load(find(bad_mask,1)));
    chks_to_load = valid_chks;
  end
  
  cur_rline = 0;
  for chunk = chks_to_load
    
    if param.load_sar_data.debug_level >= 1
      fprintf('Loading chunk %d (%s)\n', chunk, datestr(now,'HH:MM:SS'));
    end
    
    % =====================================================================
    % Load data
    
    tx = 1;
    for img = 1:length(param.load_sar_data.imgs)
      wf_adc_list = param.load_sar_data.imgs{img};
      
      % -------------------------------------------------------------------
      % Load data
      for wf_adc_idx = 1:size(wf_adc_list,1)
        wf = wf_adc_list(wf_adc_idx,1);
        adc = wf_adc_list(wf_adc_idx,2);
        
        sar_fn = fullfile(in_path,sprintf('wf_%02.0f_adc_%02.0f_chk_%03.0f.mat',wf,adc,chunk));
        if param.load_sar_data.debug_level >= 2
          fprintf(' %s (%s)\n', sar_fn, datestr(now,'HH:MM:SS'));
        end
        sar_data = load(sar_fn);

        % Only add in non-overlapping part of SAR image
        if img == 1 && wf_adc_idx == 1
          if chunk == param.load_sar_data.chunk(1)
            new_idxs = 1:length(sar_data.fcs.gps_time);
          else
            % 1e-3 added to avoid rounding errors... this is a hack
            new_idxs = find(sar_data.fcs.gps_time > fcs{img}{wf_adc_idx}.gps_time(end)+1e-3);
          end
        end
        
        % Get output image positions (not phase centers)
        if img == 1 && wf_adc_idx == 1 && subap_idx == 1
          if chunk == param.load_sar_data.chunk(1)
            metadata.lat = sar_data.lat;
            metadata.lon = sar_data.lon;
            metadata.elev = sar_data.elev;
          else
            metadata.lat = [metadata.lat sar_data.lat(new_idxs)];
            metadata.lon = [metadata.lon sar_data.lon(new_idxs)];
            metadata.elev = [metadata.elev sar_data.elev(new_idxs)];
          end
        end
        
        if subap_idx == 1
          if chunk == param.load_sar_data.chunk(1)
            fcs{img}{wf_adc_idx} = sar_data.fcs;
            
          else
            
            fcs{img}{wf_adc_idx}.gps_time = cat(2,fcs{img}{wf_adc_idx}.gps_time, ...
              sar_data.fcs.gps_time(new_idxs));
            fcs{img}{wf_adc_idx}.x = cat(2,fcs{img}{wf_adc_idx}.x, ...
              sar_data.fcs.x(:,new_idxs));
            fcs{img}{wf_adc_idx}.y = cat(2,fcs{img}{wf_adc_idx}.y, ...
              sar_data.fcs.y(:,new_idxs));
            fcs{img}{wf_adc_idx}.z = cat(2,fcs{img}{wf_adc_idx}.z, ...
              sar_data.fcs.z(:,new_idxs));
            fcs{img}{wf_adc_idx}.origin = cat(2,fcs{img}{wf_adc_idx}.origin, ...
              sar_data.fcs.origin(:,new_idxs));
            fcs{img}{wf_adc_idx}.pos = cat(2,fcs{img}{wf_adc_idx}.pos, ...
              sar_data.fcs.pos(:,new_idxs));
            fcs{img}{wf_adc_idx}.roll = cat(2,fcs{img}{wf_adc_idx}.roll, ...
              sar_data.fcs.roll(new_idxs));
            fcs{img}{wf_adc_idx}.pitch = cat(2,fcs{img}{wf_adc_idx}.pitch, ...
              sar_data.fcs.pitch(new_idxs));
            fcs{img}{wf_adc_idx}.heading = cat(2,fcs{img}{wf_adc_idx}.heading, ...
              sar_data.fcs.heading(new_idxs));
            fcs{img}{wf_adc_idx}.surface = cat(2,fcs{img}{wf_adc_idx}.surface, ...
              sar_data.fcs.surface(new_idxs));
            fcs{img}{wf_adc_idx}.bottom = cat(2,fcs{img}{wf_adc_idx}.bottom, ...
              sar_data.fcs.bottom(new_idxs));
          end
        end
        
        Nt = size(sar_data.fk_data,1);
        if ~param.load_sar_data.combine_channels
          if subap_idx == 1 && wf_adc_idx == 1
            % Allocate memory in a special way when loading 3D data
            data{img}(size(sar_data.fk_data,1), ...
              size(data{img},2)+length(new_idxs),size(wf_adc_list,1),length(param.load_sar_data.subap)) = single(0);
          end
          data{img}(:,cur_rline + (1:length(new_idxs)),wf_adc_idx,subap_idx) = sar_data.fk_data(:,new_idxs);
        else
          % When combining channels, take the mean of the data as it
          % is loaded in to reduce peak memory consumption.
          if wf_adc_idx == 1
            sar_data_data = sar_data.fk_data(:,new_idxs) / size(wf_adc_list,1);
          else
            sar_data_data = sar_data_data + sar_data.fk_data(:,new_idxs) / size(wf_adc_list,1);
          end
        end
        
      end
      
      if param.load_sar_data.combine_channels
        if param.load_sar_data.incoherent
          data{img} = [data{img} abs(sar_data_data).^2];
        else
          data{img} = [data{img} sar_data_data];
        end
      end
    end
    
    cur_rline = cur_rline + length(new_idxs);
  end
end

for img = 1:length(fcs)
  for wf_adc_idx = 1:length(fcs{img})
    [fcs{img}{wf_adc_idx}.lat,fcs{img}{wf_adc_idx}.lon,fcs{img}{wf_adc_idx}.elev] ...
      = ecef2geodetic(fcs{img}{wf_adc_idx}.origin(1,:) + sum(fcs{img}{wf_adc_idx}.x.*fcs{img}{wf_adc_idx}.pos), ...
      fcs{img}{wf_adc_idx}.origin(2,:) + sum(fcs{img}{wf_adc_idx}.y.*fcs{img}{wf_adc_idx}.pos), ...
      fcs{img}{wf_adc_idx}.origin(3,:) + sum(fcs{img}{wf_adc_idx}.z.*fcs{img}{wf_adc_idx}.pos), ...
      WGS84.ellipsoid);
    fcs{img}{wf_adc_idx}.lat = fcs{img}{wf_adc_idx}.lat * 180/pi;
    fcs{img}{wf_adc_idx}.lon = fcs{img}{wf_adc_idx}.lon * 180/pi;
  end
end

% We assume these fields are statics and don't change with each loaded file
% so we just copy them once.
metadata.fcs = fcs;
metadata.wfs = sar_data.wfs;
metadata.param_records = sar_data.param_records;
metadata.param_sar = sar_data.param_sar;

if param.load_sar_data.combine_channels && param.load_sar_data.incoherent ...
    && param.load_sar_data.combine_waveforms
  img1.Time = sar_data.wfs(param.load_sar_data.imgs{1}(1,1)).time;
  img2.Time = sar_data.wfs(param.load_sar_data.imgs{2}(1,1)).time;
  img1.Data = data{1};
  img2.Data = data{2};
  
  % ==============================================================
  % Following is copied from combine_wf_chan (except param.wf_comb):
  dt = img1.Time(2)-img1.Time(1);
  Time = img1.Time(1) : dt : img2.Time(end);
  Depth = Time * c/2;
  img2.Data= interp1(img2.Time,img2.Data,Time,'linear',0);
  
  img_bins(1) = find(Time > param.wf_comb, 1);
  img_bins(2) = img_bins(1) + 10;
  img_bin_comp = img_bins(1) + (30:40);
  
  % Combine waveforms
  difference = mean(mean(img1.Data(img_bin_comp,:))) ...
    ./ mean(mean(img2.Data(img_bin_comp,:)));
  
  trans_bins = img_bins(1)+1:img_bins(2);
  weights = 0.5+0.5*cos(pi*linspace(0,1,length(trans_bins)).');
  data = [img1.Data(1:img_bins(1),:); ...
    repmat(weights,[1 size(img1.Data,2)]).*img1.Data(trans_bins,:) ...
    + difference*repmat(1-weights,[1 size(img2.Data,2)]).*img2.Data(trans_bins,:); ...
    difference*img2.Data(img_bins(2)+1:end,:)];
  % ==============================================================
  
  data_param.time = Time;
  
end

if param.load_sar_data.combine_channels && param.load_sar_data.incoherent ...
    && param.load_sar_data.detrend.cmd ~= 5
  % Detrend data
  detrend = param.load_sar_data.detrend;
  if param.load_sar_data.combine_waveforms
    data = local_detrend(data, detrend.B_noise, ...
      detrend.B_sig, detrend.cmd, detrend.minVal);
  else
    for img = 1:length(param.imgs)
      data{img} = local_detrend(data{img}, detrend.B_noise, ...
        detrend.B_sig, detrend.cmd, detrend.minVal);
    end
  end
end

return;

