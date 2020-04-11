function [data,metadata] = load_sar_data(param)
% [data,metadata] = load_sar_data(param)
%
% Loads and concatenates sar data from a single frame. Currently requires
% all data to be in the same directory. Returns data, position, and header
% information associated with the sar data.
%
% param.load_sar_data: parameter structure controlling how data are loaded.
% See input arguments section for details
%
% Author: John Paden, Logan Smith
%
% See also: run_master.m, master.m, run_array.m, array.m, load_sar_data.m,
% array_proc.m, array_task.m, array_combine_task.m
%
% Also used in: run_load_sar_data.m

%% Check input arguments
% =====================================================================

% debug_level
if ~isfield(param.load_sar_data,'debug_level') || isempty(param.load_sar_data.debug_level)
  param.load_sar_data.debug_level = 1;
end

% combine_channels: sum wf-adc pairs together in each image
if ~isfield(param.load_sar_data,'combine_channels') || isempty(param.load_sar_data.combine_channels)
  param.load_sar_data.combine_channels = false;
end

% combine_imgs: combine images together
if ~isfield(param.load_sar_data,'combine_imgs') || isempty(param.load_sar_data.combine_imgs)
  param.load_sar_data.combine_imgs = false;
end

% chunk: two element vector specifying the start and stop chunk to load
% (minimum value is 1 and maximum value is the number of chunks in the
% frame). If the stop chunk is specified, then the stop chunk is set to the
% number of chunks. The default is to load all chunks which is [1 inf].
if ~isfield(param.load_sar_data,'chunk') || isempty(param.load_sar_data.chunk)
  param.load_sar_data.chunk = [1 inf];
end

% detrend: structure controlling detrending. Arguments are passed into
% local_detrend.m. Detrend is disabled by default (detrend.cmd == 0).
% Detrending only runs if combine_channels == true and incoherent == true.
if ~isfield(param.load_sar_data,'detrend') || isempty(param.load_sar_data.detrend)
  param.load_sar_data.detrend.cmd = [];
end
if ~isfield(param.load_sar_data.detrend,'B_noise') || isempty(param.load_sar_data.detrend.B_noise)
  param.load_sar_data.detrend.B_noise = [100 200];
end
if ~isfield(param.load_sar_data.detrend,'B_sig') || isempty(param.load_sar_data.detrend.B_sig)
  param.load_sar_data.detrend.B_sig = [10 20];
end
if ~isfield(param.load_sar_data.detrend,'cmd') || isempty(param.load_sar_data.detrend.cmd)
  param.load_sar_data.detrend.cmd = 0;
end
if ~isfield(param.load_sar_data.detrend,'minVal') || isempty(param.load_sar_data.detrend.minVal)
  param.load_sar_data.detrend.minVal = -inf;
end

% fn: Data directory. Default is "sar" which loads from "CSARP_sar".
if ~isfield(param.load_sar_data,'fn') || isempty(param.load_sar_data.fn)
  param.load_sar_data.fn = 'sar';
end

% frm: Data frame to load
if ~isfield(param.load_sar_data,'frm') || isempty(param.load_sar_data.frm)
  error('The param.load_sar_data.frm field must be set to an integer indicating which frame to load.');
end

% imgs: cell array of images to load. Each cell array contains a 2 by N
% waveform-adc pair list where N is the number of wf-adc pairs to load.
% Default is {[1 1]} which loads one image and this image is pulled from
% waveform = 1, adc = 1.
if ~isfield(param.load_sar_data,'imgs') || isempty(param.load_sar_data.imgs)
  param.load_sar_data.imgs = {[1 1]};
end

% incoherent: logical scalar, default is false, if true, the data will be
% power detected on load with abs()^2. Incoherent only runs if
% combine_channels == true.
if ~isfield(param.load_sar_data,'incoherent') || isempty(param.load_sar_data.incoherent)
  param.load_sar_data.incoherent = false;
end

% sar_type: 'fk' or 'tdbp'. Default is 'fk'.
if ~isfield(param.load_sar_data,'sar_type') || isempty(param.load_sar_data.sar_type)
  param.load_sar_data.sar_type = 'fk';
end

% subap: list of subapertures to load (array of integers, default is 1)
if ~isfield(param.load_sar_data,'subap') || isempty(param.load_sar_data.subap)
  param.load_sar_data.subap = 1;
end

% The base path for all the data
base_path = ct_filename_out(param,param.load_sar_data.fn,'');

physical_constants;

% Initialize memory for outputs
data     = cell(length(param.load_sar_data.imgs),1);
metadata = [];

%% Load subapertures
for subap_idx = 1:length(param.load_sar_data.subap)
  subap = param.load_sar_data.subap(subap_idx);
  
  in_path = fullfile(base_path, ...
    sprintf('%s_data_%03d_%02d_01/',param.load_sar_data.sar_type,param.load_sar_data.frm, subap));
  
  %% Determine which chunks are available for this subaperture
  img = 1;
  wf_adc_list = param.load_sar_data.imgs{img};
  wf_adc = 1;
  wf = wf_adc_list(wf_adc,1);
  adc = wf_adc_list(wf_adc,2);
  fns = get_filenames(in_path,sprintf('wf_%02.0f_adc_%02.0f_chk_',wf,adc),'','.mat');
  valid_chks = [];
  for fns_idx = 1:length(fns)
    [~,fn_name] = fileparts(fns{fns_idx});
    valid_chks(end+1) = str2double(fn_name(end-2:end));
  end
  param.load_sar_data.chunk(param.load_sar_data.chunk==inf) = max(valid_chks);
  chks_to_load = param.load_sar_data.chunk(1):param.load_sar_data.chunk(end);
  
  % Remove chunks that do not exist from chunks_to_load list
  [valid_chks,keep_idxs] = intersect(chks_to_load, valid_chks);
  if length(valid_chks) ~= length(chks_to_load)
    bad_mask = ones(size(chks_to_load));
    bad_mask(keep_idxs) = 0;
    warning('Nonexistent chunks specified in chks_to_load (e.g. chunk "%g" is invalid), removing these', ...
      chks_to_load(find(bad_mask,1)));
    chks_to_load = valid_chks;
  end
  
  %% Subapertures: Load chunks (blocks) of SAR data
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
      for wf_adc = 1:size(wf_adc_list,1)
        wf = wf_adc_list(wf_adc,1);
        adc = wf_adc_list(wf_adc,2);
        
        sar_fn = fullfile(in_path,sprintf('wf_%02.0f_adc_%02.0f_chk_%03.0f.mat',wf,adc,chunk));
        if param.load_sar_data.debug_level >= 2
          fprintf(' %s (%s)\n', sar_fn, datestr(now,'HH:MM:SS'));
        end
        sar_data = load(sar_fn);

        % Only add in non-overlapping part of SAR image (this is to support
        % legacy SAR data format since current format does not have
        % overlap)
        if img == 1 && wf_adc == 1
          if chunk == param.load_sar_data.chunk(1)
            new_idxs = 1:length(sar_data.fcs.gps_time);
          else
            % 1e-3 added to avoid rounding errors... this is a hack
            new_idxs = find(sar_data.fcs.gps_time > fcs{img}{wf_adc}.gps_time(end)+1e-3);
          end
        end
        
        % Get output image positions (not phase centers)
        if img == 1 && wf_adc == 1 && subap_idx == 1
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
            fcs{img}{wf_adc} = sar_data.fcs;
            
          else
            fcs{img}{wf_adc}.gps_time = cat(2,fcs{img}{wf_adc}.gps_time, ...
              sar_data.fcs.gps_time(new_idxs));
            fcs{img}{wf_adc}.x = cat(2,fcs{img}{wf_adc}.x, ...
              sar_data.fcs.x(:,new_idxs));
            fcs{img}{wf_adc}.y = cat(2,fcs{img}{wf_adc}.y, ...
              sar_data.fcs.y(:,new_idxs));
            fcs{img}{wf_adc}.z = cat(2,fcs{img}{wf_adc}.z, ...
              sar_data.fcs.z(:,new_idxs));
            fcs{img}{wf_adc}.origin = cat(2,fcs{img}{wf_adc}.origin, ...
              sar_data.fcs.origin(:,new_idxs));
            fcs{img}{wf_adc}.pos = cat(2,fcs{img}{wf_adc}.pos, ...
              sar_data.fcs.pos(:,new_idxs));
            fcs{img}{wf_adc}.roll = cat(2,fcs{img}{wf_adc}.roll, ...
              sar_data.fcs.roll(new_idxs));
            fcs{img}{wf_adc}.pitch = cat(2,fcs{img}{wf_adc}.pitch, ...
              sar_data.fcs.pitch(new_idxs));
            fcs{img}{wf_adc}.heading = cat(2,fcs{img}{wf_adc}.heading, ...
              sar_data.fcs.heading(new_idxs));
            fcs{img}{wf_adc}.surface = cat(2,fcs{img}{wf_adc}.surface, ...
              sar_data.fcs.surface(new_idxs));
            fcs{img}{wf_adc}.bottom = cat(2,fcs{img}{wf_adc}.bottom, ...
              sar_data.fcs.bottom(new_idxs));
          end
        end
        
        Nt = size(sar_data.fk_data,1);
        if ~param.load_sar_data.combine_channels
          if subap_idx == 1 && wf_adc == 1
            % Allocate memory in a special way when loading 3D data
            data{img}(size(sar_data.fk_data,1), ...
              size(data{img},2)+length(new_idxs),size(wf_adc_list,1),length(param.load_sar_data.subap)) = single(0);
          end
          data{img}(:,cur_rline + (1:length(new_idxs)),wf_adc,subap_idx) = sar_data.fk_data(:,new_idxs);
        else
          % When combining channels, take the mean of the data as it
          % is loaded in to reduce peak memory consumption.
          if wf_adc == 1
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

% Create lat/lon fields in FCS for convenience
for img = 1:length(fcs)
  for wf_adc = 1:length(fcs{img})
    [fcs{img}{wf_adc}.lat,fcs{img}{wf_adc}.lon,fcs{img}{wf_adc}.elev] ...
      = ecef2geodetic(fcs{img}{wf_adc}.origin(1,:) + sum(fcs{img}{wf_adc}.x.*fcs{img}{wf_adc}.pos), ...
      fcs{img}{wf_adc}.origin(2,:) + sum(fcs{img}{wf_adc}.y.*fcs{img}{wf_adc}.pos), ...
      fcs{img}{wf_adc}.origin(3,:) + sum(fcs{img}{wf_adc}.z.*fcs{img}{wf_adc}.pos), ...
      WGS84.ellipsoid);
    fcs{img}{wf_adc}.lat = fcs{img}{wf_adc}.lat * 180/pi;
    fcs{img}{wf_adc}.lon = fcs{img}{wf_adc}.lon * 180/pi;
  end
end

metadata.fcs = fcs;
metadata.wfs = sar_data.wfs;
metadata.param_records = sar_data.param_records;
metadata.param_sar = sar_data.param_sar;

if param.load_sar_data.combine_channels && param.load_sar_data.combine_imgs
  
end

if param.load_sar_data.combine_channels && param.load_sar_data.incoherent ...
    && param.load_sar_data.detrend.cmd
  % Detrend data
  detrend = param.load_sar_data.detrend;
  if param.load_sar_data.combine_imgs
    data = local_detrend(data, detrend.B_noise, ...
      detrend.B_sig, detrend.cmd, detrend.minVal);
  else
    for img = 1:length(param.imgs)
      data{img} = local_detrend(data{img}, detrend.B_noise, ...
        detrend.B_sig, detrend.cmd, detrend.minVal);
    end
  end
end
