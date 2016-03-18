% script compress_echogram
%
% Compresses echogram files (CSARP_*). This involves three steps:
% 1. Elevation compensation
% 2. Fast-time truncation according to param spreadsheet posting settings
% 3. Quantization of 32 bit float to user specified type (e.g. uint8)
%
% Inputs are full resolution qlook files and outputs are saved in posting
% directory.  The parameter file's command worksheet is used to determine
% which segments and frames are compressed (using the generic field).
% The posting worksheet is also used.
%
% See run_compress_echogram.m for how to run.
%
% Author: John Paden

% ====================================================================
% Automated Section
% ====================================================================

fprintf('============================================================\n');
fprintf('Compressing echograms %s\n', echogram_dir);
fprintf('============================================================\n');

physical_constants;

% Read in radar processing parameters spreadsheet
params = read_param_xls(param_fn,'','post');
day_segs = {};
for param_idx = 1:length(params)
  day_segs{end+1} = params(param_idx).day_seg;
end

if ~isempty(compress_type)
  compress_func = str2func(sprintf('@%s', compress_type));
end

for param_idx = 1:length(params)
  param = params(param_idx);
  if ~param.cmd.generic
    continue;
  end
  fprintf('Compress echogram %s (%s)\n', param.day_seg, datestr(now,'HH:MM:SS'));
  
  % Determine if segment is sea or land ice
  if ~isempty(regexpi(param.cmd.mission_names,'Sea'))
    param_post_echo_depth_override = param_post_echo_depth_override_sea;
  else
    param_post_echo_depth_override = param_post_echo_depth_override_land;
  end
  
  % Load the frames file
  frames_fn = ct_filename_support(param,'','frames');
  load(frames_fn);
  
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

  % Compress each frame in the list
  for frm = param.cmd.frms
    % Get the filename
    fn = fullfile(echogram_dir,param.day_seg,sprintf('Data_%s_%03d.mat', ...
      param.day_seg, frm));
    
    if ct_proc_frame(frames.proc_mode(frm),frm_types)
      fprintf(' Input %s (%s)\n', fn, datestr(now));
    else
      fprintf(' Skipping %s (%s)\n', fn, datestr(now));
      continue;
    end
    
    if ~exist(fn,'file')
      warning('File %s does not exist!!!!!', fn);
      continue;
    end
    
    clear minData;
    tmp = load(fn);
    if ~isfield(tmp,'minData')
      orig_size = 4*numel(tmp.Data);
      
      % Elevation compensation
      if elev_comp
        max_elev = max(tmp.Elevation);
        dRange = max_elev - tmp.Elevation;
        dt = tmp.Time(2)-tmp.Time(1);
        dBins = round(dRange / (c/2) / dt);
        zero_pad_len = max(abs(dBins));
        tmp.Data = cat(1,tmp.Data,zeros(zero_pad_len,size(tmp.Data,2)));
        tmp.Time = tmp.Time(1) + (tmp.Time(2)-tmp.Time(1)) * (0:size(tmp.Data,1)-1).';
        for rline = 1:size(tmp.Data,2)
          if dBins(rline) > 0
            tmp.Data(1+dBins(rline):end,rline) = tmp.Data(1:end-dBins(rline),rline);
            tmp.Data(1:dBins(rline),rline) = 0;
          elseif dBins(rline) < 0
            tmp.Data(1:end+dBins(rline),rline) = tmp.Data(1-dBins(rline):end,rline);
            tmp.Data(end+dBins(rline)+1:end,rline) = 0;
          end
          tmp.Elevation(rline) = tmp.Elevation(rline) + dBins(rline)*dt*c/2;
          tmp.Surface(rline) = tmp.Surface(rline) + dBins(rline)*dt;
        end
        tmp.Elevation_Correction = dBins;
      else
        dBins = zeros(1,size(tmp.Data,2));
      end
      
      if truncate_data
        % Truncate data in fast-time
        tmp.Depth = (tmp.Time-median(tmp.Surface))*c/2/sqrt(param.post.echo.er_ice);
        Surface_Depth = (tmp.Surface-median(tmp.Surface))*c/2/sqrt(param.post.echo.er_ice);
        if isempty(param.post.echo.depth)
          depth_rng = [tmp.Depth(1) tmp.Depth(end)];
        elseif ~isempty(param_post_echo_depth_override)
          depth_rng = eval(param_post_echo_depth_override);
        else
          depth_rng = eval(param.post.echo.depth);
        end
        depth_rng = depth_rng + guard_depth_rng;
        if depth_rng(1) > min_depth_rng(1)
          depth_rng(1) = min_depth_rng(1);
        end
        if depth_rng(2) < min_depth_rng(2)
          depth_rng(2) = min_depth_rng(2);
        end
        good_bins = find(tmp.Depth > depth_rng(1) & tmp.Depth < depth_rng(2));
        
        % Truncating data adds four variables
        % 1)Truncate_Bins
        tmp.Truncate_Bins = good_bins;
        tmp.Truncate_Median = NaN*zeros(1,size(tmp.Data,2));
        tmp.Truncate_Mean = NaN*zeros(1,size(tmp.Data,2));
        tmp.Truncate_Std_Dev = NaN*zeros(1,size(tmp.Data,2));
        % 2) Truncate_Median, Truncate_Mean, Truncate_Std_Dev
        %    These will have a NaN for a range line when no valid bins were
        %    truncated for that range line.  Zero padded bins from elevation
        %    compensation are not valid.
        for rline = 1:size(tmp.Data,2)
          if 1+dBins(rline) < good_bins(1)
            tmp.Truncate_Median(rline) = median(tmp.Data(1+dBins(rline):good_bins(1),rline));
            tmp.Truncate_Mean(rline) = mean(tmp.Data(1+dBins(rline):good_bins(1),rline));
            tmp.Truncate_Std_Dev(rline) = std(tmp.Data(1+dBins(rline):good_bins(1),rline));
          end
        end
        tmp.Data = tmp.Data(good_bins,:);
      end
      
      if ~isempty(compress_type)
        % Rescale float32 data to compress_type (uint8 or uint16)
        % Two new variables, minData and maxData, get
        % created in this process and Data becomes compressed
        tmp.Data = 10*log10(tmp.Data);
        tmp.minData = max(min(tmp.Data(isfinite(tmp.Data))), median(tmp.Data(isfinite(tmp.Data))) - 20);
        tmp.Data = tmp.Data - tmp.minData;
        tmp.maxData = min(max(tmp.Data(isfinite(tmp.Data))), median(tmp.Data(isfinite(tmp.Data))) + 60);
        
        tmp.Data = compress_func(round(tmp.Data*double(intmax(compress_type))/tmp.maxData));
      end
      
      if 0
        compressed_size = whos('Data');
        fprintf('Original size %d, Compressed sized %d\n', ...
          orig_size.bytes, compressed_size.bytes);
        figure(1); clf;
        imagesc([],tmp.Depth,double(tmp.Data));
        colormap(1-gray(256));
        pause;
      end
      
      out_fn = fullfile(out_dir,fn(1+length(echogram_dir):end));
      out_fn_dir = fileparts(out_fn);
      if ~exist(out_fn_dir,'dir')
        mkdir(out_fn_dir);
      end
      fprintf('  Output %s\n', out_fn);
      tmp = rmfield(tmp,'Depth');
      save('-v6',out_fn,'-struct','tmp');
    end
  end
end

return;
