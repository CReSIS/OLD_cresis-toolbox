function csarp(param,param_override)
% csarp(param,param_override)
%
% SAR processor function which breaks up the frame into chunks
% which are processed by csarp_task.m
%
% param = struct with processing parameters
%         -- OR --
%         function handle to script with processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Authors: William Blake, John Paden
%
% Example:
%  See run_csarp.m for how to run this function directly.
%  Normally this function is called from master.m using the param spreadsheet.
%
% See also: run_master.m, master.m, run_csarp.m, csarp.m,
%   csarp_task.m

% =====================================================================
% General Setup
% =====================================================================

if ~isstruct(param)
  % Functional form
  param();
end
param = merge_structs(param, param_override);

dbstack_info = dbstack;
fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', dbstack_info(1).name, param.day_seg, datestr(now,'HH:MM:SS'));
fprintf('=====================================================================\n');

% =====================================================================
% Setup processing
% =====================================================================

if ~isfield(param.csarp,'frm_overlap') || isempty(param.csarp.frm_overlap) ...
    || param.csarp.frm_overlap == 0
  param.csarp.frm_overlap = 0;
else
  error('A nonzero frame overlap is no longer allowed. Either remove the field or set to zero.');
end

% Get WGS84 ellipsoid parameters
physical_constants;

csarp_out_path = ct_filename_out(param,param.csarp.out_path,'CSARP_out');
if ~exist(csarp_out_path,'dir')
  mkdir(csarp_out_path);
end

load(ct_filename_support(param,param.records.frames_fn,'frames'));
records = load(ct_filename_support(param,param.records.records_fn,'records'));

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

if ~isfield(param.sched,'rerun_only') || isempty(param.sched.rerun_only)
  param.sched.rerun_only = false;
end

% Cleanup/remove/delete old folders
if ~param.sched.rerun_only
  for frm = param.cmd.frms
    if exist(csarp_out_path,'dir')
        if strcmpi(param.csarp.sar_type,'f-k')
          del_paths = get_filenames(csarp_out_path,sprintf('fk_data_%03d',frm),'','',struct('type','d'));
        elseif strcmpi(param.csarp.sar_type,'tdbp')
          del_paths = get_filenames(csarp_out_path,sprintf('tdpb_data_%03d',frm),'','',struct('type','d'));            
        elseif strcmpi(param.csarp.sar_type,'mltdp')
          del_paths = get_filenames(csarp_out_path,sprintf('mltdp_data_%03d',frm),'','',struct('type','d'));
        else
          error('Invalid SAR processing type (%s)\n', param.csarp.sar_type);
        end
      for idx = 1:length(del_paths)
        fprintf('Removing path: %s\n', del_paths{idx});
        rmdir(del_paths{idx},'s');
      end
    end
  end
end

if ~isfield(param.csarp,'ground_based') || isempty(param.csarp.ground_based)
  param.csarp.ground_based = false;
end

global g_data;
g_data = [];

% =====================================================================
% Collect waveform information into one structure
%  - This is used to break the frame up into chunks
% =====================================================================
if strcmpi(param.radar_name,'mcrds')
  wfs = load_mcrds_wfs(records.settings, param, ...
    records.param_records.records.file.adcs, param.csarp);
elseif any(strcmpi(param.radar_name,{'acords','mcords','mcords2','mcords3','mcords4','mcords5','seaice','accum2'}))
  wfs = load_mcords_wfs(records.settings, param, ...
    records.param_records.records.file.adcs, param.csarp);
elseif any(strcmpi(param.radar_name,{'icards'}))% add icards---qishi
  wfs = load_icards_wfs(records.settings, param, ...
    records.param_records.records.file.adcs, param.csarp);
elseif any(strcmpi(param.radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5'}))
  wfs = load_fmcw_wfs(records.settings, param, ...
    records.param_records.records.file.adcs, param.csarp);
  for wf=1:length(wfs)
    wfs(wf).time = param.csarp.time_of_full_support;
    wfs(wf).freq = 1;
  end
end

% =====================================================================
% Create reference trajectory (used for structuring processing)
% =====================================================================
trajectory_param = struct('rx_path', 0, ...
  'tx_weights', [], 'lever_arm_fh', param.csarp.lever_arm_fh);
% records = trajectory_with_leverarm(records,trajectory_param);
% Lsar = use approximate SAR aperture length
if isfield(param.csarp,'Lsar')
  Lsar = c/wfs(1).fc*(param.csarp.Lsar.agl+param.csarp.Lsar.thick/sqrt(er_ice))/(2*param.csarp.sigma_x);
else
  Lsar = c/wfs(1).fc*(500+1000/sqrt(er_ice))/(2*param.csarp.sigma_x);
end
trajectory_param = struct('gps_source',records.gps_source, ...
  'season_name',param.season_name,'radar_name',param.radar_name,'rx_path', 0, ...
  'tx_weights', [], 'lever_arm_fh', param.csarp.lever_arm_fh);
ref = trajectory_with_leverarm(records,trajectory_param);
along_track = geodetic_to_along_track(ref.lat,ref.lon,ref.elev,Lsar);

% =====================================================================
% Setup static inputs for csarp_task
% =====================================================================

% Do not apply channel equalization during csarp combine unless receivers
% are being combined at this stage (csarp-combined method)
if ~param.csarp.combine_rx
  for wf = 1:length(param.radar.wfs)
    param.radar.wfs(wf).chan_equal_dB(:) = 0;
    param.radar.wfs(wf).chan_equal_deg(:) = 0;
  end
end

task_param = param;

% =====================================================================
% Setup the scheduler
% =====================================================================
if strcmpi(param.sched.type,'custom_torque')
  global ctrl; % Make this global for convenience in debugging
  ctrl = torque_new_batch(param);
  fprintf('Torque batch: %s\n', ctrl.batch_dir);
  torque_compile('csarp_task.m',ctrl.sched.hidden_depend_funs,ctrl.sched.force_compile);
elseif ~strcmpi(param.sched.type,'no scheduler')
  % Initialize submission ctrl structure
  ctrl = [];
  ctrl.cmd = 'init';
  ctrl.sched = param.sched;
  ctrl.fd = {};
%   
%   fd = [get_filenames(param.path,'','','.m',struct('recursive',1)); ...
%   get_filenames(param.path,'','','.mexa64',struct('recursive',1))];
% fd_override = [get_filenames(param.path_override,'','','.m',struct('recursive',1)); ...
%   get_filenames(param.path_override,'','','.mexa64',struct('recursive',1))];
% 
% fd = merge_filelists(fd, fd_override);
%   ctrl.fd = fd;

  
  
  ctrl = create_task(ctrl);
  
  % Prepare submission ctrl structure for queuing jobs
  ctrl.cmd = 'task';
end

% =====================================================================
% SAR process each frame
% =====================================================================
out_recs = {};

% Determine overlap of chunks from the range to furthest target
times    = {wfs.time};
times    = cell2mat(times.');
max_time = min(max(times),param.csarp.time_of_full_support);
if param.csarp.ground_based
  chunk_overlap = (max_time*3e8/2/sqrt(param.csarp.start_eps)*c/wfs.fc)/(2*param.csarp.sigma_x); % m (maximum SAR aperture)
else
  chunk_overlap = (max_time*3e8/2)/(2*param.csarp.sigma_x); % m (maximum SAR aperture)
end

% Check to make sure the beam is not desired to steer too far
%   forward/backward
freqs       = {wfs.freq};
freqs       = cell2mat(freqs.');
v_p         = c/sqrt(er_ice);
max_k       = 2*2*pi*max(freqs)/v_p;
max_ang     = max(abs(asind(gen_kx(along_track)/max_k)));
sar_bw = (c/min(freqs)) / (2*param.csarp.sigma_x) * 180/pi;
max_des_ang = (max(abs(param.csarp.sub_aperture_steering))+0.5)*sar_bw;
% if max_des_ang > max_ang
%   error('Subaperture coefficient steers beam too far forward or backward');
% end

% =====================================================================
% Prepare task inputs to minimize raw file reading. Because some of the
% radars write multiple waveforms and adcs to the same file, it makes
% sense for a job to take care of all the wf,adc pairs from that file.
% Radars affected are mcrds, mcords2, and mcords3.
%  Two modes of operation, determined by combine_rx in param spreadsheet:
%  1. combine_rx = true
%     Combine receivers (really combines all wf_adc pairs for each image
%     in param.csarp.imgs)
%  2. combine_rx = false
%     Do not combine receivers (groups wf/adc pairs by board number).
% =====================================================================
force_one_wf_adc_pair_per_job = true;
if ~param.csarp.combine_rx && any(strcmpi(param.radar_name,{'mcrds','mcords2','mcords3'}))
  % Define adc to board mapping
  if strcmpi(param.radar_name,'mcrds')
    num_adcs_per_board = 8;
    num_boards = 1;
  elseif any(strcmpi(param.radar_name,{'mcords2','mcords3'}))
    num_adcs_per_board = 4;
    num_boards = 4;
  end
  
  % Preallocate image list
  imgs_list = cell(1, num_boards);
  if length(param.csarp.imgs) > 1
    for idx = 1:num_boards
      imgs_list{idx} = cell(1, length(param.csarp.imgs));
    end
  end
  
  for img_idx = 1:length(param.csarp.imgs)
    temp_imgs_list = param.csarp.imgs{img_idx};
    
    % Sort the list and remove duplicates
    [tmp,sorted_idxs] = unique(abs(temp_imgs_list),'rows');
    temp_imgs_list = temp_imgs_list(sorted_idxs,:);
    
    % Split the sorted list into num_boards different individual lists,
    % depending on what boards are used (see above)
    for idx = 1:num_boards
      keep_idxs = abs(temp_imgs_list(:,2))>(idx-1)*num_adcs_per_board & abs(temp_imgs_list(:,2))<(num_adcs_per_board*idx+1);
      if iscell(imgs_list{idx})
        imgs_list{idx}{img_idx} = temp_imgs_list(keep_idxs,:);
      else
        imgs_list{idx} = temp_imgs_list(keep_idxs,:);
      end
    end
  end
  % Remove empty lists from image cell array
  keep_idxs = [];
  for task_idx = 1:length(imgs_list)
    for img_idx = 1:length(param.csarp.imgs)
      if iscell(imgs_list{task_idx})
        if ~isempty(imgs_list{task_idx}{img_idx})
          keep_idxs = cat(2, keep_idxs, img_idx);
        end
      else
        if ~isempty(imgs_list{task_idx})
          keep_idxs = cat(2, keep_idxs, task_idx);
        end
      end
    end
    if iscell(imgs_list{1})
      imgs_list{task_idx} = imgs_list{task_idx}(keep_idxs);
      keep_idxs = [];
    end
  end
  if iscell(imgs_list{1})
    for task_idx = 1:length(imgs_list)
      if ~isempty(imgs_list{task_idx})
        keep_idxs = cat(2, keep_idxs, task_idx);
      end
    end
    imgs_list = imgs_list(keep_idxs);
  else
    imgs_list = imgs_list(keep_idxs);
  end
  
else
  % Don't try to recombine the images... just do it the normal way:
  imgs_list = param.csarp.imgs;
end

% If the number of jobs is going to be too small to fully utilize the cluster,
% then break up the jobs as much as possible (or break the jobs up if
% force_one_wf_adc_pair_per_job is set to true)
if ~param.csarp.combine_rx && ...
    (force_one_wf_adc_pair_per_job || ~strcmpi(param.sched.type,'no scheduler'))
  % THIS IS A HACK... 2*64 should come from ClusterSize information
  %   6 IS FROM TYPICAL NUMBER OF CHUNKS PER FRAME
  if force_one_wf_adc_pair_per_job ...
      || length(param.cmd.frms) * 6 * length(imgs_list) < 2*64
    imgs_list = {};
    for img_idx = 1:length(param.csarp.imgs)
      for wf_adc_idx = 1: size(param.csarp.imgs{img_idx},1)
        imgs_list{end+1} = param.csarp.imgs{img_idx}(wf_adc_idx,:);
      end
    end
  end
end

retry_fields = {};
for frm_idx = 1:length(param.cmd.frms)
  frm = param.cmd.frms(frm_idx);
  if ct_proc_frame(frames.proc_mode(frm),param.csarp.frm_types)
    fprintf('csarp %s_%03i (%i of %i) %s\n', param.day_seg, frm, frm_idx, length(param.cmd.frms), datestr(now,'HH:MM:SS'));
  else
    fprintf('Skipping frame %s_%03i (no process frame)\n', param.day_seg, frm);
    continue;
  end
  
  task_param.proc.frm = frm;
  
  % Current frame goes from the start record specified in the frames file
  % to the record just before the start record of the next frame.  For
  % the last frame, the stop record is just the last record in the segment.
  start_rec = frames.frame_idxs(frm);
  if frm < length(frames.frame_idxs)
    stop_rec = frames.frame_idxs(frm+1)-1;
  else
    stop_rec = length(along_track);
  end
  
  % Each frame is processed independently of all the other frames
  % The output data is at along-track positions [0,sigma_x,2*sigma_x, ...]
  %   - where 0 aligns with the first record of the frame
  %   - last output occurs before the last record of the frame, but
  %     will not in general be aligned with it
  output_along_track = along_track(start_rec) : param.csarp.sigma_x : along_track(stop_rec);
  
  % Add additional output_along_track bins for frame overlap at the start
  % and end of the frame, unless this is the first or last frame
  frame_overlap = param.csarp.sigma_x:param.csarp.sigma_x:param.csarp.frm_overlap/2;
  if frm == 1 && frm == length(frames.frame_idxs)
    % Don't add any frame overlap
  elseif frm == 1
    % Only add frame overlap at the end of the frame
    output_along_track = [output_along_track output_along_track(end)+frame_overlap];
  elseif frm == length(frames.frame_idxs)
    % Old add frame overlap at the start of the frame
    output_along_track = [output_along_track(1)+fliplr(-frame_overlap) output_along_track];
  else
    % Add frame overlap at the start and end of the frame
    output_along_track = [output_along_track(1)+fliplr(-frame_overlap) output_along_track output_along_track(end)+frame_overlap];
  end
  
  % Break this frame up into chunks, each chunk will be processed separately
  output_chunk_idxs = get_equal_alongtrack_spacing_idxs(output_along_track,param.csarp.chunk_len);
  
  % If last chunk is less than half the desired chunk size, combine with the
  % previous chunk
  if output_along_track(end) - output_along_track(output_chunk_idxs(end)) < param.csarp.chunk_len/2 ...
      && length(output_chunk_idxs) > 1
    output_chunk_idxs = output_chunk_idxs(1:end-1);
  end
  
  %% SAR process each chunk of data in the frame
  if strcmpi(param.csarp.sar_type,'mltdp') && 0    % set to 1 for surface fit over whole frame
    task_param.proc.along_track_frm = along_track(start_rec:stop_rec);
    [B,A] = butter(4,0.1);
    
    % Force elevation to be smooth (might be required for refraction)
    task_param.proc.smoothed_elevation = filtfilt(B,A,records.elev(start_rec:stop_rec));
    %smoothed_elevation = records.elev;
    
    % Fit surface to polynomial to force it to be smooth (required for refraction)
    %  - Fit is done with special x-axis to prevent bad conditioning
    smoothed_surface = filtfilt(B,A,records.surface(start_rec:stop_rec));
    sz_poly_order = 21;
    xfit = linspace(-1,1,length(smoothed_surface));
    task_param.proc.smoothed_surface = polyval(polyfit(xfit,smoothed_surface,sz_poly_order),xfit);
    if 0  % check if the fit is good 
      figure(1);plot(records.elev(start_rec:stop_rec)); hold on;plot(smoothed_elevation,'r--');
      figure(2);plot(records.surface(start_rec:stop_rec)*150e6); hold on;plot(smoothed_surface*150e6,'r--');
      figure(3);plot(records.elev(start_rec:stop_rec)-records.surface(start_rec:stop_rec)*150e6);hold on;plot(smoothed_elevation-smoothed_surface*150e6,'r--');
    end
  end
  if ~isfield(param.csarp,'start_chk') | isempty(param.csarp.start_chk)
    start_chk = 1;
  else
    start_chk = param.csarp.start_chk;
  end
  if ~isfield(param.csarp,'end_chk') | isempty(param.csarp.end_chk)  
    end_chk = length(output_chunk_idxs);
  else
    end_chk = param.csarp.end_chk; 
  end
  for chunk_idx = start_chk:end_chk 
    % This chunk will process these along track outputs
    task_param.csarp.chunk_id = chunk_idx;
    start_x = output_along_track(output_chunk_idxs(chunk_idx));
    if chunk_idx < length(output_chunk_idxs)
      % Chunk overlap allows combine_wf_chan to work on neighborhoods of pixels
      % (a chunk overlap of 10 allows a neighborhood of 10 along-track pixels)
      stop_x = output_along_track(output_chunk_idxs(chunk_idx+1)-1 + param.csarp.chunk_overlap);
    else
      stop_x = output_along_track(end);
    end
    
    % These are the records which will be used
    cur_recs = [find(along_track > start_x-chunk_overlap,1) ...
      find(along_track < stop_x+chunk_overlap, 1, 'last')];
    task_param.load.recs = cur_recs;
    
    % Along-track information required to register output properly
    task_param.proc.output_along_track_offset = start_x - along_track(cur_recs(1));
    task_param.proc.output_along_track_Nx = round((stop_x - start_x)/param.csarp.sigma_x) + 1;
    if strcmpi(param.csarp.sar_type,'tdbp') 
      output_along_track_idx_step = round(param.csarp.sigma_x/mean(diff(along_track(cur_recs(1):cur_recs(2)))));
      load_idxs = 1:output_along_track_idx_step:cur_recs(2)-cur_recs(1)+1;
      if chunk_idx < length(output_chunk_idxs)
        task_param.proc.output_along_track_idxs = load_idxs(find(along_track(load_idxs + cur_recs(1)-1)>=start_x...
          & along_track(load_idxs + cur_recs(1)-1)<start_x + param.csarp.chunk_len));
      else
        task_param.proc.output_along_track_idxs = load_idxs(find(along_track(load_idxs + cur_recs(1)-1)>=start_x...
          & along_track(load_idxs + cur_recs(1)-1)<=min(start_x + param.csarp.chunk_len,stop_x)));
      end
    end

    fprintf('  Processing chunk %s_%03d:%d records %d to %d, %.3f to %.3f km (%s)\n', ...
      param.day_seg, frm, chunk_idx, cur_recs(1), cur_recs(2), ...
      ([start_x stop_x]-along_track(start_rec))/1000, datestr(now));
    for imgs_list_idx = 1:length(imgs_list)
      if iscell(imgs_list{imgs_list_idx})
        task_param.load.imgs = imgs_list{imgs_list_idx};
      else
        task_param.load.imgs = imgs_list(imgs_list_idx);
      end
      
      % Create a list of waveform/adc pairs for displaying to the string.
      wf_adc_str = '';
      for img = 1:length(task_param.load.imgs)
        for wf_adc_idx = 1:size(task_param.load.imgs{img},1)
          wf = abs(task_param.load.imgs{img}(wf_adc_idx,1));
          adc = abs(task_param.load.imgs{img}(wf_adc_idx,2));
          if isempty(wf_adc_str)
            wf_adc_str = [wf_adc_str sprintf('wf,adc: %d,%d',wf,adc)];
          else
            wf_adc_str = [wf_adc_str sprintf(' %d,%d',wf,adc)];
          end
        end
      end
      
      sub_band_idx = 1;
      
      % =================================================================
      % Rerun only mode: Test to see if we need to run this task
      if param.sched.rerun_only
        file_exists = true;
        for subap = 1:length(param.csarp.sub_aperture_steering)
          % If we are in rerun only mode AND all the csarp task output files
          % already exists, then we do not run the task
          if strcmpi(param.csarp.sar_type,'f-k')
            out_path = fullfile(ct_filename_out(param, ...
              param.csarp.out_path, 'CSARP_out'), ...
              sprintf('fk_data_%03d_%02d_%02d',frm, ...
              subap, sub_band_idx));
          elseif strcmpi(param.csarp.sar_type,'tdbp')
            out_path = fullfile(ct_filename_out(param, ...
              param.csarp.out_path, 'CSARP_out'), ...
              sprintf('tdbp_data_%03d_%02d_%02d',frm, ...
              subap, sub_band_idx));
          elseif strcmpi(param.csarp.sar_type,'mltdp')
            out_path = fullfile(ct_filename_out(param, ...
              param.csarp.out_path, 'CSARP_out'), ...
              sprintf('mltdp_data_%03d_%02d_%02d',frm, ...
              subap, sub_band_idx));
          end
          for rerun_img_idx = 1:length(task_param.load.imgs)
            % Check for the existence of every file that will be created by this job
            if param.csarp.combine_rx
              % Only one is created when combining and it is named using
              % the first wf/adc pair in the list (even though multiple
              % wf/adc pairs are included in the file).
              wf  = abs(task_param.load.imgs{rerun_img_idx}(1,1));
              adc = abs(task_param.load.imgs{rerun_img_idx}(1,2));
              out_fn = sprintf('wf_%02d_adc_%02d_chk_%03d', abs(wf), abs(adc), chunk_idx);
              out_full_fn = fullfile(out_path,[out_fn '.mat']);
              if ~exist(out_full_fn,'file')
                file_exists = false;
              end
            else
              for wf_adc_idx = size(task_param.load.imgs{rerun_img_idx},1)
                wf  = abs(task_param.load.imgs{rerun_img_idx}(wf_adc_idx,1));
                adc = abs(task_param.load.imgs{rerun_img_idx}(wf_adc_idx,2));
                out_fn = sprintf('wf_%02d_adc_%02d_chk_%03d', abs(wf), abs(adc), chunk_idx);
                out_full_fn = fullfile(out_path,[out_fn '.mat']);
                if ~exist(out_full_fn,'file')
                  file_exists = false;
                end
              end
            end
          end
        end
        if file_exists
          fprintf('  %d/%d: Already exists %s [rerun_only skipping], combine_rx %d (%s)\n', ...
            frm, chunk_idx, wf_adc_str, param.csarp.combine_rx, datestr(now));
          continue;
        end
      end
      task_param.proc.sub_band_idx = sub_band_idx;
      task_param.proc.along_track = along_track(task_param.load.recs(1):task_param.load.recs(end));
      
      % =================================================================
      % Execute tasks/jobs
      fh = @csarp_task;
      arg{1} = task_param;
      
      if strcmp(param.sched.type,'custom_torque')
        create_task_param.conforming = true;
        create_task_param.notes = sprintf('%s %d (%d of %d)/%d of %d %s combine_rx %d', ...
          param.day_seg, frm, frm_idx, length(param.cmd.frms), chunk_idx, length(output_chunk_idxs), wf_adc_str, param.csarp.combine_rx);
        ctrl = torque_create_task(ctrl,fh,1,arg,create_task_param);
        
      elseif ~strcmp(param.sched.type,'no scheduler')
        [ctrl,job_id,task_id] = create_task(ctrl,fh,1,arg);
        fprintf('  %d/%d: %s in job,task %d,%d, combine_rx %d (%s)\n', ...
          frm, chunk_idx, wf_adc_str, job_id, task_id, param.csarp.combine_rx, datestr(now));
        retry_fields{job_id,task_id}.wf_adc_str = wf_adc_str;
        retry_fields{job_id,task_id}.arg = arg;
        retry_fields{job_id,task_id}.frm = frm;
        retry_fields{job_id,task_id}.chunk_idx = chunk_idx;
        if ctrl.error_mask ~= 0 && ctrl.error_mask ~= 2
          % Quit if a bad error occurred
          fprintf('Bad errors occurred, quitting (%s)\n\n', datestr(now));
          ctrl.cmd = 'done';
          ctrl = create_task(ctrl);
          return;
        end
      else
        fprintf('  %s, %d (%d of %d)/%d of %d %s combine_rx %d (%s)\n', param.day_seg, frm, frm_idx, ...
          length(param.cmd.frms), chunk_idx, length(output_chunk_idxs), ...
          wf_adc_str, param.csarp.combine_rx, datestr(now));
        success = fh(arg{1});
      end
    end
  end
end

% =======================================================================
% Wait for jobs to complete if a scheduler was used
% =======================================================================
if strcmpi(param.sched.type,'custom_torque')
  % Wait until all submitted jobs to complete
  ctrl = torque_rerun(ctrl);
  if ~all(ctrl.error_mask == 0)
    if ctrl.sched.stop_on_fail
      torque_cleanup(ctrl);
      error('Not all jobs completed, but out of retries (%s)', datestr(now));
    else
      warning('Not all jobs completed, but out of retries (%s)', datestr(now));
      keyboard;
    end
  else
    fprintf('Jobs completed (%s)\n\n', datestr(now));
  end
  torque_cleanup(ctrl);
  
elseif ~strcmpi(param.sched.type,'no scheduler')
  ctrl.cmd = 'done';
  ctrl = create_task(ctrl);
  if ctrl.error_mask ~= 0 && ctrl.error_mask ~= 2
    % Quit if a bad error occurred
    fprintf('Bad errors occurred, quitting (%s)\n', datestr(now));
    if strcmp(ctrl.sched.type,'torque')
      fprintf('Often on the Torque scheduler, these are not bad errors\n');
      fprintf('because of system instabilities (e.g. file IO failure)\n');
      fprintf('and the task simply needs to be resubmitted. If this is the case,\n');
      fprintf('run "ctrl.error_mask = 2" and then run "dbcont".\n');
      keyboard
      if ctrl.error_mask ~= 0 && ctrl.error_mask ~= 2
        return;
      end
    else
      return
    end
  end
  
  retry = 1;
  while ctrl.error_mask == 2 && retry <= param.sched.max_retries
    fprintf('Tasks failed, retry %d of max %d\n', retry, param.sched.max_retries);
    
    % Bookkeeping (move previous run info to "old_" variables)
    old_ctrl = ctrl;
    old_retry_fields = retry_fields;
    retry_fields = {};
    
    % Initialize submission ctrl structure
    ctrl = [];
    ctrl.cmd = 'init';
    ctrl.sched = param.sched;
    ctrl.fd = {};
    ctrl = create_task(ctrl);
    
    % Prepare submission ctrl structure for queuing jobs
    ctrl.cmd = 'task';
    
    for job_idx = 1:length(old_ctrl.jobs)
      for task_idx = old_ctrl.jobs{job_idx}.error_idxs
        [ctrl,job_id,task_id] = create_task(ctrl,fh,1,old_retry_fields{job_idx,task_idx}.arg);
        fprintf('  %d/%d %s in job,task %d,%d, combine_rx %d (%s)\n', ...
          old_retry_fields{job_idx,task_idx}.frm, old_retry_fields{job_idx,task_idx}.chunk_idx, ...
          old_retry_fields{job_idx,task_idx}.wf_adc_str, ...
          job_id, task_id, param.csarp.combine_rx, datestr(now));
        retry_fields{job_id,task_id} = old_retry_fields{job_idx,task_idx};
      end
    end
    ctrl.cmd = 'done';
    ctrl = create_task(ctrl);
    retry = retry + 1;
    
    if ctrl.error_mask ~= 0 && ctrl.error_mask ~= 2
      % Quit if a bad error occurred
      fprintf('Bad errors occurred, quitting (%s)\n\n', datestr(now));
      if strcmp(ctrl.sched.type,'torque')
        fprintf('Often on the Torque scheduler, these are not bad errors\n');
        fprintf('because of system instabilities (e.g. file IO failure)\n');
        fprintf('and the task simply needs to be resubmitted. If this is the case,\n');
        fprintf('run "ctrl.error_mask = 2" and then run "dbcont".\n');
        keyboard
        if ctrl.error_mask ~= 0 && ctrl.error_mask ~= 2
          return;
        end
      else
        return
      end
    end
  end
  if ctrl.error_mask ~= 0
    fprintf('Not all jobs completed, but out of retries (%s)\n', datestr(now));
    return;
  else
    fprintf('Jobs completed (%s)\n\n', datestr(now));
  end
end

fprintf('Done %s\n', datestr(now));

return;

