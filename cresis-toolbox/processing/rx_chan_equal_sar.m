function rx_chan_equal_sar(param,param_override)
% rx_chan_equal_sar(param,param_override)
%
% RUN THIS FUNCTION FROM "run_rx_chan_equal_sar"
%
% This function is for computing the receiver coefficients
% from f-k processed SAR data. It requires loading SAR data from N receive
% channels obtained using either waveform 1 or 2 (using surface or bottom
% as reference) then use these data to compute the receiver coefficients.
%
% param = struct with processing parameters
%         -- OR --
%         function handle to script with processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Authors: Peng Seng Tan, John Paden
%
% See also: rx_chan_equal_sar_task.m

% =====================================================================
% General Setup
% =====================================================================

dbstack_info = dbstack;
if ~exist('param','var') || isempty(param) || length(dbstack_info) == 1
  % =====================================================================
  % Debug Setup
  % =====================================================================
  param = read_param_xls(ct_filename_param('rds_param_2009_Greenland_TO.xls'),'20090411_01','equal');
  
  clear('param_override');
  param_override.sched.type = 'no scheduler';
  param_override.sched.rerun_only = false;

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

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', dbstack_info(1).name, param.day_seg, datestr(now,'HH:MM:SS'));
fprintf('=====================================================================\n');

% =====================================================================
% Setup processing
% =====================================================================

if ~isfield(param,'debug_level')
  param.debug_level = 1;
end

in_path = ct_filename_out(param, ...
  param.equal.in_path, 'CSARP_out');

out_path = ct_filename_out(param, ...
  param.equal.out_path, 'CSARP_equalsar');

% =====================================================================
% Setup the scheduler
% =====================================================================

if ~strcmp(param.sched.type,'no scheduler')
  global ctrl; % Make this global for convenience in debugging
  ctrl = torque_new_batch(param);
  fprintf('Torque batch: %s\n', ctrl.batch_dir);
  torque_compile('rx_chan_equal_sar_task.m',ctrl.sched.hidden_depend_funs,ctrl.sched.force_compile);
end

% =====================================================================
% Loop through all the frame directories and process the fk
% chunks in those directories
% =====================================================================

frames_fn = ct_filename_support(param,'','frames');
load(frames_fn);
if isempty(param.cmd.frms)
  param.cmd.frms = 1:length(frames.frame_idxs);
end

for frame_idx = 1:length(param.cmd.frms)
  frm = param.cmd.frms(frame_idx);
  
  fprintf('rx_chan_equal_sar frame %s_%03d %s\n', param.day_seg, frm, datestr(now));

  if strcmpi(param.csarp.sar_type,'f-k')
    sar_type = 'fk';
  else
    sar_type = 'tdc';
  end
  
  % Check for existing equalization files and remove
  out_fn_dir = fullfile(ct_filename_out(param,param.equal.out_path,'CSARP_out'), ...
    sprintf('equal_%03d', frm));
  if exist(out_fn_dir,'dir') && ~param.sched.rerun_only
    fprintf('  Cleaning equal directory %s\n', out_fn_dir);
    rmdir(out_fn_dir,'s');
  end
  
  % Get the list of SAR files for one image (we assume all files and images
  % are uniform and SAR processed with the same settings)
  img = 1;
  wf = param.equal.imgs{img}(1,1);
  adc = param.equal.imgs{img}(1,2);
  fns = get_filenames(fullfile(in_path,sprintf('%s_data_%03d_01_01',sar_type,frm)), ...
    sprintf('wf_%02d',wf),sprintf('adc_%02d',adc),'.mat');
  
  % Look at the first SAR processed data file and verify parameters will
  % work
  load(fns{1},'param_csarp');
  if max(param.combine.rline_rng) - min(param.combine.rline_rng) > param_csarp.csarp.chunk_overlap
    error('SAR processing chunks will not align properly, chunk_overlap too small');
  end
  
  num_chunks_per_task = 1;
  for chunk = 1:num_chunks_per_task:length(fns)
    
    if param.sched.rerun_only
      
      % If we are in rerun only mode AND all the combine_wf_chan task output files
      % already exists, then we do not run the task
      file_exists = true;
      for img = 1:length(param.equal.imgs)
        equal_fn = fullfile(out_fn_dir, sprintf('chk_%03d_img_%02d.mat',chunk,img));
        if ~exist(equal_fn,'file')
          file_exists = false;
        end
      end
      if file_exists
        fprintf('  %d already exists [rerun_only skipping] (%s)\n', ...
          chunk, datestr(now));
        continue;
      end
    end
    
    %% To make the SAR processed chunks fit together seamlessly without
    % having to resample, we determine the start range line output for
    % each chunk.
    % chunks: SAR chunks that will be processed by this task,
    %   plus one additional one for calculating rlines(2)
    chunks = chunk + (0:num_chunks_per_task);
    % chunk_Nx: the number of non-overlapping SAR chunk outputs
    chunk_Nx = floor(param_csarp.csarp.chunk_len / param_csarp.csarp.sigma_x);
    % min_offset: the minimum offset into the SAR chunk which array_proc can
    %   output a full support estimate (since the output uses a neighborhood
    %   of points around the pixel in question, the first output line generally
    %   be from the first input line)
    min_offset = -min(param.combine.rline_rng);
    % rlines(1,:): this will be the first range line output by array_proc
    %   for each SAR chunk this task is array processing
    rlines = [];
    rlines(1,:) = 1+mod(min_offset+param.combine.dline-(chunks-1)*chunk_Nx, param.combine.dline);
    rlines(rlines<min_offset) = rlines(rlines<min_offset) + ceil(param.combine.dline/(1+min_offset)) * param.combine.dline;
    rlines(2,1:end-1) = chunk_Nx + rlines(1,2:end) - param.combine.dline;
    rlines = rlines(:,1:end-1);
    
    % Check if this is the last chunk. This last chunk could have variable
    % length and we want to return all of the data from this chunk. To tell
    % combine_task to do this, we set rlines(2) to infinity for
    % this chunk
    if chunk+num_chunks_per_task-1 >= length(fns)
      rlines(2,end) = inf;
    end
    
    %% Get the filetimes that this task will process
    chunk_last = min(chunk+num_chunks_per_task-1, length(fns));
    param.equal.frm = frm;
    param.equal.chunks = chunk:chunk_last;
    param.equal.rlines = rlines;
    
    %% Execute tasks/jobs
    fh = @rx_chan_equal_sar_task;
    arg{1} = param;
    
    if ~strcmp(param.sched.type,'no scheduler')
      create_task_param.conforming = true;
      create_task_param.notes = sprintf('%d/%d', frm, chunk);
      ctrl = torque_create_task(ctrl,fh,1,arg,create_task_param);
      
    else
      success = fh(arg{1});
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
end

% =====================================================================
% Loop through all the array_path directories and combine
% =====================================================================
for frame_idx = 1:length(param.cmd.frms)
  frm = param.cmd.frms(frame_idx);
  
  out_fn_dir = fullfile(ct_filename_out(param,param.equal.out_path,'CSARP_out'), ...
    sprintf('equal_%03d', frm));
  
  for img = 1:length(param.equal.imgs)
    fns = get_filenames(out_fn_dir,'chk_','',sprintf('_img_%02d.mat', img));
    
    num_chunks_per_task = 1;
    roll_est = [];
    chan_equal_est = [];
    GPS_time= [];
    Latitude = [];
    Longitude = [];
    Elevation = [];
    Roll = [];
    Pitch = [];
    Heading = [];
    Surface = [];
    Bottom = [];
    for chunk = 1:num_chunks_per_task:length(fns)
      equal_fn_name = sprintf('chk_%03d_img_%02d.mat',chunk,img);
      equal_fn = fullfile(out_fn_dir, equal_fn_name);
      fprintf('  Combining %s (%s)\n', equal_fn_name, datestr(now));
      equal = load(equal_fn);
      roll_est = [roll_est equal.roll_est];
      chan_equal_est = [chan_equal_est equal.chan_equal_est];
      
      GPS_time = [GPS_time equal.GPS_time];
      Latitude = [Latitude equal.Latitude];
      Longitude = [Longitude equal.Longitude];
      Elevation = [Elevation equal.Elevation];
      Roll = [Roll equal.Roll];
      Pitch = [Pitch equal.Pitch];
      Heading = [Heading equal.Heading];
      Surface = [Surface equal.Surface];
      Bottom = [Bottom equal.Bottom];
      
      %delete(equal_fn);
    end
    param_equal = equal.param_equal;
    param_records = equal.param_records;
    param_csarp = equal.param_csarp;
    
    out_fn = fullfile(ct_filename_out(param,param.equal.out_path,'CSARP_out'), ...
      sprintf('Equal_%s_%03d_%02d.mat', param.day_seg, frm, img));
    fprintf('  Writing output to %s\n', out_fn);
    save(out_fn, 'roll_est', 'chan_equal_est', 'GPS_time', 'Latitude', ...
      'Longitude', 'Elevation', 'Roll', 'Pitch', 'Heading', 'Surface', ...
      'Bottom', 'param_equal', 'param_records', 'param_csarp');
  end
end

return;

img = 1;
param.radar_name = 'rds';
param.season_name = '2008_Greenland_TO';
param.day_seg = '20080717_07';
load(ct_filename_support(param,'','frames'));
for frm = [1:length(frames.frame_idxs)]
equal = load(fullfile(ct_filename_out(param,'','CSARP_out'),sprintf('Equal_%s_%03d_%02d.mat',param.day_seg,frm,img)));
[B,A] = butter(2,0.2);
roll_est_filt = filtfilt(B,A,equal.roll_est);
figure(1); clf;
plot(equal.GPS_time, -roll_est_filt*180/pi)
hold on;
plot(equal.GPS_time, equal.Roll*180/pi,'r');
hold off;
legend('Estimated','Actual');
xlabel('Range lines');
ylabel('Roll (deg)');
title(sprintf('%s %s %s_%03d img %02d',equal.param_equal.radar_name, equal.param_equal.season_name, ...
  equal.param_equal.day_seg,equal.param_equal.equal.frm,img),'interpreter','none');

% Normalize to reference
equal.chan_equal_est = equal.chan_equal_est ...
  ./ repmat(equal.chan_equal_est(equal.param_equal.equal.ref_wf_adc_idx,:),[size(equal.chan_equal_est,1) 1]);

fprintf('chan_equal_deg:\n');
fprintf('%.2f\t', angle(mean(equal.chan_equal_est,2))*180/pi); fprintf('\n');

pow_equal = lp(mean(abs(equal.chan_equal_est).^2,2));
pow_equal' - pow_equal(3);
fprintf('chan_equal_dB:\n');
fprintf('%.2f\t', pow_equal); fprintf('\n');
figure(2); clf;
%h = plot(angle(filtfilt(B,A,equal.chan_equal_est.'))*180/pi);
h = plot(angle(equal.chan_equal_est.')*180/pi);
hold on
plot(equal.Roll*180/pi,'k');
hold off;
legend(h,{'1','2','3','4','5','6'})

pause
end


plot_geotiff('/N/dc2/projects/cresis/GIS_data/greenland/Landsat-7/Greenland_natural.tif',equal.Latitude,equal.Longitude)


plot(angle(equal.chan_equal_est)*180/pi)
plot(angle(equal.chan_equal_est(1,:))*180/pi)

colorbar

