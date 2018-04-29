function create_records_fmcw_accum(param,param_override)
% create_records_fmcw_accum(param,param_override)
%
% Function for creating records file for some versions of the accumulation
% radar, kuband, and snow radars. The function is usually called from
% master.m but can also be called from run_create_records_mcords2.m.
%
% Corrects jumps in utc time, typically one second jumps. This script
% obtains fmcw headers from data indicated in the vectors param
% spreadsheets. After loading the files using basic_load*.m, the records files
% are saved in the support directories for the specific radar in .mat form
% for quicker access.
%
% Output file contains:
% hdr: structure with the following fields
%    .utc_time_sod: vector of corrected utc times
%    .seconds: vector
%    .fraction: vector
%
% Inputs:
% param = struct with processing parameters
%         -- OR --
%         function handle to script with processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Authors: Aric Beaver, John Paden
%
% See also: run_master.m, master.m, run_create_records_fmcw_accum.m, create_records_fmcw_accum.m,
%   create_records_fmcw_accum_sync.m, check_records.m

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
%% Check input parameters
% =====================================================================

[output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name);

if ~isfield(param.records,'force_all') || isempty(param.records.force_all)
  param.records.force_all = false;
end

if ~isfield(param.records,'file_version') || isempty(param.records.file_version)
  param.records.file_version = 1;
end

if ~isfield(param.records,'tmp_fn_uses_adc_folder_name') || isempty(param.records.tmp_fn_uses_adc_folder_name)
  param.records.tmp_fn_uses_adc_folder_name = 1;
end

if ~isfield(param.vectors.gps,'utc_time_halved') || isempty(param.vectors.gps.utc_time_halved)
  param.vectors.gps.utc_time_halved = [];
end

% =====================================================================
%% Setup the scheduler
% =====================================================================
if strcmp(radar_name,'accum')
  fh = @basic_load_accum;
elseif any(strcmp(radar_name,{'snow','kuband'}))
  fh = @basic_load_fmcw;
elseif any(strcmp(radar_name,{'snow2','kuband2'}))
  fh = @basic_load_fmcw2;
elseif any(strcmp(radar_name,{'snow3','kuband3','kaband3'}))
  if param.records.file_version == 4
    fh = @basic_load_fmcw2;
  elseif param.records.file_version == 6
    fh = @basic_load_fmcw4;
  else
    fh = @basic_load_fmcw3;
  end
elseif any(strcmp(radar_name,{'snow5'}))
  fh = @basic_load;
elseif any(strcmp(radar_name,{'snow8'}))
  fh = @basic_load_fmcw8;
end

param.sched.type = 'no scheduler'; % Scheduler no longer used, force to no scheduler
if strcmpi(param.sched.type,'custom_torque')
  global ctrl; % Make this global for convenience in debugging
  ctrl = torque_new_batch(param);
  fprintf('Torque batch: %s\n', ctrl.batch_dir);
  torque_compile([func2str(fh) '.m'],ctrl.sched.hidden_depend_funs,ctrl.sched.force_compile);
end

% =====================================================================
%% Determine which adcs to load
% =====================================================================
if ~isfield(param.records.file,'adc_headers') || isempty(param.records.file.adc_headers)
  param.records.file.adc_headers = param.records.file.adcs;
end
boards = unique(param.records.file.adc_headers);

% =====================================================================
%% Loop through each ADC and load
% =====================================================================

clear wfs board_hdrs;
board_hdrs = {};
records = [];
for board_idx = 1:length(boards)
  board = boards(board_idx);
  
  fprintf('Getting files for board %d (%d of %d) (%s)\n', ...
    board, board_idx, length(boards), datestr(now));
  
  % =====================================================================
  %% Get the list of raw files to load in
  % =====================================================================
  
  [base_dir,adc_folder_name,board_fns{board_idx},file_idxs] = get_segment_file_list(param,board);
  
  % =====================================================================
  %% Load headers from radar data files
  % =====================================================================
  
  % Load header from radar data file
  fprintf('Loading raw files %i to %i\n',file_idxs([1 end]));
  
  board_hdrs{board_idx}.seconds = zeros([0 0],'uint32');
  board_hdrs{board_idx}.fraction = zeros([0 0],'uint32');
  if param.records.file_version ~= 101
    board_hdrs{board_idx}.epri = zeros([0 0],'uint32');
  end
  board_hdrs{board_idx}.offset = zeros([0 0],'uint32');
  board_hdrs{board_idx}.file_idx = zeros([0 0],'uint32');
  if param.records.file_version == 2
    board_hdrs{board_idx}.nyquist_zone = zeros([0 0],'uint8');
    board_hdrs{board_idx}.loopback_mode = zeros([0 0],'uint8');
  end
  % records.relative_rec_num = This variable contains the first record
  % number of each file. After the loop runs there will always be one
  % to many elements (161 files will mean 162 elements in the array)
  % and the last entry is the record number that would have been next
  % so that length(hdr.utc_time_sod) = records.relative_rec_num(end)-1
  records.relative_rec_num{board_idx} = 1;
  
  init_EPRI_estimate = create_records_epri_estimate(param,file_idxs,board_fns{board_idx});
  
  % load('/home/cresis1/tmp/test_flight_settings_file_20133615.mat');
  %fmcw_settings = fmcw_settings * 86400;
  for file_idx = 1:length(file_idxs)
    file_num = file_idxs(file_idx);
    fn = board_fns{board_idx}{file_num};
    
    %% Prepare arguments based on param.radar_name
    first_byte = 0;
    arg{1} = fn;
    arg{2} = struct('clk',param.radar.fs,'utc_time_halved',param.vectors.gps.utc_time_halved, ...
      'first_byte',first_byte, 'file_version', param.records.file_version, ...
      'records',struct('en',1,'epri',init_EPRI_estimate,'force_all',param.records.force_all));
    if strcmp(radar_name,'accum')
      fh = @basic_load_accum;
      finfo = fname_info_accum(fn);
    elseif any(strcmp(radar_name,{'snow','kuband'}))
      fh = @basic_load_fmcw;
      finfo = fname_info_fmcw(fn);
    elseif any(strcmp(radar_name,{'snow2','kuband2'}))
      fh = @basic_load_fmcw2;
      finfo = fname_info_fmcw(fn);
    elseif any(strcmp(radar_name,{'snow3','kuband3','kaband3'}))
      finfo = fname_info_fmcw(fn);
      if file_idx < length(file_idxs)
        next_finfo = fname_info_fmcw(board_fns{board_idx}{file_idxs(file_idx+1)});
      else
        next_finfo.datenum = inf;
      end
      if param.records.file_version == 4
        fh = @basic_load_fmcw2;
      elseif param.records.file_version == 6
        fh = @basic_load_fmcw4;
      else
        fh = @basic_load_fmcw3;
        if 0
          settings_changed_guard_time = 5;
          if any(fmcw_settings > finfo.datenum*86400-settings_changed_guard_time ...
              & fmcw_settings < next_finfo.datenum*86400+settings_changed_guard_time)
            arg{2}.records.force_all = true;
          else
            arg{2}.records.force_all = param.records.force_all || false;
          end
        else
          arg{2}.records.force_all = true;
        end
      end
    elseif any(strcmp(radar_name,{'snow5'}))
      fh = @basic_load;
      finfo = fname_info_fmcw(fn);
    elseif any(strcmp(radar_name,{'snow8'}))
      fh = @basic_load;
      finfo = fname_info_fmcw(fn);
    end
    
    %% Create/run task for each file
    if strcmpi(param.sched.type,'custom_torque')
      create_task_param.conforming = true;
      create_task_param.notes = sprintf('%i/%i %s', ...
        file_idx,length(file_idxs), fn);
      ctrl = torque_create_task(ctrl,fh,2,arg,create_task_param);
    else
      fprintf('  %i/%i %s (%s)\n', ...
        file_idx,length(file_idxs), fn, datestr(now,'HH:MM:SS'));
      
      [~,fn_name] = fileparts(fn);
      
      if param.records.tmp_fn_uses_adc_folder_name
        tmp_hdr_fn = ct_filename_ct_tmp(rmfield(param,'day_seg'),'','headers', ...
          fullfile(adc_folder_name, [fn_name '.mat']));
      else
        tmp_hdr_fn = ct_filename_ct_tmp(rmfield(param,'day_seg'),'','headers',[fn_name '.mat']);
      end
      if exist(tmp_hdr_fn,'file')
        hdr_tmp = load(tmp_hdr_fn);
      elseif any(strcmp(radar_name,{'snow5','snow8'}))
        error('Temporary header file (%s) not found. Have you run run_create_segment_raw_file_list_v2.m?', tmp_hdr_fn);
      else
        [success,hdr_tmp] = fh(arg{1},arg{2});
      end
      % Concatenate hdr_tmp fields
      board_hdrs{board_idx}.seconds(end+1:end+length(hdr_tmp.seconds)) = hdr_tmp.seconds;
      board_hdrs{board_idx}.fraction(end+1:end+length(hdr_tmp.fraction)) = hdr_tmp.fraction;
      board_hdrs{board_idx}.file_idx(end+1:end+length(hdr_tmp.epri)) = file_idx;
      board_hdrs{board_idx}.offset(end+1:end+length(hdr_tmp.offset)) = hdr_tmp.offset;
      
      if length(board_hdrs{board_idx}.seconds) ~= length(board_hdrs{board_idx}.fraction)
        fprintf('Bad file: seconds and fraction fields do not match.\n');
        keyboard
      end
      
      if param.records.file_version ~= 101
        if isfield(hdr_tmp,'epri')
          board_hdrs{board_idx}.epri(end+1:end+length(hdr_tmp.epri)) = hdr_tmp.epri;
        end
        if param.records.file_version == 2
          board_hdrs{board_idx}.nyquist_zone(end+1:end+length(hdr_tmp.epri)) = hdr_tmp.nyquist_zone;
          board_hdrs{board_idx}.loopback_mode(end+1:end+length(hdr_tmp.epri)) = hdr_tmp.loopback_mode;
        end
        board_hdrs{board_idx}.wfs{file_idx} = hdr_tmp.wfs;
      end
      
      %% Create records and file numbers
      records.relative_rec_num{board_idx}(file_idx+1) = length(hdr_tmp.seconds)+records.relative_rec_num{board_idx}(file_idx);
      [fn_dir fn_name fn_ext] = fileparts(fn);
      records.relative_filename{board_idx}{file_idx} = [fn_name fn_ext];
    end
  end
  
  if ~strcmpi(param.sched.type,'no scheduler')
    %% Wait until all submitted jobs complete
    ctrl = torque_rerun(ctrl);
    if ~all(ctrl.error_mask == 0)
      if ctrl.sched.stop_on_fail
        torque_cleanup(ctrl);
        error('Not all jobs completed, but out of retries (%s)', datestr(now,'HH:MM:SS'));
      else
        warning('Not all jobs completed, but out of retries (%s)', datestr(now,'HH:MM:SS'));
        keyboard;
      end
    else
      fprintf('Jobs completed (%s)\n\n', datestr(now,'HH:MM:SS'));
    end
    
    % ======================================================================
    %% Copy all successful task header outputs
    hdrs = {};
    for file_idx = 1:length(ctrl.out)
      hdrs{file_idx} = ctrl.out{file_idx}.argsout{2};
    end
    
    % ======================================================================
    %% Copy all successful task header outputs to final destination
    for file_idx = 1:length(file_idxs)
      file_num = file_idxs(file_idx);
      fn = board_fns{board_idx}{file_num};
      % Concatenate hdr_tmp fields
      board_hdrs{board_idx}.seconds = cat(2,board_hdrs{board_idx}.seconds,hdrs{file_idx}.seconds);
      board_hdrs{board_idx}.fraction = cat(2,board_hdrs{board_idx}.fraction,hdrs{file_idx}.fraction);
      board_hdrs{board_idx}.epri = cat(2,board_hdrs{board_idx}.epri,hdrs{file_idx}.epri);
      board_hdrs{board_idx}.offset = cat(2,board_hdrs{board_idx}.offset,hdrs{file_idx}.offset);
      board_hdrs{board_idx}.file_idx(end+(1:length(hdrs{file_idx}.seconds))) = file_idx;
      if param.records.file_version == 2
        board_hdrs{board_idx}.loopback_mode = cat(2,board_hdrs{board_idx}.loopback_mode,hdrs{file_idx}.loopback_mode);
        board_hdrs{board_idx}.nyquist_zone = cat(2,board_hdrs{board_idx}.nyquist_zone,hdrs{file_idx}.nyquist_zone);
      end
      
      % Create records and file numbers
      records.relative_rec_num{board_idx}(file_idx+1) = length(hdrs{file_idx}.seconds)+records.relative_rec_num{board_idx}(file_idx);
      board_hdrs{board_idx}.wfs{file_idx} = hdrs{file_idx}.wfs(1);
      [fn_dir fn_name fn_ext] = fileparts(fn);
      records.relative_filename{board_idx}{file_idx} = [fn_name fn_ext];
    end
    
    %% Cleanup torque jobs
    torque_cleanup(ctrl);
    
  end
  
end

%% Save workspace in case there is a failure
create_records_save_workspace;

%% Correct time, sync GPS data, and save records
create_records_fmcw_accum_sync;

return;

