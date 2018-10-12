function create_records(param,param_override)
% create_records(param,param_override)
%
% Function for creating records file for all radar systems. The function is
% usually called from master.m but can also be called from
% run_create_records.m.
%
% Corrects jumps in utc time, typically one second jumps. This script
% obtains headers from data indicated in the vectors param spreadsheets.
% After loading the files using basic_load*.m, the records files are saved
% in the support directories for the specific radar in .mat form for
% quicker access.
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
% Authors: John Paden
%
% See also: run_master.m, master.m, run_create_records.m, create_records.m,
%   create_records_sync.m, check_records.m

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Input Checks
% =====================================================================

[output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name);

if ~isfield(param.records.file,'boards') || isempty(param.records.file.boards)
  % Assume a single channel system
  param.records.file.boards = {''};
end

% boards: List of file groupings based on how ADC channels are stored in
%   the files.
if any(param.records.file.version == [1:10 101:102 401 404:409 411])
  % Each channel has its own file
  boards = param.records.file.boards;
elseif any(param.records.file.version == [410])
  % MCRDS: All channels in the same file
  boards = 1;
elseif any(param.records.file.version == [402 403])
  % NI MCRDS: 4 channels per board
  boards = unique(floor((param.records.file.adcs-1)/4));
elseif any(param.records.file.version == [9:10 103 412])
  % RSS: Complicated mapping that must be manually specified
  boards = param.records.file.boards;
else
  error('Unsupported file version\n');
end

if ~isfield(param.records,'epri_jump_threshold') || isempty(param.records.epri_jump_threshold)
  param.records.epri_jump_threshold = 10000;
end

if ~isfield(param.records,'file') || isempty(param.records.file)
  param.records.file = [];
end
if ~isfield(param.records.file,'version') || isempty(param.records.file.version)
  error('The param.records.file.version field must be specified.');
end

if ~isfield(param.records,'gps') || isempty(param.records.gps)
  param.records.gps = [];
end
if ~isfield(param.records.gps,'en') || isempty(param.records.gps.en)
  % Assume that GPS synchronization is enabled
  param.records.gps.en = true;
end

if ~isfield(param.records,'presum_mode') || isempty(param.records.presum_mode)
  if any(param.records.file.version == [402 403 404])
    % 8-channel DDS board had a waveform bug that required an extra waveform
    % to be sent and discarded (not used for presumming)
    param.records.presum_mode = 1;
  elseif any(param.records.file.version == [407 408])
    error('The param.records.presum_mode must be specified. Set to 0 if not using the 8-channel 1 GSPS DDS by Ledford/Leuschen. Set to 1 if it was used.');;
  else
    param.records.presum_mode = 0;
  end
end

if ~isfield(param.records,'use_ideal_epri') || isempty(param.records.use_ideal_epri)
  param.records.use_ideal_epri = false;
end

if ~isfield(param.records.gps,'utc_time_halved') || isempty(param.records.gps.utc_time_halved)
  if param.records.file.version == 1
    error('The param.records.gps.utc_time_halved must be specified. Set to 0 if the UTC time in the data files is half of what it should have been.');
  else
    param.records.gps.utc_time_halved = 0;
  end
end

if ~isfield(param.records,'frames') || isempty(param.records.frames)
  param.records.frames = [];
end

if ~isfield(param.records.frames,'mode') || isempty(param.records.frames.mode)
  param.records.frames.mode = 0;
end

%% Load headers from each board
% =====================================================================
clear board_hdrs;
board_hdrs = {};
records = [];
for board_idx = 1:length(boards)
  board = boards{board_idx};
  
  fprintf('Getting files for board %s (%d of %d) (%s)\n', ...
    board, board_idx, length(boards), datestr(now));
  
  %% Load headers: get files
  % =====================================================================
  [base_dir,adc_folder_name,board_fns{board_idx},file_idxs] = get_segment_file_list(param,board_idx);
  
  % Load header from radar data file
  fprintf('Loading raw files %i to %i\n',file_idxs([1 end]));
  
  % Initialize variables to store header fields
  if any(param.records.file.version == [9 10 103 412])
    % Arena based systems
    board_hdrs{board_idx}.file_size = zeros([0 0],'int32');
    board_hdrs{board_idx}.file_idx = zeros([0 0],'int32');
    board_hdrs{board_idx}.offset = zeros([0 0],'int32');
    board_hdrs{board_idx}.mode_latch = zeros([0 0],'int32');
    board_hdrs{board_idx}.subchannel = zeros([0 0],'int32');
    board_hdrs{board_idx}.pps_cntr_latch = zeros([0 0],'double'); % Time when PPS edge came in
    board_hdrs{board_idx}.pps_ftime_cntr_latch = zeros([0 0],'double'); % Fractional time since edge
    board_hdrs{board_idx}.profile_cntr_latch = zeros([0 0],'double'); % PRI counter
    board_hdrs{board_idx}.rel_time_cntr_latch = zeros([0 0],'double'); % 10 MHz counts counter
    cur_idx = 0;
    
  else
    % NI, Rink, Paden, Leuschen, and Ledford systems
    board_hdrs{board_idx}.seconds = zeros([0 0],'uint32');
    board_hdrs{board_idx}.fraction = zeros([0 0],'uint32');
    board_hdrs{board_idx}.offset = zeros([0 0],'int32');
    board_hdrs{board_idx}.file_idx = zeros([0 0],'uint32');
    % records.relative_rec_num = This variable contains the first record
    % number of each file. After the loop runs there will always be one
    % to many elements (161 files will mean 162 elements in the array)
    % and the last entry is the record number that would have been next
    % so that length(hdr.utc_time_sod) = records.relative_rec_num(end)-1
    records.relative_rec_num{board_idx} = 1;
    
    if param.records.file.version ~= 101
      board_hdrs{board_idx}.epri = zeros([0 0],'uint32');
    end
    if param.records.file.version == 2
      board_hdrs{board_idx}.nyquist_zone = zeros([0 0],'uint8');
      board_hdrs{board_idx}.loopback_mode = zeros([0 0],'uint8');
    end
  end
  
  for file_idx = 1:length(file_idxs)
    % Get the temporary filename from the filename list
    file_num = file_idxs(file_idx);
    [~,fn_name] = fileparts(board_fns{board_idx}{file_num});
    if any(param.records.file.version == [9 10 103 412])
      % Update the filename to point to the packet stripped files
      board_fns{board_idx}{file_num} = ct_filename_ct_tmp(rmfield(param,'day_seg'),'','headers', ...
        fullfile(adc_folder_name, [fn_name '.dat']));
    end
    fn = board_fns{board_idx}{file_num};
    dir_info = dir(fn);
    
    fprintf('  %i/%i %s (%s)\n', ...
      file_idx,length(file_idxs), fn, datestr(now,'HH:MM:SS'));
    
    % Load temporary files
    tmp_hdr_fn = ct_filename_ct_tmp(rmfield(param,'day_seg'),'','headers', ...
      fullfile(adc_folder_name, [fn_name '.mat']));
    hdr_tmp = load(tmp_hdr_fn);
    
    %% Concatenate all the fields together
    %  - Note that all fields from the file should have the same hdr_tmp
    %  length. Error in the file if not.
    if any(param.records.file.version == [9 10 103 412])
      % Arena based systems
      if hdr_tmp.offset > 1e9
        keyboard
      end

      board_hdrs{board_idx}.file_size(cur_idx + (1:length(hdr_tmp.mode_latch))) = dir_info.bytes;
      board_hdrs{board_idx}.file_idx(cur_idx + (1:length(hdr_tmp.mode_latch))) = file_num;
      
      board_hdrs{board_idx}.mode_latch(cur_idx + (1:length(hdr_tmp.mode_latch))) = hdr_tmp.mode_latch;
      board_hdrs{board_idx}.offset(cur_idx + (1:length(hdr_tmp.mode_latch))) = hdr_tmp.offset;
      board_hdrs{board_idx}.subchannel(cur_idx + (1:length(hdr_tmp.mode_latch))) = hdr_tmp.subchannel;
      
      board_hdrs{board_idx}.pps_cntr_latch(cur_idx + (1:length(hdr_tmp.mode_latch))) = hdr_tmp.pps_cntr_latch;
      board_hdrs{board_idx}.pps_ftime_cntr_latch(cur_idx + (1:length(hdr_tmp.mode_latch))) = hdr_tmp.pps_ftime_cntr_latch;
      board_hdrs{board_idx}.profile_cntr_latch(cur_idx + (1:length(hdr_tmp.mode_latch))) = hdr_tmp.profile_cntr_latch;
      board_hdrs{board_idx}.rel_time_cntr_latch(cur_idx + (1:length(hdr_tmp.mode_latch))) = hdr_tmp.rel_time_cntr_latch;
      
      cur_idx = cur_idx + length(hdr_tmp.mode_latch);
      
    else
      % NI, Rink, Paden, Leuschen, and Ledford systems
      board_hdrs{board_idx}.seconds(end+1:end+length(hdr_tmp.seconds)) = hdr_tmp.seconds;
      board_hdrs{board_idx}.fraction(end+1:end+length(hdr_tmp.seconds)) = hdr_tmp.fraction;
      board_hdrs{board_idx}.file_idx(end+1:end+length(hdr_tmp.seconds)) = file_num;
      board_hdrs{board_idx}.offset(end+1:end+length(hdr_tmp.seconds)) = int32(hdr_tmp.offset);
      
      if any(param.records.file.version == [1:8 102 401:404 407:408])
        % Ledford, Rink and NI systems have EPRI field
        board_hdrs{board_idx}.epri(end+1:end+length(hdr_tmp.seconds)) = hdr_tmp.epri;
      end
      
      % Copy the waveform structure
      wfs = hdr_tmp.wfs;
      
      % Create records and file numbers
      records.relative_rec_num{board_idx}(file_idx+1) = length(hdr_tmp.seconds)+records.relative_rec_num{board_idx}(file_idx);
      [fn_dir fn_name fn_ext] = fileparts(fn);
      records.relative_filename{board_idx}{file_idx} = [fn_name fn_ext];
      
      % Handle records that span two files
      if file_idx ~= length(file_idxs)
        % The last record in a file is generally incomplete and continues
        % in the next file. This incomplete record is marked as being
        % in the next file (so we increment file_idx) and we use a negative
        % index to indicate that it actually started in this file where the
        % negative index is relative to the end of this file.
        board_hdrs{board_idx}.file_idx(end) = board_hdrs{board_idx}.file_idx(end) + 1;
        file_size = dir(fn);
        board_hdrs{board_idx}.offset(end) = board_hdrs{board_idx}.offset(end) - file_size.bytes;
      else
        % Drop the last record of the last file since it is generally not a
        % complete record and there is no additional file to load which
        % contains the remainder of the record.
        if any(param.records.file.version == [1:8 102 401:404 407:408])
          board_hdrs{board_idx}.epri = board_hdrs{board_idx}.epri(1:end-1);
        end
        board_hdrs{board_idx}.seconds = board_hdrs{board_idx}.seconds(1:end-1);
        board_hdrs{board_idx}.fraction = board_hdrs{board_idx}.fraction(1:end-1);
        board_hdrs{board_idx}.file_idx = board_hdrs{board_idx}.file_idx(1:end-1);
        board_hdrs{board_idx}.offset = board_hdrs{board_idx}.offset(1:end-1);
      end
    end
  end
  
end

%% Correct EPRI for each board individually
% ======================================================================

if any(param.records.file.version == [9 10 103 412])
  % Arena based systems

  % Load XML configs file
  config_fn = ct_filename_ct_tmp(rmfield(param,'day_seg'),'','headers', ...
    fullfile(param.records.config_fn));
  configs = read_arena_xml(config_fn);
  
  for board_idx = 1:length(boards)
    %% Correct EPRI/Arena: Find the PRIs associated with the EPRI profile
    % =====================================================================
    if size(param.records.data_map{board_idx},2) == 4
      % No Profile Processor Digital System (use mode_latch,subchannel instead)
      % Each row of param.records.data_map{board_idx} = [mode_latch channel wf adc]
      
      % Get the first row (which is always the EPRI row)
      epri_mode = param.records.data_map{board_idx}(1,1);
      epri_subchannel = param.records.data_map{board_idx}(1,2);
      
      mask = board_hdrs{board_idx}.mode_latch == epri_mode & board_hdrs{board_idx}.subchannel == epri_subchannel;
    else
      error('Profile mode not supported.');
      % Profile Processor Digital System
      % Each row of param.records.data_map{board_idx} = [profile wf adc]
      epri_profile = param.records.data_map{board_idx}(1,1);
      
      mask = profile == epri_profile;
    end
    % Data before first PPS reset may be incorrectly time tagged so we do not
    % use it
    first_reset = find(diff(board_hdrs{board_idx}.pps_cntr_latch)>0,1);
    if board_hdrs{board_idx}.pps_ftime_cntr_latch(first_reset) > 10e6
      mask(1:first_reset) = 0;
    end
    
    epri_pris = board_hdrs{board_idx}.profile_cntr_latch(mask);
    
    %% Correct EPRI/Arena: Find EPRI jumps and mask out
    % =====================================================================
    jump_idxs = find( abs(diff(double(epri_pris))/configs.total_presums - 1) > 0.1);
    
    bad_mask = zeros(size(epri_pris));
    for jump_idx = jump_idxs
      jump = (epri_pris(jump_idx+1)-epri_pris(jump_idx))/configs.total_presums - 1;
      fprintf('jump_idx: %d, jump: %d\n', jump_idx, jump);
      fprintf('epri_pris: %d to %d\n', epri_pris(jump_idx), epri_pris(jump_idx+1));
      if jump < -0.1
        fprintf('Negative or zero time jump\n');
        keyboard
        bad_mask(jump_idx+1:end) = epri_pris(jump_idx+1:end) < epri_pris(jump_idx);
      elseif jump > 0.1 && jump < 100
        fprintf('Dropped some records\n');
        %keyboard
      elseif jump > 50000
        fprintf('Record header error or dropped > 50000 records\n');
        keyboard
        epri_pris(jump_idx+1) = epri_pris(jump_idx);
        bad_mask(jump_idx+1) = 1;
      elseif jump > 100
        fprintf('Dropped many records or record header error\n');
        %       keyboard
      end
    end
    
    mask(mask) = ~bad_mask;
    
    %% Correct EPRI/Arena: Update epri_pri_idxs
    % Find the first subrecord in each epri_pri and update epri_pri_idxs
    epri_pri_idxs = find(mask);
    for idx = 1:length(epri_pri_idxs)
      %if ~mod(idx-1,10000)
      %  fprintf('%d\n', idx);
      %end
      epri_pri_idx = epri_pri_idxs(idx);
      pri = board_hdrs{board_idx}.profile_cntr_latch(epri_pri_idx);
      if board_hdrs{board_idx}.profile_cntr_latch(1) == pri
        % Special case for first PRI
        epri_pri_idxs(idx) = 1;
      else
        % For second and later PRIs
        while board_hdrs{board_idx}.profile_cntr_latch(epri_pri_idx) == pri
          epri_pri_idx = epri_pri_idx-1; % Search backwards until we find the previous PRI
        end
        epri_pri_idxs(idx) = epri_pri_idx+1; % The subrecord after this is the first PRI in the current EPRI
      end
    end
    
    %% Correct EPRI/Arena: Ensure all records are valid in each EPRI
    % =====================================================================
    mask = zeros(size(board_hdrs{board_idx}.pps_cntr_latch));
    mask(epri_pri_idxs) = 1;
    for idx = 1:length(epri_pri_idxs)-1
      %if ~mod(idx-1,1000000)
      %  fprintf('%d\n', idx);
      %end
      for pri_idx = epri_pri_idxs(idx):epri_pri_idxs(idx+1)-1
        if board_hdrs{board_idx}.mode_latch(pri_idx) >= size(configs.adc,2) ...
            || board_hdrs{board_idx}.subchannel(pri_idx) >= size(configs.adc,3) ...
            || ~isfield(configs.adc{board_idx,board_hdrs{board_idx}.mode_latch(pri_idx)+1,board_hdrs{board_idx}.subchannel(pri_idx)+1},'rg') ...
            || isempty(configs.adc{board_idx,board_hdrs{board_idx}.mode_latch(pri_idx)+1,board_hdrs{board_idx}.subchannel(pri_idx)+1}.rg)
          fprintf('Bad record %d\n', pri_idx);
          mask(pri_idx) = 0;
          break;
        end
      end
    end
    mask(find(mask,1,'last')) = 0; % Remove last potentially incomplete record
    epri_pri_idxs = find(mask);
    
    %% Correct EPRI/Arena: Fix records split between two files
    % Find records that are split between two files and use a negative
    % offset to indicate this.
    for idx = 2:length(epri_pri_idxs)
      if board_hdrs{board_idx}.file_idx(epri_pri_idxs(idx)) > board_hdrs{board_idx}.file_idx(epri_pri_idxs(idx-1))
        % This record is split between files
        if board_hdrs{board_idx}.offset(epri_pri_idxs(idx)) >= 2^31
          % This split is already handled by arena_packet_strip, but we need
          % to convert uint32 number to negative int32.
          board_hdrs{board_idx}.offset(epri_pri_idxs(idx)) = board_hdrs{board_idx}.offset(epri_pri_idxs(idx)) - 2^32;
        else
          % Convert last record from previous file to negative offset in the
          % new file
          board_hdrs{board_idx}.offset(epri_pri_idxs(idx-1)) = board_hdrs{board_idx}.offset(epri_pri_idxs(idx-1)) - board_hdrs{board_idx}.file_size(epri_pri_idxs(idx-1));
          board_hdrs{board_idx}.file_idx(epri_pri_idxs(idx-1)) = board_hdrs{board_idx}.file_idx(epri_pri_idxs(idx-1)) + 1;
        end
      end
    end
    
    % Store new outputs
    board_hdrs{board_idx}.epri_pri_idxs = epri_pri_idxs;
    board_hdrs{board_idx}.mask = logical(mask);
  end

  %% Correct EPRI/Arena: Mask outputs
  for board_idx = 1:length(boards)
    board_hdrs{board_idx}.offset = board_hdrs{board_idx}.offset(board_hdrs{board_idx}.mask);
    board_hdrs{board_idx}.file_idx = board_hdrs{board_idx}.file_idx(board_hdrs{board_idx}.mask);
    board_hdrs{board_idx}.profile_cntr_latch = board_hdrs{board_idx}.profile_cntr_latch(board_hdrs{board_idx}.mask);
    board_hdrs{board_idx}.rel_time_cntr_latch = board_hdrs{board_idx}.rel_time_cntr_latch(board_hdrs{board_idx}.mask);
    board_hdrs{board_idx}.pps_cntr_latch = board_hdrs{board_idx}.pps_cntr_latch(board_hdrs{board_idx}.mask);
    board_hdrs{board_idx}.pps_ftime_cntr_latch = board_hdrs{board_idx}.pps_ftime_cntr_latch(board_hdrs{board_idx}.mask);
    board_hdrs{board_idx} = rmfield(board_hdrs{board_idx},{'file_size','mode_latch','subchannel','mask'});
  end
  
  %% Correct EPRI/Arena: Populate wfs structure
  for board_idx = 1:length(boards)
    board = boards(board_idx);
    
    for map_idx = 1:size(param.records.data_map{board_idx},1)
      wf = param.records.data_map{board_idx}(map_idx,3);
      % adc = param.records.data_map{board_idx}(map_idx,4);
      mode_latch = param.records.data_map{board_idx}(map_idx,1);
      subchannel = param.records.data_map{board_idx}(map_idx,2);
      
      wfs(wf).num_sam = configs.adc{board_idx,mode_latch+1,subchannel+1}.num_sam;
      wfs(wf).bit_shifts = param.radar.wfs(wf).bit_shifts;
      wfs(wf).t0 = param.radar.wfs(wf).Tadc;
      wfs(wf).presums = configs.adc{board_idx,mode_latch+1,subchannel+1}.presums;
    end
  end

else
  % NI, Rink, Paden, Leuschen, and Ledford systems

end

%% Save workspace in case there is a failure
create_records_save_workspace;

%% Correct time, sync GPS data, and save records
create_records_sync;

%% Create frames
% param.records.frames.mode == 0: Do nothing
if param.records.frames.mode == 1
  frames_fn = ct_filename_support(param,'','frames');
  if ~exist(frames_fn,'file')
    autogenerate_frames(param,param_override);
  end
  create_frames(param,param_override);
elseif param.records.frames.mode == 2
  autogenerate_frames(param,param_override);
end
