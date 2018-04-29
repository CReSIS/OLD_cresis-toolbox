function create_records_mcords2(param, param_override)
% create_records_mcords2(param, param_override)
%
% Function for creating records file for MCORDS2 data. The function is
% usually called from master.m but can also be called from
% run_create_records_mcords2.m.
%
% This function should be run after the GPS file has been created.
% For example, cresis-toobox/gps/missions/make_gps_2009_antarctica_DC8.m
%
% This function's output file is used by all other parts of the processing.
%
% param = struct with processing parameters
%           -- OR --
%         function handle to script which creates a structure with the
%           processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Author: John Paden
%
% See also: run_master.m, master.m, run_create_records_mcords2.m, create_records_mcords2.m,
%   create_records_mcords2_sync.m, check_records.m

%% General Setup
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

%% Prep work
% =====================================================================

% Get WGS84 ellipsoid parameters
physical_constants

if ~isfield(param.records.file,'adc_headers') || isempty(param.records.file.adc_headers)
  param.records.file.adc_headers = param.records.file.adcs;
end

if any(param.records.file_version == [404, 407, 408, 411])
  boards = unique(param.records.file.adc_headers);
  FRAME_SYNC = '1ACFFC1D';
elseif any(param.records.file_version == [402, 403])
  boards = unique(floor((param.records.file.adc_headers-1)/4));
  FRAME_SYNC = 'BADA55E5';
else
  error('Unsupported file version\n');
end

if ~isfield(param.records,'use_ideal_epri') || isempty(param.records.use_ideal_epri)
  param.records.use_ideal_epri = false;
end

if ~isfield(param.records,'gps') || isempty(param.records.gps.en)
  param.records.gps.en = true;
end

if ~isfield(param.records,'presum_bug_fixed') || isempty(param.records.presum_bug_fixed)
  param.records.presum_bug_fixed = false;
end

if ~isfield(param.records,'tmp_fn_uses_adc_folder_name') || isempty(param.records.tmp_fn_uses_adc_folder_name)
  param.records.tmp_fn_uses_adc_folder_name = true;
end

%% Get the files
% =====================================================================
clear wfs board_hdrs;
for board_idx = 1:length(boards)
  board = boards(board_idx);
  
  fprintf('Getting files for board %d (%d of %d) (%s)\n', ...
    board, board_idx, length(boards), datestr(now));
  
  %% Get the list of files to include in this records file
  % =====================================================================
  if any(param.records.file_version == [404, 407, 408, 411])
    [base_dir,adc_folder_name,board_fns{board_idx},file_idxs] = get_segment_file_list(param,board);
  elseif any(param.records.file_version == [402, 403])
    [base_dir,adc_folder_name,board_fns{board_idx},file_idxs] = get_segment_file_list(param,(board+1)*4);
  end
  
  if board_idx == 1
    if strcmp(param.season_name,'2012_Greenland_P3') ...
        && strcmp(param.day_seg,'20120416_02')
      % ERROR CAUSED BY DIGITAL SYSTEM: DDS loaded with 2 waveforms,
      % but NI system set to 3
      hdr_master = basic_load_mcords2(board_fns{board_idx}{file_idxs(1)}, ...
        struct('clk',param.radar.fs,'first_byte',2^26));
      init_EPRI_estimate = 40/12137.55; % 12000 Hz PRF with 7+33=40 presums
      rec_size = median(diff(hdr_master.sync_offsets));
      hdr_master.wfs(1).presums = 6;
      hdr_master.wfs(2).presums = 32;
      hdr_master.wfs(3).presums = -1;
      param.radar.prf = 12137.55;
    else
      [init_EPRI_estimate,hdr_master] = create_records_epri_estimate(param,file_idxs,board_fns{board_idx},adc_folder_name);
      if isfield(param.records,'use_ideal_epri') && ~isempty(param.records.use_ideal_epri) && param.records.use_ideal_epri
        % Count the presums
        num_presum = 0;
        for wf = 1:length(hdr_master.wfs)
          if param.records.presum_bug_fixed
            num_presum = num_presum + hdr_master.wfs(wf).presums;
          else
            num_presum = num_presum + hdr_master.wfs(wf).presums + 1;
          end
        end
        init_EPRI_estimate = num_presum/param.radar.prf;
      end
    end
  
    expected_sync_offset = hdr_master.sync_offsets(1);
    epri_est = hdr_master.epri(1);
  end
  
  board_hdrs{board_idx}.epri = [];
  board_hdrs{board_idx}.seconds = [];
  board_hdrs{board_idx}.fractions = [];
  board_hdrs{board_idx}.file_idx = [];
  board_hdrs{board_idx}.offset = [];
  first_record = true;
  
  % Parse files
  for file_idx = file_idxs
    fn = board_fns{board_idx}{file_idx};
    
    [~,fn_name] = fileparts(fn);
    fprintf('  Parsing file %s (%s)\n', fn_name, datestr(now))
    
    %% Check for temporary files
    if param.records.tmp_fn_uses_adc_folder_name
      tmp_hdr_fn = ct_filename_ct_tmp(rmfield(param,'day_seg'),'','headers', ...
        fullfile(adc_folder_name, [fn_name '.mat']));
    else
      tmp_hdr_fn = ct_filename_ct_tmp(rmfield(param,'day_seg'),'','headers',[fn_name '.mat']);
    end
    if exist(tmp_hdr_fn,'file')
      %% There are temporary files, just load these
      hdr_tmp = load(tmp_hdr_fn);
      
      % Concatenate hdr_tmp fields
      board_hdrs{board_idx}.epri(end+1:end+length(hdr_tmp.epri)) = hdr_tmp.epri;
      board_hdrs{board_idx}.seconds(end+1:end+length(hdr_tmp.seconds)) = hdr_tmp.seconds;
      board_hdrs{board_idx}.fractions(end+1:end+length(hdr_tmp.fraction)) = hdr_tmp.fraction;
      board_hdrs{board_idx}.file_idx(end+1:end+length(hdr_tmp.epri)) = file_idx;
      board_hdrs{board_idx}.offset(end+1:end+length(hdr_tmp.offset)) = hdr_tmp.offset;
      
      if file_idx ~= file_idxs(end)
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
        board_hdrs{board_idx}.epri = board_hdrs{board_idx}.epri(1:end-1);
        board_hdrs{board_idx}.seconds = board_hdrs{board_idx}.seconds(1:end-1);
        board_hdrs{board_idx}.fractions = board_hdrs{board_idx}.fractions(1:end-1);
        board_hdrs{board_idx}.file_idx = board_hdrs{board_idx}.file_idx(1:end-1);
        board_hdrs{board_idx}.offset = board_hdrs{board_idx}.offset(1:end-1);
      end
      
      continue;
    end
    
    %% No Temporary Files: Create records and file numbers
    try
      rec_size = 2*hdr_master.rec_size;
      % =====================================================================
      % Check first record
      % =====================================================================
      if param.records.file_version == 404
        hdr = basic_load_mcords4(fn, struct('clk',param.radar.fs/4));
      elseif param.records.file_version == 403
        hdr = basic_load_mcords3(fn, struct('clk',param.radar.fs));
      elseif param.records.file_version == 402
        hdr = basic_load_mcords2(fn, struct('clk',param.radar.fs));
      else
        error('createrecords:fileversion', ...
          'The file version %d must use create_segment_raw_file_list_v2.m or is not supported.', ...
          param.records.file_version);
      end

      if hdr.sync_offsets(1) ~= expected_sync_offset
        error('create_records:framesync_wrong','Framesync wrong');
      end
      if hdr.epri(1) ~= epri_est
        error('create_records:EPRI_wrong','EPRI offset actual %d estimated %d\n', hdr.epri(1), epri_est);
      end

      num_rec = floor((hdr.file_size - hdr.sync_offsets(1)) / rec_size);
      last_rec_offset = hdr.sync_offsets(1) + (num_rec-1)*rec_size;
      epri_est = hdr.epri(1) + num_rec-1;

      if first_record
        board_hdrs{board_idx}.epri = hdr.epri(1) + (0:num_rec-1);
        board_hdrs{board_idx}.seconds = NaN*zeros(1,num_rec);
        board_hdrs{board_idx}.seconds(1) = hdr.seconds(1);
        board_hdrs{board_idx}.fractions = NaN*zeros(1,num_rec);
        board_hdrs{board_idx}.fractions(1) = hdr.fractions(1);
        board_hdrs{board_idx}.file_idx = file_idx*ones(1,num_rec);
        board_hdrs{board_idx}.offset = hdr.sync_offsets(1) + rec_size*(0:num_rec-1);
      else
        if hdr.sync_offsets(1) == 0
          % File starts with this record (rare)
          start_rec = 0;
        else
          % File starts with a truncated record from the previous file
          start_rec = -1;
        end
        board_hdrs{board_idx}.epri = cat(2,board_hdrs{board_idx}.epri, ...
          hdr.epri + (start_rec:num_rec-1));
        cur_rec = length(board_hdrs{board_idx}.seconds) + 1 - start_rec;
        
        board_hdrs{board_idx}.seconds = cat(2,board_hdrs{board_idx}.seconds, ...
          NaN*zeros(1,num_rec-start_rec));
        board_hdrs{board_idx}.seconds(cur_rec) = hdr.seconds(1);
        
        board_hdrs{board_idx}.fractions = cat(2,board_hdrs{board_idx}.fractions, ...
          NaN*zeros(1,num_rec-start_rec));
        board_hdrs{board_idx}.fractions(cur_rec) = hdr.fractions(1);
        
        board_hdrs{board_idx}.file_idx = cat(2,board_hdrs{board_idx}.file_idx, ...
          file_idx*ones(1,num_rec-start_rec));
        board_hdrs{board_idx}.offset = cat(2,board_hdrs{board_idx}.offset, ...
          hdr.sync_offsets(1) + rec_size*(start_rec:num_rec-1));
      end
      first_record = false;

      % =====================================================================
      % Check if file is consistent, by looking at last complete record in
      % the file
      % =====================================================================
      [fid,msg] = fopen(fn,'r','ieee-be');
      if fid < 1
        fprintf('Could not open file %s\n', fn);
        error(msg);
      end    

      % Seek to last record
      fseek(fid,last_rec_offset,-1);

      % Check last record
      frame_sync = fread(fid,1,'uint32');
      if frame_sync ~= hdr_master.frame_sync
        error('create_records:EOF_framesync_wrong','EOF framesync wrong');
      end
      epri = fread(fid,1,'uint32');
      if epri ~= epri_est
        error('create_records:EOF_EPRI_wrong', ...
          'EOF EPRI offset actual %d estimated %d\n', epri, epri_est);
      end
      if param.records.file_version == 404
        fseek(fid,8,0);
        seconds_BCD = fread(fid,1,'uint32').'; % From NMEA string converted to DCB
        board_hdrs{board_idx}.seconds(end) = ...
          3600*(10*mod(floor(seconds_BCD/2^8),2^4) + mod(floor(seconds_BCD/2^12),2^4)) ...
          + 60*(10*mod(floor(seconds_BCD/2^16),2^4) + mod(floor(seconds_BCD/2^20),2^4)) ...
          + (10*mod(floor(seconds_BCD/2^24),2^4) + mod(floor(seconds_BCD/2^28),2^4));
        board_hdrs{board_idx}.fractions(end) = fread(fid,1,'uint32');
      elseif param.records.file_version == 403
        seconds_BCD = fread(fid,1,'uint32').'; % From NMEA string converted to DCB
        board_hdrs{board_idx}.seconds(end) = ...
          3600*(10*mod(floor(seconds_BCD/2^8),2^4) + mod(floor(seconds_BCD/2^12),2^4)) ...
          + 60*(10*mod(floor(seconds_BCD/2^16),2^4) + mod(floor(seconds_BCD/2^20),2^4)) ...
          + (10*mod(floor(seconds_BCD/2^24),2^4) + mod(floor(seconds_BCD/2^28),2^4));
        board_hdrs{board_idx}.fractions(end) = fread(fid,1,'uint32');
      elseif param.records.file_version == 402
        board_hdrs{board_idx}.seconds(end) = fread(fid,1,'uint32');
        board_hdrs{board_idx}.fractions(end) = fread(fid,1,'uint32');
      end
      
      num_rec_jump = 1 + (hdr.file_size - last_rec_offset - rec_size > 0);
      expected_sync_offset = mod(rec_size - (hdr.file_size - last_rec_offset - rec_size), rec_size);
      epri_est = epri_est + num_rec_jump;

      % Close file
      fclose(fid);
      
    catch ME
      if strcmp(ME.identifier,'createrecords:fileversion')
        rethrow(ME)
      end
      try
        fclose(fid);
      end
      fprintf('File error, slow load required: %s, line %d\n', ME.message, ME.stack(1).line);
      [finfo] = frame_sync_info(fn,struct('cont_mode',false,'sync',FRAME_SYNC));
      
      if ~isempty(strfind(ME.identifier,'EOF'))
        num_records_added = num_rec;
      else
        num_records_added = 0;
      end
      cur_rec = length(board_hdrs{board_idx}.epri) - num_records_added;
      num_records_to_add = finfo.num_rec-num_records_added;
      if num_records_to_add < 0
        % Truncate memory
        board_hdrs{board_idx}.epri = board_hdrs{board_idx}.epri(1:end+num_records_to_add);
        board_hdrs{board_idx}.seconds = board_hdrs{board_idx}.seconds(1:end+num_records_to_add);
        board_hdrs{board_idx}.fractions = board_hdrs{board_idx}.fractions(1:end+num_records_to_add);
        board_hdrs{board_idx}.file_idx = board_hdrs{board_idx}.file_idx(1:end+num_records_to_add);
        board_hdrs{board_idx}.offset = board_hdrs{board_idx}.offset(1:end+num_records_to_add);
      elseif num_records_to_add > 0
        % Preallocate memory
        board_hdrs{board_idx}.epri(end+num_records_to_add) = 0;
        board_hdrs{board_idx}.seconds(end+num_records_to_add) = 0;
        board_hdrs{board_idx}.fractions(end+num_records_to_add) = 0;
        board_hdrs{board_idx}.file_idx(end+num_records_to_add) = 0;
        board_hdrs{board_idx}.offset(end+num_records_to_add) = 0;
      end
      
      fid = fopen(fn,'r','ieee-be');
      if param.records.file_version == 404
        fseek(fid,8,0);
        for rec = 1:finfo.num_rec
          fseek(fid,finfo.syncs(rec)+4,-1);
          board_hdrs{board_idx}.epri(cur_rec+rec) = fread(fid,1,'uint32');
          fseek(fid,8,0);
          seconds_BCD = fread(fid,1,'uint32').'; % From NMEA string converted to DCB
          board_hdrs{board_idx}.seconds(cur_rec+rec) = ...
            3600*(10*mod(floor(seconds_BCD/2^8),2^4) + mod(floor(seconds_BCD/2^12),2^4)) ...
            + 60*(10*mod(floor(seconds_BCD/2^16),2^4) + mod(floor(seconds_BCD/2^20),2^4)) ...
            + (10*mod(floor(seconds_BCD/2^24),2^4) + mod(floor(seconds_BCD/2^28),2^4));
          board_hdrs{board_idx}.fractions(cur_rec+rec) = fread(fid,1,'uint32');
          board_hdrs{board_idx}.file_idx(cur_rec+rec) = file_idx;
          board_hdrs{board_idx}.offset(cur_rec+rec) = finfo.syncs(rec);
        end
      elseif param.records.file_version == 403
        for rec = 1:finfo.num_rec
          fseek(fid,finfo.syncs(rec)+4,-1);
          board_hdrs{board_idx}.epri(cur_rec+rec) = fread(fid,1,'uint32');
          seconds_BCD = fread(fid,1,'uint32').'; % From NMEA string converted to DCB
          board_hdrs{board_idx}.seconds(cur_rec+rec) = ...
            3600*(10*mod(floor(seconds_BCD/2^8),2^4) + mod(floor(seconds_BCD/2^12),2^4)) ...
            + 60*(10*mod(floor(seconds_BCD/2^16),2^4) + mod(floor(seconds_BCD/2^20),2^4)) ...
            + (10*mod(floor(seconds_BCD/2^24),2^4) + mod(floor(seconds_BCD/2^28),2^4));
          board_hdrs{board_idx}.fractions(cur_rec+rec) = fread(fid,1,'uint32');
          board_hdrs{board_idx}.file_idx(cur_rec+rec) = file_idx;
          board_hdrs{board_idx}.offset(cur_rec+rec) = finfo.syncs(rec);
        end
      elseif param.records.file_version == 402
        for rec = 1:finfo.num_rec
          fseek(fid,finfo.syncs(rec)+4,-1);
          board_hdrs{board_idx}.epri(cur_rec+rec) = fread(fid,1,'uint32');
          board_hdrs{board_idx}.seconds(cur_rec+rec) = fread(fid,1,'uint32');
          board_hdrs{board_idx}.fractions(cur_rec+rec) = fread(fid,1,'uint32');
          board_hdrs{board_idx}.file_idx(cur_rec+rec) = file_idx;
          board_hdrs{board_idx}.offset(cur_rec+rec) = finfo.syncs(rec);
        end
      end
      fclose(fid);
      num_rec_jump = 1 + length(finfo.syncs) - finfo.num_rec;
      expected_sync_offset = rec_size - (hdr.file_size - finfo.syncs(end));
      epri_est = board_hdrs{board_idx}.epri(end) + num_rec_jump;
    end
  end
  
end

%% Save workspace in case there is a failure
create_records_save_workspace;

%% Correct time, sync GPS data, and save records
create_records_mcords2_sync;

return;
