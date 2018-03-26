function create_records_arena(param, param_override)
% create_records_arena(param, param_override)
%
% Function for creating records file for Arena data. The function is
% usually called from master.m but can also be called from
% run_create_records_arena.m.
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

if param.records.file_version == 412
  boards = param.records.file.boards;
else
  error('Unsupported file version\n');
end

if ~isfield(param.records,'use_ideal_epri') || isempty(param.records.use_ideal_epri)
  param.records.use_ideal_epri = false;
end

if ~isfield(param.records,'gps') || isempty(param.records.gps.en)
  param.records.gps.en = true;
end

%% HACK TO BE FIXED LATER:
param.records.clk = 10e6;
param.radar.prf = 20e3;
param.radar.wfs(1).presums = 1;
param.radar.wfs(2).presums = 1;
param.radar.wfs(1).bit_shifts = 1;
param.radar.wfs(2).bit_shifts = 1;
param.radar.wfs(1).t0 = 0;
param.radar.wfs(2).t0 = 0;

total_presums = 0;
for wf = 1:length(param.radar.wfs)
  total_presums = total_presums + param.radar.wfs(wf).presums;
end

EPRI = total_presums / param.radar.prf;

xml_fn = fullfile(param.vectors.file.base_dir,param.vectors.xml_fn);
settings = read_arena_xml(xml_fn);

%% Get the files
% =====================================================================
clear board_hdrs;
for board_idx = 1:length(boards)
  board = boards(board_idx);
  
  fprintf('Getting files for board %d (%d of %d) (%s)\n', ...
    board, board_idx, length(boards), datestr(now));
  
  %% Get the list of files to include in this records file
  % =====================================================================
  [base_dir,adc_folder_name,board_fns{board_idx},file_idxs] = get_segment_file_list(param,board);
  
  %% Load temporary files
  % =====================================================================
  file_size = [];
  file_idx_array = [];
  mode_latch = [];
  offset = [];
  subchannel = [];
  pps_cntr_latch = [];
  pps_ftime_cntr_latch = []; %
  profile_cntr_latch = []; % PRI counter
  rel_time_cntr_latch = []; % 10 MHz counts counter
  cur_idx = 0;
  for file_idx = file_idxs
    fn = board_fns{board_idx}{file_idx};
    
    [~,fn_name] = fileparts(fn);
    dir_info = dir(fn);
    tmp_hdr_fn = ct_filename_ct_tmp(rmfield(param,'day_seg'),'','headers', ...
      fullfile(adc_folder_name, [fn_name '.mat']));
    
    fprintf('  Parsing file %s (%s)\n', tmp_hdr_fn, datestr(now))
    
    hdr_tmp = load(tmp_hdr_fn);
    if hdr_tmp.offset > 1e9
      keyboard
    end
    
    % Concatenate all the fields together
    file_size(cur_idx + (1:length(hdr_tmp.mode_latch))) = dir_info.bytes;
    file_idx_array(cur_idx + (1:length(hdr_tmp.mode_latch))) = file_idx;
    
    mode_latch(cur_idx + (1:length(hdr_tmp.mode_latch))) = hdr_tmp.mode_latch;
    offset(cur_idx + (1:length(hdr_tmp.offset))) = hdr_tmp.offset;
    subchannel(cur_idx + (1:length(hdr_tmp.subchannel))) = hdr_tmp.subchannel;
    
    pps_cntr_latch(cur_idx + (1:length(hdr_tmp.pps_cntr_latch))) = hdr_tmp.pps_cntr_latch;
    pps_ftime_cntr_latch(cur_idx + (1:length(hdr_tmp.pps_ftime_cntr_latch))) = hdr_tmp.pps_ftime_cntr_latch;
    profile_cntr_latch(cur_idx + (1:length(hdr_tmp.profile_cntr_latch))) = hdr_tmp.profile_cntr_latch;
    rel_time_cntr_latch(cur_idx + (1:length(hdr_tmp.rel_time_cntr_latch))) = hdr_tmp.rel_time_cntr_latch;
    
    cur_idx = cur_idx + length(hdr_tmp.mode_latch);
  end
  
  %% HACKS TO BE FIXED LATER:
  subchannel(:) = 0;
  profile = zeros(size(mode_latch));
  
  %% Find the PRIs associated with the EPRI profile
  % =====================================================================
  epri_profile_idx = find(param.records.profiles{board}(:,3) == param.records.epri_profile(board));
  
  if isfinite(param.records.profiles{board}(1,1))
    % No Profile Processor Digital System (use mode,subchannel instead)
    
    epri_mode = param.records.profiles{board}(epri_profile_idx,1);
    epri_subchannel = param.records.profiles{board}(epri_profile_idx,2);
    
    mask = mode_latch == epri_mode & subchannel == epri_subchannel;
  else
    % Profile Processor Digital System
    %  (indicated by mode and subchannel columns in param.records.profiles{board} not being finite)
    mask = profile == param.records.epri_profile{board};
  end
  
  epri_pris = profile_cntr_latch(mask);
  
  jump_idxs = find( abs(diff(double(epri_pris))/total_presums - 1) > 0.1);
  
  bad_mask = zeros(size(epri_pris));
  for jump_idx = jump_idxs
    jump = (epri_pris(jump_idx+1)-epri_pris(jump_idx))/total_presums - 1;
    fprintf('jump_idx: %d, jump: %d\n', jump_idx, jump);
    fprintf('epri_pris: %d %d\n', epri_pris(jump_idx+1), epri_pris(jump_idx));
    if jump < -0.1
      fprintf('Negative or zero time jump\n');
      keyboard
      bad_mask(jump_idx+1:end) = epri_pris(jump_idx+1:end) < epri_pris(jump_idx);
    elseif jump > 0.1 && jump < 100
      fprintf('Dropped some records\n');
      %keyboard
    elseif jump > 50000
      fprintf('Record header error\n');
      keyboard
      epri_pris(jump_idx+1) = epri_pris(jump_idx);
      bad_mask(jump_idx+1) = 1;
    elseif jump > 100
      fprintf('Dropped many records or record header error\n');
%       keyboard
    end
  end
  
  mask(mask) = ~bad_mask;
  
  % Now find the first subrecord in each epri_pri and update epri_pri_idxs
  epri_pri_idxs = find(mask);
  for idx = 1:length(epri_pri_idxs)
    %if ~mod(idx-1,10000)
    %  fprintf('%d\n', idx);
    %end
    epri_pri_idx = epri_pri_idxs(idx);
    pri = profile_cntr_latch(epri_pri_idx);
    if profile_cntr_latch(1) == pri
      % Special case for first PRI
      epri_pri_idxs(idx) = 1;
    else
      % For second and later PRIs
      while profile_cntr_latch(epri_pri_idx) == pri
        epri_pri_idx = epri_pri_idx-1; % Search backwards until we find the previous PRI
      end
      epri_pri_idxs(idx) = epri_pri_idx+1; % The subrecord after this is the first PRI in the current EPRI
    end
  end
  
  % Make sure all records in each epri are valid
  mask = zeros(size(pps_cntr_latch));
  mask(epri_pri_idxs) = 1;
  for idx = 1:length(epri_pri_idxs)-1
    if ~mod(idx-1,1000000)
     fprintf('%d\n', idx);
    end
    for pri_idx = epri_pri_idxs(idx):epri_pri_idxs(idx+1)-1
      if mode_latch(pri_idx) >= size(settings.hdrs,1) ...
          || subchannel(pri_idx) >= size(settings.hdrs,2) ...
          || isempty(settings.hdrs{mode_latch(pri_idx)+1,subchannel(pri_idx)+1})
        fprintf('Bad record %d\n', pri_idx);
        mask(pri_idx) = 0;
        break;
      end
    end
  end
  epri_pri_idxs = find(mask);
  epri_pri_idxs = epri_pri_idxs(1:end-1); % Remove last potentially incomplete record
  mask(epri_pri_idxs) = 1;
  
  % Find records that are split between two files and use a negative
  % offset to indicate this.
  for idx = 2:length(epri_pri_idxs)
    if file_idx_array(epri_pri_idxs(idx)) > file_idx_array(epri_pri_idxs(idx-1))
      % This record is split between files
      if offset(epri_pri_idxs(idx)) >= 2^31
        % This split is already handled by arena_packet_strip, but we need
        % to convert uint32 number to negative int32.
        offset(epri_pri_idxs(idx)) = offset(epri_pri_idxs(idx)) - 2^32;
      else
        % Convert last record from previous file to negative offset in the
        % new file
        offset(epri_pri_idxs(idx-1)) = offset(epri_pri_idxs(idx-1)) - file_size(epri_pri_idxs(idx-1));
        file_idx_array(epri_pri_idxs(idx-1)) = file_idx_array(epri_pri_idxs(idx-1)) + 1;
      end
    end
  end

  %% Store outputs
  board_hdrs{board_idx}.mask = logical(mask);
  board_hdrs{board_idx}.file_idx = file_idx_array;
  board_hdrs{board_idx}.mode_latch = mode_latch;
  board_hdrs{board_idx}.subchannel = subchannel;
  board_hdrs{board_idx}.profile = profile;
  board_hdrs{board_idx}.offset = offset;
  board_hdrs{board_idx}.pps_cntr_latch = pps_cntr_latch;
  board_hdrs{board_idx}.pps_ftime_cntr_latch = pps_ftime_cntr_latch;
  board_hdrs{board_idx}.profile_cntr_latch = profile_cntr_latch;
  board_hdrs{board_idx}.rel_time_cntr_latch = rel_time_cntr_latch;
end

%% Populate wfs structure
for board_idx = 1:length(boards)
  board = boards(board_idx);
  
  for file_idx = file_idxs(1)
    fn = board_fns{board_idx}{file_idx};
    [hdr,data] = basic_load_arena(fn,param);
    
    wf_adc_idx = find(param.records.wf_adc_profiles == param.records.epri_profile(board) ...
      & param.records.wf_adc_boards == board);
    if isempty(wf_adc_idx)
      error('Could not find (profile,board) = (%d,%d).', ...
        param.records.epri_profile(board), board);
    end
    epri_adc = ceil(wf_adc_idx / size(param.records.wf_adc_profiles,1));
    for profile_idx = 1:size(param.records.profiles{board},1)
      wf_adc_idx = find(param.records.wf_adc_profiles == param.records.profiles{board}(profile_idx,3) ...
        & param.records.wf_adc_boards == board);
      if isempty(wf_adc_idx)
        error('Could not find (profile,board) = (%d,%d).', ...
          param.records.profiles{board}(profile_idx,3), board);
      end
      
      wf = 1+mod(wf_adc_idx-1,size(param.records.wf_adc_profiles,1));
      mode = param.records.profiles{board}(profile_idx,1);
      subchannel = param.records.profiles{board}(profile_idx,2);
      
      wfs(wf).num_sam = size(data{subchannel+1,mode+1},1);
      wfs(wf).bit_shifts = param.radar.wfs(wf).bit_shifts;
      wfs(wf).t0 = param.radar.wfs(wf).Tadc;
      wfs(wf).presums = param.radar.wfs(wf).presums;
    end
  end
end

%% SYNCHRONIZE ALL BOARDS, OMIT BAD RECORDS, AND CORRECT HEADER ERRORS
%==========================================================================
% To Do... This section not done yet

% jump_idxs = find(abs((diff(records.gps_time)-EPRI)) / EPRI > 0.1);
% 
% bad_mask = zeros(size(records.gps_time));
% for jump_idx = jump_idxs
%   jump = (records.gps_time(jump_idx+1)-records.gps_time(jump_idx)-EPRI) / EPRI;
%   fprintf('jump_idx: %d, jump: %d\n', jump_idx, jump);
%   fprintf('records.gps_time: %d %d\n', records.gps_time(jump_idx+1), records.gps_time(jump_idx));
%   if jump <= 0
%     fprintf('Negative or zero time jump\n');
%     keyboard
%     bad_mask(jump_idx+1:end) = records.gps_time(jump_idx+1:end) < records.gps_time(jump_idx);
%   elseif jump > 0.1 && jump < 100
%     fprintf('Dropped some records\n');
%     %keyboard
%   elseif jump > 100
%     fprintf('Dropped many records or record header error\n');
%     keyboard
%   elseif jump > 50000
%     fprintf('Record header error\n');
%     keyboard
%     records.gps_time(jump_idx+1) = records.gps_time(jump_idx);
%     bad_mask(jump_idx+1) = 1;
%   end
% end

% Don't need this here, but this is how wf,adc is determined for a
% particular profile.
%   wf_adc_idx = find(param.records.wf_adc == param.records.epri_profile{board});
%
%   epri_wf = 1+mod(wf_adc_idx-1,size(param.records.wf_adc,1));
%   epri_adc = ceil(wf_adc_idx / size(param.records.wf_adc,1));

%% Create Outputs

for board_idx = 1:length(boards)
  board_hdrs{board_idx}.offset = board_hdrs{board_idx}.offset(board_hdrs{board_idx}.mask);
  board_hdrs{board_idx}.file_idx = board_hdrs{board_idx}.file_idx(board_hdrs{board_idx}.mask);
end

for board_idx = 1
  board_hdrs{board_idx}.profile_cntr_latch = board_hdrs{board_idx}.profile_cntr_latch(board_hdrs{board_idx}.mask);
  board_hdrs{board_idx}.rel_time_cntr_latch = board_hdrs{board_idx}.rel_time_cntr_latch(board_hdrs{board_idx}.mask);
  board_hdrs{board_idx}.pps_cntr_latch = board_hdrs{board_idx}.pps_cntr_latch(board_hdrs{board_idx}.mask);
  board_hdrs{board_idx}.pps_ftime_cntr_latch = board_hdrs{board_idx}.pps_ftime_cntr_latch(board_hdrs{board_idx}.mask);
end

records.gps_time = double(board_hdrs{1}.pps_cntr_latch) ...
  + double(board_hdrs{1}.pps_ftime_cntr_latch)/param.records.clk;

records.raw.profile_cntr_latch = board_hdrs{1}.profile_cntr_latch;
records.raw.rel_time_cntr_latch = board_hdrs{1}.rel_time_cntr_latch;
records.raw.pps_cntr_latch = board_hdrs{1}.pps_cntr_latch;
records.raw.pps_ftime_cntr_latch = board_hdrs{1}.pps_ftime_cntr_latch;

utc_time_sod_corrected = records.gps_time - utc_leap_seconds(records.gps_time(1));

%% Correlate GPS with radar data
% ===================================================================
fprintf('Loading GPS data (%s)\n', datestr(now));

if param.records.gps.en
  records = sync_radar_to_gps(param,records,utc_time_sod_corrected);
  
else
  records.lat = NaN*zeros(size(utc_time_sod_corrected));
  records.lon = NaN*zeros(size(utc_time_sod_corrected));
  records.elev = NaN*zeros(size(utc_time_sod_corrected));
  records.gps_time = NaN*zeros(size(utc_time_sod_corrected));
  records.roll = NaN*zeros(size(utc_time_sod_corrected));
  records.pitch = NaN*zeros(size(utc_time_sod_corrected));
  records.heading = NaN*zeros(size(utc_time_sod_corrected));
  records.gps_source = 'NA';
end

%% Save records files
% =====================================================================

records_fn = ct_filename_support(param,'','records');
[records_fn_dir records_fn_name] = fileparts(records_fn);
if ~exist(records_fn_dir,'dir')
  fprintf('Output directory %s does not exist, creating...\n', records_fn_dir);
  mkdir(records_fn_dir);
end

% Standard Fields
% records
%  .lat
%  .lon
%  .elev
%  .roll
%  .pitch
%  .heading
%  .gps_time
%  .gps_source
records.surface = NaN*zeros(size(records.lat));
records.relative_filename = [];
records.relative_rec_num = [];
records.offset = [];
for board_idx = 1:length(board_hdrs)
  unique_file_idxs = unique(board_hdrs{board_idx}.file_idx);
  for file_idx = 1:length(unique_file_idxs)
    records.relative_rec_num{board_idx}(file_idx) = find(board_hdrs{board_idx}.file_idx == unique_file_idxs(file_idx),1);
    [fn_dir fn_name fn_ext] = fileparts(board_fns{board_idx}{unique_file_idxs(file_idx)});
    records.relative_filename{board_idx}{file_idx} = [fn_name fn_ext];
  end
  records.offset(board_idx,:) = board_hdrs{board_idx}.offset;
end
records.radar_name = param.radar_name;
records.ver = 3;
%  records.raw.epri(1 .. Nx)
%  records.raw.seconds(1 .. Nx)
%  records.raw.fraction(1 .. Nx)

records.settings = [];

records.settings.wfs_records = 1;
records.settings.wfs = wfs;

records.notes = '';
records.param_records = param;

fprintf('Saving records file %s (%s)\n',records_fn,datestr(now));
save(records_fn,'-v7.3','-struct','records'); % Handle large file sizes, so use v7.3

%% Create record aux files for faster loading times
% =====================================================================

fprintf('Creating auxiliary records files %s (%s)\n',records_fn,datestr(now));
create_records_aux_files(records_fn);

fprintf('Done (%s)\n\n', datestr(now));
