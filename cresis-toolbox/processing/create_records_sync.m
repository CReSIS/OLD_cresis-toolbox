% create_records_sync
%
% Script called from create_records
%
% To run this in debug mode, you need to set the debug setup section in
% this file and then just run this as a script (must be run from ">>"
% prompt and not a debug "K>>" prompt).
%
% Author: John Paden

%% Debug Setup (for running create_records_sync directly)
% =====================================================================
dbstack_info = dbstack;
if ~exist('param','var') || isempty(param) || length(dbstack_info) == 1
  new_param = read_param_xls(ct_filename_param('accum_param_2018_Antarctica_TObas.xls'),'20180817_01');

  fn = ct_filename_ct_tmp(new_param,'','records','workspace');
  fn = [fn '.mat'];
  fprintf('Loading workspace %s (%s)\n', fn, datestr(now));
  if exist(fn,'file')
    load(fn);
  else
    error('Temporary records file does not exist');
  end
  param = new_param;
  clear new_param;
  
  clear('param_override');
  param_override.sched.type = 'no scheduler';
  
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
  param = merge_structs(param, param_override);
end

fprintf('Running %s correction and gps sync (%s)\n', param.day_seg, datestr(now));

%% Input checks
% ======================================================================

if ~isfield(param.records,'manual_time_correct') || isempty(param.records.manual_time_correct)
  param.records.manual_time_correct = 0;
end

% Initialize notes string that will be stored in records file
radar_time_notes = '';
epri_notes = '';
clock_notes = '';

%% Correct GPS time using EPRI
% ======================================================================

%% Align all boards using EPRI
% ======================================================================
if any(param.records.file.version == [9 10 103 412])
  % Arena based systems
  
  %% Align/Arena: Create output EPRI vector
  min_epri = inf;
  max_epri = -inf;
  epri_list = [];
  for board_idx = 1:length(boards)
    % Cluster EPRI values
    
    epri_pri_idxs = board_hdrs{board_idx}.profile_cntr_latch;
    epri_jump_threshold = 10000;
    
    [A,B] = sort(epri_pri_idxs);
    med = median(A);
    [~,med_idx] = min(abs(A-med));
    A = A-med;
    dA = diff(A);
    bad_mask = zeros(size(B));
    bad_start_idx = find(dA(med_idx:length(dA)) > epri_jump_threshold,1);
    if ~isempty(bad_start_idx)
      bad_mask(med_idx+bad_start_idx:end) = true;
    end
    bad_start_idx = find(dA(med_idx-1:-1:1) > epri_jump_threshold,1);
    if ~isempty(bad_start_idx)
      bad_mask(med_idx-bad_start_idx:-1:1) = true;
    end
    back_idxs = 1:length(B);
    back_idxs = back_idxs(B);
    if sum(bad_mask) > 0
      warning('%d of %d records show bad out of range EPRI values.', sum(bad_mask), length(epri_pri_idxs));
    end
    epri_pri_idxs = epri_pri_idxs(back_idxs(logical(~bad_mask)));
    
    % Remove isolated EPRI values
    min_epri = min(min_epri,min(epri_pri_idxs));
    max_epri = max(max_epri,max(epri_pri_idxs));
    epri_list(end+(1:length(epri_pri_idxs))) = epri_pri_idxs;
    diff_epri_pri_idxs = diff(epri_pri_idxs);
    diff_epri(board_idx) = median(diff_epri_pri_idxs);
    min_score = inf;
    for offset = 0:diff_epri(board_idx)
      score = sum(mod((epri_pri_idxs - offset)/diff_epri(board_idx),1) ~= 0);
      if score < min_score
        min_score = score;
      end
    end
    if min_score > 0
      warning('%d of %d records show slipped EPRI values.', min_score, length(epri_pri_idxs));
    end
  end
  master_epri = mode(epri_list);
  if any(diff_epri ~= diff_epri(1))
    error('Inconsistent EPRI step size between boards. Should all be the same: %s', mat2str_generic(diff_epri));
  end
  epri = [fliplr(master_epri:-diff_epri(1):min_epri), master_epri+diff_epri:diff_epri:max_epri];
  
  %% Align/Arena: Fill in missing records from each board
  records.raw.pps_cntr_latch = nan(size(epri));
  records.raw.pps_ftime_cntr_latch = nan(size(epri));
  records.raw.profile_cntr_latch = nan(size(epri));
  records.raw.rel_time_cntr_latch = nan(size(epri));
  for board_idx = 1:length(boards)
    [~,out_idxs] = intersect(epri,board_hdrs{board_idx}.profile_cntr_latch);
    fprintf('Board %d is missing %d of %d records.\n', board_idx, length(epri)-length(out_idxs), length(epri));
    
    % offset: Missing records filled in with -2^31
    offset = zeros(size(epri),'int32');
    offset(:) = -2^31;
    offset(out_idxs) = board_hdrs{board_idx}.offset;
    board_hdrs{board_idx}.offset = offset;
    
    % file_idx: Missing records filled in with NaN
    file_idx = nan(size(epri));
    file_idx(out_idxs) = board_hdrs{board_idx}.file_idx;
    board_hdrs{board_idx}.file_idx = interp_finite(file_idx,[],'nearest');
    
    % Time stamps are assumed to be the same from each board so each board
    % just writes all of its time stamps to the output records fields.
    records.raw.pps_cntr_latch(out_idxs) = board_hdrs{board_idx}.pps_cntr_latch;
    records.raw.pps_ftime_cntr_latch(out_idxs) = board_hdrs{board_idx}.pps_ftime_cntr_latch;
    records.raw.profile_cntr_latch(out_idxs) = board_hdrs{board_idx}.profile_cntr_latch;
    records.raw.rel_time_cntr_latch(out_idxs) = board_hdrs{board_idx}.rel_time_cntr_latch;
  end
  records.raw.pps_cntr_latch = interp_finite(records.raw.pps_cntr_latch);
  records.raw.pps_ftime_cntr_latch = interp_finite(records.raw.pps_ftime_cntr_latch);
  records.raw.profile_cntr_latch = interp_finite(records.raw.profile_cntr_latch);
  records.raw.rel_time_cntr_latch = interp_finite(records.raw.rel_time_cntr_latch);

  radar_time = double(records.raw.pps_cntr_latch) ...
    + double(records.raw.pps_ftime_cntr_latch)/param.records.file.clk;
else
  
end

%% Correct radar time with EPRI
% ===================================================================
if any(param.records.file.version == [9 10 103 412])
  epri_time = epri/param.radar.prf;
  radar_time_error = epri_time - radar_time;
  epri_time = epri_time - median(radar_time_error);
  radar_time_error = radar_time_error - median(radar_time_error);
  
  if 0
    plot(radar_time_error);
    ylim([-1 1])
  end
  
  bad_idxs = find(abs(radar_time_error)>0.2);

  if 0
    bad_idxs
    test_time = radar_time;
    test_time(bad_idxs) = epri_time(bad_idxs);
    figure(1); clf;
    plot(test_time);
    ylims = ylim;
    hold on;
    plot(radar_time)
    ylim(ylims)
  end
  
  radar_time(bad_idxs) = epri_time(bad_idxs);
  
end

%% Correlate GPS with radar data
% ===================================================================
fprintf('Loading GPS data (%s)\n', datestr(now));

if param.records.gps.en
  records = sync_radar_to_gps(param,records,radar_time);
  
else
  records.gps_time = NaN*zeros(size(radar_time));
  records.lat = NaN*zeros(size(radar_time));
  records.lon = NaN*zeros(size(radar_time));
  records.elev = NaN*zeros(size(radar_time));
  records.roll = NaN*zeros(size(radar_time));
  records.pitch = NaN*zeros(size(radar_time));
  records.heading = NaN*zeros(size(radar_time));
  records.gps_source = 'NA';
end

%% Save records files
% =====================================================================

% Save concatenated files in records directories after time fix
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
if param.ct_file_lock
  records.file_version = '1L';
else
  records.file_version = '1';
end

records.settings = [];

% Create the first entry in the records.settings field
records.settings.wfs_records = 1;
records.settings.wfs = wfs;

records.notes = cat(2,sprintf('\nEPRI NOTES\n%s',epri_notes), ...
  sprintf('\nCLOCK NOTES\n%s',clock_notes), ...
  sprintf('\nRADAR TIME NOTES\n%s',radar_time_notes));
records.param_records = param;

fprintf('Saving records file %s (%s)\n',records_fn,datestr(now));
ct_file_lock_check(records_fn,3);
save(records_fn,'-v7.3','-struct','records');

%% Create record aux files for faster loading times
% =====================================================================

fprintf('Creating auxiliary records files %s (%s)\n',records_fn,datestr(now));
create_records_aux_files(records_fn);

