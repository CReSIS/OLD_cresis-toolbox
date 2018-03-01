% create_records_mcords2_sync
%
% Script called from create_records_mcords2
%
% To run this in debug mode, you need to set the debug setup section in
% this file and then just run this as a script (must be run from ">>"
% prompt and not a debug "K>>" prompt).
%
% Author: John Paden

dbstack_info = dbstack;
if ~exist('param','var') || isempty(param) || length(dbstack_info) == 1
  % =====================================================================
  % Debug Setup
  % =====================================================================
  new_param = read_param_xls(ct_filename_param('rds_param_2016_Greenland_Polar6.xls'),'20160426_05');
  
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

if ~isfield(param.records,'presum_bug_fixed') || isempty(param.records.presum_bug_fixed)
  param.records.presum_bug_fixed = false;
end

% =====================================================================
%% Synchronize Between Boards
% =====================================================================
fprintf('Fix all the EPRIs (%s)\n', datestr(now));
if 1
  % Special hack code to fix bursts of bad data in mcords4: 20140102,
  % 20140103
  
  epri_notes = [];
    
  %% Remove sections of EPRIs that have only sparse good data
  for board_idx = 1:length(board_hdrs)
    diff_epri = diff(board_hdrs{board_idx}.epri);
    diff_epri_good = diff_epri == 1;
    diff_epri_good(end+1) = 0;
    diff_epri_good(2:end) = diff_epri_good(2:end) | diff_epri == 1;
    
    board_hdrs{board_idx}.epri = board_hdrs{board_idx}.epri(diff_epri_good);
    board_hdrs{board_idx}.seconds = board_hdrs{board_idx}.seconds(diff_epri_good);
    board_hdrs{board_idx}.fractions  = board_hdrs{board_idx}.fractions(diff_epri_good);
    board_hdrs{board_idx}.file_idx = board_hdrs{board_idx}.file_idx(diff_epri_good);
    board_hdrs{board_idx}.offset = board_hdrs{board_idx}.offset(diff_epri_good);
    epri_notes = cat(2,epri_notes,sprintf('Board idx %d: %d bad epri\n', board_idx, sum(~diff_epri_good)));
    
    diff_epri = diff(board_hdrs{board_idx}.epri);
    
    diff_epri_bad = diff_epri > 1000;
    diff_epri_bad_idxs = find(diff_epri_bad)
    for idx = 1:length(diff_epri_bad_idxs)
      diff_epri_bad_idx = diff_epri_bad_idxs(idx);
      epri1 = board_hdrs{board_idx}.epri(diff_epri_bad_idx);
      epri2 = board_hdrs{board_idx}.epri(diff_epri_bad_idx+1);
      fidx1 = board_hdrs{board_idx}.file_idx(diff_epri_bad_idx);
      fidx2 = board_hdrs{board_idx}.file_idx(diff_epri_bad_idx+1);
      fn1 = board_fns{board_idx}{fidx1};
      fn2 = board_fns{board_idx}{fidx2};
      epri_note = sprintf('Board idx %d: EPRI %d-%d: %d records dropped\n  %s (%d)\n  %s (%d)\n', ...
        board_idx, epri1, epri2, diff_epri(diff_epri_bad_idx), fn1, fidx1, fn2, fidx2);
      fprintf(epri_note);
      epri_notes = cat(2,epri_notes,epri_note);
    end
    if ~isempty(diff_epri_bad_idxs)
      warning('Large numbers of EPRIs dropped (dbcont to ignore). Recommend creating new segments for drops greater than 1 km (dbquit).');
%       keyboard
    end
  end
  
  %% Find the first and last EPRI that is common for all channels
  A = {};
  for board_idx = 1:length(board_hdrs)
    A{board_idx} = board_hdrs{board_idx}.epri;
  end
  [min_idx max_idx] = find_common_min_max_elements(A);
  start_epri = A{1}(min_idx);
  stop_epri = A{1}(max_idx);
  
  if stop_epri-start_epri > 100e6
    warning('Too many epri...\n');
    keyboard
  end
  epris = start_epri:stop_epri;

  %% Only use records (epris) where every board contains the record (epri)
  good_mask = zeros(size(epris));
  for board_idx = 1:length(board_hdrs)
    % Rename variable for convenience
    raw_epris = board_hdrs{board_idx}.epri;
    if strcmpi(param.season_name,'2014_Antarctica_DC8') & (strcmpi(param.day_seg(1:8),'20141018') ... 
        | strcmpi(param.day_seg(1:8),'20141023'))
      bad_idxs0 = find(diff(raw_epris,2) ~= 0);
      bad_idxs1 = bad_idxs0 - 5;
      bad_idxs2 = bad_idxs0 + 5;
      bad_idxs = union(union(bad_idxs0,bad_idxs1),bad_idxs2);
      raw_epris = raw_epris(setdiff([1:length(raw_epris)],bad_idxs));  % cleared most of the sparse digital errors, but not all for those in big chunks
    end
      % Add 1 to good_mask for each epri from this board
    good_mask(raw_epris(raw_epris >= start_epri & raw_epris <=stop_epri) - start_epri + 1) ...
      = good_mask(raw_epris(raw_epris >= start_epri & raw_epris <= stop_epri) - start_epri + 1) + 1;
  end
  % An entry in good_mask equal to the number of boards means that every
  % board has that epri.
  good_mask = good_mask == length(board_hdrs);
  epris = epris(good_mask);
  
  % Print out some debug information and keep track in epri_notes
  for board_idx = 1:length(board_hdrs)
    [~,good_idxs] = intersect(board_hdrs{board_idx}.epri,epris);
    dropped_epris = length(board_hdrs{board_idx}.epri) - length(good_idxs);
    epri_notes = cat(2,epri_notes,sprintf('Board idx %d: %d uncommon epri\n', board_idx, dropped_epris));
    board_hdrs{board_idx}.epri = board_hdrs{board_idx}.epri(good_idxs);
    board_hdrs{board_idx}.seconds = board_hdrs{board_idx}.seconds(good_idxs);
    board_hdrs{board_idx}.fractions  = board_hdrs{board_idx}.fractions(good_idxs);
    board_hdrs{board_idx}.file_idx = board_hdrs{board_idx}.file_idx(good_idxs);
    board_hdrs{board_idx}.offset = board_hdrs{board_idx}.offset(good_idxs);
  end
  
else
  
  for board_idx = 1:length(board_hdrs)
  B = board_hdrs{board_idx}.epri;
  
  done = false;
  fprintf('Fixing %d\n', board_idx);
  iter=1;
  while ~done
    fprintf('  Iteration %d\n', iter);
    A = medfilt1(B,11);
    A(1:5) = median(B(1:6));
    A(end-4:end) = median(B(end-5:end));
    bad_mask = abs(A-B) > 50;
    if bad_mask(1)
      first_good_idx = find(~bad_mask,1);
      bad_mask(1:first_good_idx-1) = 0;
      B(1:first_good_idx-1) = B(first_good_idx) + (-first_good_idx+1:-1);
    end
    bad_idxs = find(bad_mask);
    good_idxs = find(bad_mask)-1;
    B(bad_idxs) = B(good_idxs) + 1;
    if isempty(bad_idxs)
      done = true;
    end
    iter = iter+1;
  end
  
  board_hdrs{board_idx}.epri = B;
  end
  
  %% Find the first and last EPRI that is common for all channels
  A = {};
  for board_idx = 1:length(board_hdrs)
  A{board_idx} = board_hdrs{board_idx}.epri;
  end
  [min_idx max_idx] = find_common_min_max_elements(A);
  start_epri = A{1}(min_idx);
  stop_epri = A{1}(max_idx);
  
  epri_notes = [];
  for board_idx = 1:length(board_hdrs)
  fprintf('Board idx %d\n', board_idx);
  epri_notes = cat(2,epri_notes,sprintf('Board idx %d\n', board_idx));
  good_records = find(board_hdrs{board_idx}.epri >= start_epri ...
    & board_hdrs{board_idx}.epri <= stop_epri);
    
  num_trim = length(board_hdrs{board_idx}.epri) - length(good_records);
  if num_trim > 0
    fprintf('  Trimming %d records\n', num_trim);
    epri_notes = cat(2,epri_notes,sprintf('  Trimming %d records\n', num_trim));
    board_hdrs{board_idx}.epri = board_hdrs{board_idx}.epri(good_records);
    board_hdrs{board_idx}.seconds = board_hdrs{board_idx}.seconds(good_records);
    board_hdrs{board_idx}.fractions = board_hdrs{board_idx}.fractions(good_records);
    board_hdrs{board_idx}.file_idx = board_hdrs{board_idx}.file_idx(good_records);
    board_hdrs{board_idx}.offset = board_hdrs{board_idx}.offset(good_records);
  end
  
  cur_epri = board_hdrs{board_idx}.epri(1);
  for idx = 2:length(board_hdrs{board_idx}.epri)
    cur_epri = cur_epri + 1;
    if board_hdrs{board_idx}.epri(idx) < cur_epri
      num_repeated= board_hdrs{board_idx}.epri(idx) - cur_epri;
      fprintf('  Repeated %d records at %d (code tries to fix unless over 3000)\n', -num_repeated, cur_epri);
      epri_notes = cat(2,epri_notes,sprintf('Repeated %d records at %d\n', -num_repeated, cur_epri));
      if num_repeated <= -3000
        fprintf('  This segment should be broken up at epri %d\n', cur_epri);
        fprintf('  This is file idx %d (%s)\n', ...
          board_hdrs{board_idx}.file_idx(board_idx), ...
          board_fns{board_idx}{board_hdrs{board_idx}.file_idx(board_idx)});
        keyboard
      end
      board_hdrs{board_idx}.epri = cat(2, ...
        board_hdrs{board_idx}.epri(1:idx-1), ...
        board_hdrs{board_idx}.epri(idx-num_repeated:end));
      board_hdrs{board_idx}.seconds = cat(2, ...
        board_hdrs{board_idx}.seconds(1:idx-1), ...
        board_hdrs{board_idx}.seconds(idx-num_repeated:end));
      board_hdrs{board_idx}.fractions = cat(2, ...
        board_hdrs{board_idx}.fractions(1:idx-1), ...
        board_hdrs{board_idx}.fractions(idx-num_repeated:end));
      board_hdrs{board_idx}.file_idx = cat(2, ...
        board_hdrs{board_idx}.file_idx(1:idx-1), ...
        board_hdrs{board_idx}.file_idx(idx-num_repeated:end));
      board_hdrs{board_idx}.offset = cat(2, ...
        board_hdrs{board_idx}.offset(1:idx-1), ...
        board_hdrs{board_idx}.offset(idx-num_repeated:end));
    elseif board_hdrs{board_idx}.epri(idx) > cur_epri
      num_skipped = board_hdrs{board_idx}.epri(idx) - cur_epri;
      fprintf('  Skipped %d records at %d (code tries to fix unless over 3000)\n', num_skipped, cur_epri);
      epri_notes = cat(2,epri_notes,sprintf('Skipped %d records at %d\n', num_skipped, cur_epri));
      if num_skipped >= 3000
        fprintf('  This segment should be broken up at epri %d\n', cur_epri);
        fprintf('  This is file idx %d (%s)\n', ...
          board_hdrs{board_idx}.file_idx(idx), ...
          board_fns{board_idx}{board_hdrs{board_idx}.file_idx(idx)});
        keyboard
      end
      board_hdrs{board_idx}.epri = cat(2, ...
        board_hdrs{board_idx}.epri(1:idx-1), ...
        cur_epri:board_hdrs{board_idx}.epri(idx)-1, ...
        board_hdrs{board_idx}.epri(idx:end));
      board_hdrs{board_idx}.seconds = cat(2, ...
        board_hdrs{board_idx}.seconds(1:idx-1), ...
        NaN*ones(1,num_skipped), ...
        board_hdrs{board_idx}.seconds(idx:end));
      board_hdrs{board_idx}.fractions = cat(2, ...
        board_hdrs{board_idx}.fractions(1:idx-1), ...
        NaN*ones(1,num_skipped), ...
        board_hdrs{board_idx}.fractions(idx:end));
      board_hdrs{board_idx}.file_idx = cat(2, ...
        board_hdrs{board_idx}.file_idx(1:idx-1), ...
        board_hdrs{board_idx}.file_idx(idx-1) * ones(1,num_skipped), ...
        board_hdrs{board_idx}.file_idx(idx:end));
      board_hdrs{board_idx}.offset = cat(2, ...
        board_hdrs{board_idx}.offset(1:idx-1), ...
        (-2^31) * ones(1,num_skipped), ...
        board_hdrs{board_idx}.offset(idx:end));
    end
  end
  end
end

% =====================================================================
%% Basic Checks For Timing
% =====================================================================
fprintf('Determine if in a 1 second jump state at start of segment\n');
fprintf('  Note: this determination is not conclusive and additional\n');
fprintf('  external checks on gps-sync are still necessary. (%s)\n', datestr(now));

time = board_hdrs{1}.seconds + board_hdrs{1}.fractions/param.radar.fs;
good_idx = find(~isnan(time),1);
time(1) = time(good_idx) - (good_idx-1)*init_EPRI_estimate;

for time_idx = 2:length(time)
  if isnan(time(time_idx))
    time(time_idx) = time(time_idx-1) + init_EPRI_estimate;
  end
end

PRI_GUARD = 1.25;
time_diff = diff(time);
pos_sec_jumps = abs(time_diff - 1 - init_EPRI_estimate) < init_EPRI_estimate*PRI_GUARD;

neg_sec_jumps = abs(time_diff + 1 - init_EPRI_estimate) < init_EPRI_estimate*PRI_GUARD;

pos_sec_jump_idxs = find(pos_sec_jumps);
neg_sec_jump_idxs = find(neg_sec_jumps);

if isempty(pos_sec_jump_idxs) & isempty(neg_sec_jump_idxs)
  one_second_jump_state = 0;
elseif ~isempty(pos_sec_jump_idxs) & isempty(neg_sec_jump_idxs)
  one_second_jump_state = 0;
elseif isempty(pos_sec_jump_idxs) & ~isempty(neg_sec_jump_idxs)
  one_second_jump_state = 1;
elseif ~isempty(pos_sec_jump_idxs) & ~isempty(neg_sec_jump_idxs)
  one_second_jump_state = neg_sec_jump_idxs(1) < pos_sec_jump_idxs(1);
end

fprintf('Rebuilding time vector from EPRI vector and measuring frequency offset (%s)\n', datestr(now));
records.raw.epri = uint32(board_hdrs{1}.epri);
records.raw.seconds = board_hdrs{1}.seconds;
records.raw.fractions = board_hdrs{1}.fractions;

if strcmpi(param.day_seg,'20120416_02')
  fprintf('IMPLEMENTING HACK: This segment was programmed one way on the DDS\n');
  fprintf('and another way on the NI digital system.\n');
  % First waveform is bad?, presums are set to give a correct EPRI below
  hdr_master.wfs(1).num_wfs = 2;
  hdr_master.wfs(1).presums = 22;
  % Second waveform is good
  hdr_master.wfs(2).num_wfs = 2;
  % Third waveform does not exist in the file (except the headers indicate
  % it exists)
  hdr_master.wfs = hdr_master.wfs(1:2);
  
  init_EPRI_estimate = 0.003332854351482;
end

% Count the presums
num_presum = 0;
for wf = 1:length(hdr_master.wfs)
  if param.records.presum_bug_fixed
    num_presum = num_presum + hdr_master.wfs(wf).presums;
  else
    num_presum = num_presum + hdr_master.wfs(wf).presums + 1;
  end
end
EPRI = num_presum/param.radar.prf;

if param.records.use_ideal_epri
  param.records.use_ideal_epri = EPRI;
end

gps_clk_notes = '';
if all(board_hdrs{1}.seconds(~isnan(board_hdrs{1}.seconds)) == 0)
  warning('NMEA or 1 PPS signal appears to have been disconnected so that the seconds field is zero for all records. Verify this is the case and then type dbcont.')
  keyboard
  % Using EPRI to recreate time
  utc_time_sod_measured = double(records.raw.epri) * init_EPRI_estimate;
  gps_clk_notes = sprintf('  No 1 PPS or NMEA signal, reconstructing time from EPRI.\n');
else
  if param.records.file_version == 402 || param.records.file_version == 403
    utc_time_sod_measured = records.raw.seconds ...
      + records.raw.fractions/param.radar.fs - one_second_jump_state;
  elseif param.records.file_version == 404
    utc_time_sod_measured = records.raw.seconds ...
      + records.raw.fractions/(param.radar.fs/4) - one_second_jump_state;
  elseif param.records.file_version == 407 || param.records.file_version == 408
    utc_time_sod_measured = records.raw.seconds ...
      + records.raw.fractions/(param.radar.fs/8) - one_second_jump_state;
  end
end

good_idxs = find(~isnan(utc_time_sod_measured));
utc_time_sod_measured = utc_time_sod_measured(good_idxs);

%% Check for day wrap
day_wrap_idx = find(diff(utc_time_sod_measured) < -85400 & diff(utc_time_sod_measured) > -87400);
if ~isempty(day_wrap_idx)
  warning('A jump in time that looks like a day wrap has been found. Check plot. Run "dbcont" to fix the seconds of day wrap and continue.');
  figure(1); clf;
  plot(utc_time_sod_measured);
  ylabel('UTC time seconds of day');
  keyboard;
  utc_time_sod_measured(day_wrap_idx+1:end) = utc_time_sod_measured(day_wrap_idx+1:end) + 86400;
end

utc_time_sod_expected = double(records.raw.epri(good_idxs)) * init_EPRI_estimate;
utc_time_sod_expected = utc_time_sod_expected - utc_time_sod_expected(1) + utc_time_sod_measured(1) - one_second_jump_state;

%% Compare expected and measured GPS times
clock_notes = '';
if strcmpi(param.day_seg,'20110418_03')
  fprintf('IMPLEMENTING HACK: This segment is known to have bad EPRI values in files 37, 58, and 113\n');
  fprintf('because the time stamps in the filenames on the computer indicate that the\n');
  fprintf('UTC time in the files is correct, but the EPRI values are not.\n');
  
  % We do need to correct the 1 second jump offsets though
  jumps = [0 diff(utc_time_sod_expected-utc_time_sod_measured)];
  jump_mask = abs(jumps) > 0.5;
  jumps(~jump_mask) = 0;
  jump_correction = cumsum(jumps);
  jump_correction(jump_correction > -0.5) = 0;
  jump_correction(jump_correction <= -0.5) = -1;
  fixed_time = utc_time_sod_measured + jump_correction;
  
  utc_time_sod_corrected = zeros(size(records.raw.epri));
  utc_time_sod_corrected(good_idxs) = fixed_time;
  for idx = 1:length(utc_time_sod_corrected)
    if utc_time_sod_corrected(idx) == 0
      utc_time_sod_corrected(idx) = utc_time_sod_corrected(idx-1) + init_EPRI_estimate;
    end
  end
  
  p = [NaN NaN];
  clock_notes = cat(2,clock_notes,sprintf('  Clock error %.12f\n  fs %.2f Hz\n  PRF %.6f Hz\n  max_error %.2f ms\n', ...
    p(1)/EPRI, param.radar.fs / (p(1)/EPRI), param.radar.prf / (p(1)/EPRI), ...
    max(abs(utc_time_sod_corrected(good_idxs)-utc_time_sod_measured)) * 1000));
  
else
  
  if strcmpi(param.day_seg,'20111218_02') && strcmpi(param.season_name,'2011_Antarctica_TO')
    % Hack to fix 1 second offset that is not properly corrected
    one_second_jump_idxs = abs(utc_time_sod_expected-utc_time_sod_measured-1) < 0.13;
    utc_time_sod_measured(one_second_jump_idxs) = utc_time_sod_measured(one_second_jump_idxs) + 1;
  end
  
  if all(abs(utc_time_sod_expected-utc_time_sod_measured) < 0.13 ...
      | abs(utc_time_sod_expected-utc_time_sod_measured+1) < 0.13)
    %% This handles the case where all times are good EXCEPT for 1 second
    % jumps which it corrects.
    one_second_jump_idxs = abs(utc_time_sod_expected-utc_time_sod_measured+1) < 0.13;
    utc_time_sod_measured(one_second_jump_idxs) = utc_time_sod_measured(one_second_jump_idxs) - 1;
    figure(1); clf;
    plot(utc_time_sod_expected-utc_time_sod_measured);
    ylabel('Mismatch (sec)');
    title('Mismatch between expected and measured UTC time');
    xlabel('Records');
    good_mask = abs(utc_time_sod_expected-utc_time_sod_measured) <= 0.1;
    good_idxs = good_idxs(good_mask);
    utc_time_sod_measured = utc_time_sod_measured(good_mask);
  elseif any(abs(utc_time_sod_expected-utc_time_sod_measured) > 0.1)
    figure(1); clf;
    plot(utc_time_sod_expected-utc_time_sod_measured);
    ylabel('Mismatch (sec)');
    title('Mismatch between expected and measured UTC time');
    xlabel('Records');
    figure(2); clf;
    plot(diff(utc_time_sod_measured) ./ diff(double(records.raw.epri(good_idxs))),'.');
    ylabel('Estimated EPRI (sec)');
    xlabel('Records');
    epri_double = double(records.raw.epri(good_idxs));
    final_EPRI_estimate = (utc_time_sod_measured(end)-utc_time_sod_measured(1)) / (epri_double(end)-epri_double(1));
    title(sprintf('param.radar EPRI %.8f ms\ninit_EPRI_estimate %.8f ms\nfinal_EPRI_estimate %.8f ms', ...
      EPRI*1e3, init_EPRI_estimate*1e3, final_EPRI_estimate*1e3),'interpreter','none');
    warning('The expected and measured times are off by > 0.1.');
    fprintf('Verify in the plot that all differences are less than 0.1 seconds\n');
    fprintf('except a few outliers. dbcont replaces these outliers. However, if\n');
    fprintf('there are large sections that are greater than 0.1 seconds (and not\n');
    fprintf('just a few spikes from bad record headers), there might be a systematic\n');
    fprintf('problem that should be debugged. Some situations (including total\n');
    fprintf('failure of the 1 PPS) is solved by using the expected time completely):\n');
    fprintf('    utc_time_sod_measured = utc_time_sod_expected;\n');
    fprintf('For total failure of 1 PPS, you also have to set \n');
    fprintf('param.vectors.gps.time_offset properly.\n');
    fprintf('If the init_EPRI_estimate is off by a little, there will be\n');
    fprintf('a slowly growing error between expected and measured UTC time.\n');
    fprintf('You may need to set UTC_MAX_ERROR to a larger value to allow for this\n');
    fprintf('or use param.records.use_ideal_epri to help find a better EPRI.\n');
    UTC_MAX_ERROR = 0.1;
    good_mask = abs(utc_time_sod_expected-utc_time_sod_measured) <= UTC_MAX_ERROR;
    good_percent = sum(good_mask)/length(good_mask);
    clock_notes = cat(2,clock_notes,sprintf('Mismatch between expected and measured UTC time:\n'));
    clock_notes = cat(2,clock_notes,sprintf('  %.4f%%\n', good_percent*100));
    fn_fig = ct_filename_ct_tmp(param,'','records', ['mismatch_UTC_time.fig']);
    fprintf('Saving %s\n', fn_fig);
    [fn_fig_dir,name] = fileparts(fn_fig);
    if ~exist(fn_fig_dir,'dir')
      mkdir(fn_fig_dir);
    end
    saveas(1,fn_fig);
    if good_percent < 0.995
      keyboard
    end
    good_mask = abs(utc_time_sod_expected-utc_time_sod_measured) <= UTC_MAX_ERROR;
    good_idxs = good_idxs(good_mask);
    utc_time_sod_measured = utc_time_sod_measured(good_mask);
  end
  
  % utc_time_sod_measured = records.raw.seconds ...
  %   + records.raw.fractions/param.radar.fs - one_second_jump_state;
  % utc_time_sod_measured = utc_time_sod_measured(good_idxs);
  
  p = polyfit(double(records.raw.epri(good_idxs)),utc_time_sod_measured,1);
  utc_time_sod_corrected = double(records.raw.epri)*p(1) + p(2);
  clock_notes = cat(2,clock_notes,sprintf('  Clock error %.12f\n  fs %.2f Hz\n  PRF %.6f Hz\n  max_error %.2f ms\n', ...
    p(1)/EPRI, param.radar.fs / (p(1)/EPRI), param.radar.prf / (p(1)/EPRI), ...
    max(abs(utc_time_sod_corrected(good_idxs)-utc_time_sod_measured)) * 1000));
  fprintf('%s',clock_notes);
  plot(utc_time_sod_corrected(good_idxs), utc_time_sod_corrected(good_idxs)-utc_time_sod_measured);
  xlabel('UTC time seconds of day (sec)');
  ylabel('Mismatch (sec)');
  title('Mismatch between corrected and measured UTC time','FontSize',10,'FontWeight','normal');
  pause(0.1);
  drawnow;
  if isfield(param.records,'debug_level') && ~isempty(param.records.debug_level) ...
      && param.records.debug_level > 1
    fprintf('Please review the time clock comparison and then run "dbcont"\n');
    keyboard;
  end
  
end

% =====================================================================
%% Apply time correction to UTC time read in from radar files
% =====================================================================

% ===================================================================
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

% =====================================================================
% Save records files
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
records.ver = 3;
%  records.raw.epri(1 .. Nx)
%  records.raw.seconds(1 .. Nx)
%  records.raw.fraction(1 .. Nx)

records.settings = [];

records.settings.wfs_records = 1;
records.settings.wfs = hdr_master.wfs;

records.notes = cat(2,sprintf('\nEPRI NOTES\n%s',epri_notes), ...
  sprintf('\nGPS_CLK NOTES\n%s',gps_clk_notes), ...
  sprintf('\nCLOCK NOTES\n%s',clock_notes));
records.param_records = param;

fprintf('Saving records file %s (%s)\n',records_fn,datestr(now));
save(records_fn,'-v7.3','-struct','records'); % Handle large file sizes, so use v7.3

% =====================================================================
% Create record aux files for faster loading times
% =====================================================================

fprintf('Creating auxiliary records files %s (%s)\n',records_fn,datestr(now));
create_records_aux_files(records_fn);

fprintf('Done (%s)\n\n', datestr(now));

return;

