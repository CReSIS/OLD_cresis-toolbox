% create_records_fmcw_accum_sync
%
% Script called from create_records_fmcw_accum
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
  
  new_param = read_param_xls(ct_filename_param('snow_param_2015_Greenland_Polar6.xls'),'20150911_02');

  fn = ct_filename_tmp(new_param,new_param.records.records_fn,'records','workspace');
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

% ======================================================================

if ~isfield(param.records,'manual_time_correct') || isempty(param.records.manual_time_correct)
  param.records.manual_time_correct = 0;
end

% Initialize notes string that will be stored in records file
radar_time_notes = '';
epri_notes = '';
clock_notes = '';

%% Loop through each board header
for board_idx = 1:length(board_hdrs)
  % Create an alias/pointer of the data so that code is more readable... at the end
  % of the loop this will be written back into board_hdrs.
  hdr = board_hdrs{board_idx};
  radar_time_notes = cat(2,radar_time_notes, sprintf('Board Index %d\n', board_idx));
  epri_notes = cat(2,epri_notes, sprintf('Board Index %d\n', board_idx));
  clock_notes = cat(2,clock_notes, sprintf('Board Index %d\n', board_idx));
    
    %% Check to see if there are big seconds jumps (those encountered in 2016_Greenland_P3,
    %% jump to a large number and then drop back to correct values)?
  big_sec_jump_idxs = find(abs(diff(double(hdr.seconds)))>1e5);
  fraction_wrap_idxs = find(diff(double(hdr.fraction))<0);
  mid_jumps = find(~ismember(big_sec_jump_idxs,fraction_wrap_idxs));
  if ~isempty(mid_jumps)
      [tmp,tmp_idxs] = min(abs(fraction_wrap_idxs - big_sec_jump_idxs(mid_jumps)));
      big_sec_jump_idxs(mid_jumps) = fraction_wrap_idxs(tmp_idxs);
  end
  big_sec_jump_idxs = big_sec_jump_idxs(ismember(big_sec_jump_idxs,fraction_wrap_idxs));
  if ~isempty(big_sec_jump_idxs)
    warning('Header seconds jump more than 1e5 sec, correcting jumps');
    if 0
      figure(101);plot(hdr.fraction);
      max_fraction = max(hdr.fraction);
      for idx = 1:length(big_sec_jump_idxs)
        figure(101);hold on;plot([big_sec_jump_idxs(idx),big_sec_jump_idxs(idx)]+1,[0,1.2*max_fraction],'r--');
      end
    end
    for idx = 1:2:length(big_sec_jump_idxs) - mod(length(big_sec_jump_idxs),2)
      jump_start = big_sec_jump_idxs(idx)+1;
      wrap_idxs = [find(fraction_wrap_idxs == big_sec_jump_idxs(idx)):find(fraction_wrap_idxs == big_sec_jump_idxs(idx+1))];
      wrap_idxs(1) = [];
      for wrap_idx = wrap_idxs
        hdr.seconds(jump_start:fraction_wrap_idxs(wrap_idx)) = hdr.seconds(jump_start-1)+1;
        jump_start = fraction_wrap_idxs(wrap_idx) + 1;
      end
    end
    if mod(length(big_sec_jump_idxs),2)
      jump_start = big_sec_jump_idxs(length(big_sec_jump_idxs))+1;
      hdr.seconds(jump_start:length(hdr.seconds)) = hdr.seconds(jump_start-1)+1; 
    end
  end
  
  %% Create hdr.utc_time_sod vector
  if any(param.records.file_version == [8 101])
    hdr.utc_time_sod = double(hdr.seconds) + 2*double(hdr.fraction)/param.radar.fs;
  else
    hdr.utc_time_sod = double(hdr.seconds) + double(hdr.fraction)/param.radar.fs;
  end
  if strcmpi(param.season_name,'2014_Greenland_P3') && (strcmpi(param.day_seg,'20140421_02') ||...
      strcmpi(param.day_seg,'20140423_01') || strcmpi(param.day_seg,'20140502_00') ||...
      strcmpi(param.day_seg,'20140508_02'))
    first_file = board_fns{board_idx}{file_idxs(1)};
    hdr.utc_time_sod = hdr.utc_time_sod - hdr.utc_time_sod(1);
    hdr.utc_time_sod = hdr.utc_time_sod + str2num(first_file(end-14:end-13))*3600 +...
      str2num(first_file(end-12:end-11))*60 + str2num(first_file(end-10:end-9));
  end
  
  
  %% Find EPRI (effective pulse repetition interval) by radar name
  if strcmp(param.radar_name,'accum')
    num_bands = 16;
    hdr.wfs{1}.presums = 8;
  elseif any(strcmp(param.radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5','snow8'}))
    num_bands = 1;
  end
  
  if ~isfield(param.radar.wfs(1),'presum_override') || isempty(param.radar.wfs(1).presum_override)
    presums = 0;
    for wf = 1:length(hdr.wfs{1})
      presums = presums + hdr.wfs{1}(wf).presums;
    end
    EPRI = num_bands*presums/param.radar.prf;
  else
    warning('Presum override turned on (file presums %i, override %i)\n', hdr.wfs{1}.presums, param.radar.wfs(1).presum_override);
    EPRI = num_bands*param.radar.wfs(1).presum_override/param.radar.prf;
  end
  
  %% Check to see if 1 PPS signal was connected?
  if all(hdr.seconds == 0)
    warning('Header seconds field is all zero. This often means no 1 PPS signal. It can also mean there was no NMEA data.');
    plot(hdr.utc_time_sod);
    no_1PPS_signal = false;
    no_NMEA_signal = false;
    fprintf('Set "no_1PPS_signal = true" or "no_NMEA_signal = true" and then run "dbcont" to do correction\n');
    keyboard
    
    TIME_JUMP_THRESHOLD = 1.25; % In units of EPRI (1 = EPRI, 2 = 2*EPRI, etc)
    time_jumps_relative = (diff(hdr.utc_time_sod) ./ diff(double(hdr.epri)) - EPRI)/EPRI;
    time_jump_idxs = find(abs(time_jumps_relative) > TIME_JUMP_THRESHOLD);
    
    if no_1PPS_signal
      radar_time_notes = cat(2,radar_time_notes, ...
        sprintf('Fixing no 1 PPS signal.\n', length(time_jump_idxs)));
      % There was no 1 PPS signal so the fractions counter would keep wrapping
      % through out the day
      EPRI_estimated = median(diff(hdr.utc_time_sod));
      % 1. Find all the places that the fractions wrapped (32 bit roll over)
      diff_fraction = diff(hdr.fraction);
      diff_utc_time_sod = diff(hdr.utc_time_sod);
      wrap_idxs = find(diff_fraction < -2^31);
      % 2. Unwrap these
      for wrap_idx = wrap_idxs
        hdr.utc_time_sod(wrap_idx:end) = hdr.utc_time_sod(wrap_idx:end) - diff_utc_time_sod(wrap_idx) + EPRI_estimated;
      end
    end
    if no_NMEA_signal
      radar_time_notes = cat(2,radar_time_notes, ...
        sprintf('Fixing no NMEA signal.\n', length(time_jump_idxs)));
      % There was no NMEA signal so the fractions counter would get reset for
      % each 1 PPS, but seconds was not getting incremented.
      EPRI_estimated = median(diff(hdr.utc_time_sod));
      % 1. Find all the places that the 1 PPS caused a fraction reset without a second increase
      diff_utc_time_sod = diff(hdr.utc_time_sod);
      wrap_idxs = find(diff_utc_time_sod < -0.75);
      % 2. Unwrap these
      for wrap_idx = wrap_idxs
        hdr.utc_time_sod(wrap_idx+1:end) = hdr.utc_time_sod(wrap_idx+1:end) - diff_utc_time_sod(wrap_idx) + EPRI_estimated;
      end
    end
  end
  

  % =====================================================================
  %% Find bad UTC time SOD and EPRI entries
  % =====================================================================
  if board_idx == 1
    orig_hdr = hdr;
  end
  
  if param.records.file_version == 101
    hdr.epri = round((hdr.utc_time_sod - min(hdr.utc_time_sod)) / EPRI);
  end
  
  %% Quick Load using just first and last record
  if any(strcmpi(param.season_name,{'2013_Greenland_P3','2013_Antarctica_Basler','2013_Antarctica_P3','2014_Alaska_TOnrl','2014_Greenland_P3'}))
    warning('Applying 2013 Greenland P3 EPRI HACK');
    % 2013 Greenland P3 EPRI HACK: If this is a problem, a new column
    % should be added to records worksheet to enable this bug fix. This
    % is a temporary solution until then.
    % Determine if first EPRI is behind by 1
    bad_idx = find(diff(double(hdr.epri)) == 0);
    hdr.epri(bad_idx) = hdr.epri(bad_idx) - 1;
    epri_notes = cat(2,epri_notes,'Applied hack to fix epri counting 0,0,2,2,4,4,6, etc\n');
  elseif strcmpi(param.season_name,'2010_Greenland_DC8') ...
      && strcmpi(param.radar_name,'snow') ...
      && strcmpi(param.day_seg,'20100405_02')
    % 2010 Greenland DC8 segment 20100405_02 has a large negative EPRI
    % jump that appears to be due to a random jump in the EPRI since the time
    % vectors shows no problems
    bad_idx = find(diff(double(hdr.epri)) < -1e6,1) + 1;
    hdr.epri(bad_idx:end) = hdr.epri(bad_idx:end) + hdr.epri(bad_idx-1) - hdr.epri(bad_idx) + 1;
    epri_notes = cat(2,epri_notes,'Applied hack to remove random jump in 20100405_02\n');
  end
  
  %% Correct for digital errors in epri
  fprintf('Searching for and correcting jumps in EPRI (%s)\n', datestr(now));
  epri_diff = diff(double(hdr.epri));
  drop_last_files = 0;
  if any(epri_diff > 1000)
    warning('Large positive jump. Sometimes these are single header errors (one record has an EPRI header error that can be fixed with "random_large_epri_offset=true" and then run "dbcont", possibly adjusting random_large_epri_offset_tresholds = [1 3] (largest and smallest allowable epri jumps, usually adjust to [-M N] so that M/N exclude all EPRI jumps that are real jumps rather than just header errors), and a note in cmd.notes) and sometimes it means the segment needs to be broken into two at the file idx listed because there is a data gap > 1km (typically the file with the epri jump in it is left out so that segment from 1 to N is broken into 1 to L-1 and L+1 to N). If the EPRI jump is real and is small enough to be ignored, just run dbcont.');
    figure(1); clf;
    plot(epri_diff);
    title('diff of EPRI');
    a1 = gca;
    figure(2); clf;
    plot(hdr.epri);
    title('EPRI');
    a2 = gca;
    figure(3); clf;
    plot(hdr.utc_time_sod);
    title('UTC time');
    a3 = gca;
    linkaxes([a1 a2 a3],'x');
    
    random_large_epri_offset_tresholds = [1 3];
    epri_diff_big_jumps = find(~(epri_diff >= random_large_epri_offset_tresholds(1) & epri_diff <= random_large_epri_offset_tresholds(2)));
    fprintf('Printing %d of %d epri jumps outside thresholds\n', min(1000,length(epri_diff_big_jumps)), length(epri_diff_big_jumps));
    for idx = 1:min(1000,length(epri_diff_big_jumps))
      file_idx = find(records.relative_rec_num{board_idx} < epri_diff_big_jumps(idx),1,'last');
      num_jumped = epri_diff(epri_diff_big_jumps(idx));
      fprintf('rec idx %d, file idx %d %s jump of %d at %d\n', ...
        idx, file_idxs(file_idx), records.relative_filename{board_idx}{file_idx}, ...
        num_jumped, epri_diff_big_jumps(idx));
    end
    
    random_large_epri_offset = 0;
    keyboard
    
    if random_large_epri_offset
      % In case of large offset to epri for short segments of data
      % Used to fix kuband2: 20120319_02
      % Used to fix snow2: 20120319_02
      
      epri_diff_big_jumps = find(~(epri_diff >= random_large_epri_offset_tresholds(1) & epri_diff <= random_large_epri_offset_tresholds(2)));
      fprintf('Fixing %d epri jumps outside thresholds\n', length(epri_diff_big_jumps));
      for idx = 1:min(100,length(epri_diff_big_jumps))
        file_idx = find(records.relative_rec_num{board_idx} < epri_diff_big_jumps(idx),1,'last');
        num_jumped = epri_diff(epri_diff_big_jumps(idx));
        epri_notes = cat(2,epri_notes, ...
          sprintf('EPRI jump fixed rec idx %d, file idx %d %s jump of %d at %d\n', ...
          idx, file_idxs(file_idx), records.relative_filename{board_idx}{file_idx}, ...
          num_jumped, epri_diff_big_jumps(idx)));
        fprintf('rec idx %d, file idx %d %s jump of %d at %d\n', ...
          idx, file_idxs(file_idx), records.relative_filename{board_idx}{file_idx}, ...
          num_jumped, epri_diff_big_jumps(idx));
      end
      
      epri_diff(epri_diff_big_jumps) = 1;
      epri_new = cumsum([double(hdr.epri(1)) epri_diff]);
      correction = double(hdr.epri) - epri_new;
      epri_new(abs(correction)<10) = hdr.epri(abs(correction)<10);
      figure(1); clf;
      set(1,'Renderer','painters');
      plot(double(hdr.epri) - epri_new,'x');
      ylim([-20 20]);
      title('Is it all zeros for the good records?');
      figure(2); clf;
      plot(epri_new);
      title('Is it reasonable for EPRI?');
      figure(3); clf;
      plot(diff(epri_new));
      title('Are the differences all 1,2,3?');
      warning('Please look at plots to make sure correction worked before running "dbcont"');
      keyboard
      hdr.epri = int32(epri_new);
    end
    if exist('drop_last_files','var') && drop_last_files == 1
      % Rerun the function with the correct variables in the param sheet
      % and then run "drop_last_files = 1".
      keyboard
      bad_files = file_idxs > param.vectors.file.stop_idx;
      num_dropped = sum(bad_files);
      records.relative_filename{board_idx} = records.relative_filename{board_idx}(1:end-num_dropped);
      records.relative_rec_num{board_idx} = records.relative_rec_num{board_idx}(1:end-num_dropped);
      hdr.wfs = hdr.wfs(1:end-num_dropped);
      hdr.seconds = hdr.seconds(1:records.relative_rec_num{board_idx}(end)-1);
      hdr.fraction = hdr.fraction(1:records.relative_rec_num{board_idx}(end)-1);
      hdr.epri = hdr.epri(1:records.relative_rec_num{board_idx}(end)-1);
      hdr.offset = hdr.offset(1:records.relative_rec_num{board_idx}(end)-1);
      hdr.utc_time_sod = hdr.utc_time_sod(1:records.relative_rec_num{board_idx}(end)-1);
    end
  end
  epri_diff = diff(double(hdr.epri));
  if (any(strcmp(param.day_seg,{'20111012_04','20111023_08'})) ...
      && strcmp(param.radar_name,'kuband')) ...
      || (any(strcmp(param.day_seg,{'20090406_03','20090501_05','20090502_03'})) ...
      && strcmp(param.radar_name,'snow'))
    % HACK to fix random large negative jump in EPRI that should not be there.
    % utc_time_sod is fine, but EPRI makes large negative jump at one position
    % as if timing was not met in the FPGA and the register for EPRI messed up.
    warning('HACK for fixing epri with large negative EPRI jump');
    epri_notes = cat(2,epri_notes,'Applied hack to fix large negative EPRI jump\n');
    figure(1); clf;
    plot(hdr.epri,'r--')
    neg_jump_idxs = find(epri_diff < 0)
    meas_epri_diff = diff(double(hdr.epri(neg_jump_idxs + (-1:1))))
    actual_epri_diff = round(diff(double(hdr.utc_time_sod(neg_jump_idxs + (-1:1))) / EPRI))
    hdr.epri(neg_jump_idxs + 1 : end) = hdr.epri(neg_jump_idxs + 1 : end) + actual_epri_diff(2) - meas_epri_diff(2)
    hold on;
    plot(hdr.epri)
    hold off;
    legend('Old','New','Location','Northwest')
    figure(2); clf;
    plot(diff(hdr.utc_time_sod))
    fprintf('If hdr.epri and hdr.utc_time_sod look correct, just run "dbcont"\n');
    epri_diff = diff(double(hdr.epri));
    keyboard
  end
  if any(epri_diff < 1)
    warning('Nonpositive EPRI jumps');
    epri_diff_neg_jumps = find(epri_diff < 1);
    fprintf('Printing %d of %d nonpositive jumps. Run "dbcont" to throw out all repeated epri associated with these jumps.\n', ...
      min(10,length(epri_diff_neg_jumps)), length(epri_diff_neg_jumps));
    for idx = 1:min(10,length(epri_diff_neg_jumps))
      file_idx = find(records.relative_rec_num{board_idx} < epri_diff_neg_jumps(idx),1,'last');
      num_jump = epri_diff(epri_diff_neg_jumps(idx));
      fprintf('rec idx %d, file idx %d %s nonpositive jump of %d at %d\n', ...
        idx, file_idxs(file_idx), records.relative_filename{board_idx}{file_idx}, ...
        num_jump, epri_diff_neg_jumps(idx));
    end
    keyboard
    % First write notes to epri_notes
    epri_diff_neg_jumps = find(epri_diff < 1);
    for idx = 1:length(epri_diff_neg_jumps)
      file_idx = find(records.relative_rec_num{board_idx} < epri_diff_neg_jumps(idx),1,'last');
      epri_notes = cat(2,epri_notes,sprintf('rec idx %d, file idx %d %s\n', idx, file_idxs(file_idx), records.relative_filename{board_idx}{file_idx}));
    end
    % Fix negative jumps, by throwing out all records with these jumps
    cur_epri = hdr.epri(1);
    epri_mask = ones(size(hdr.epri));
    for idx = 2:length(hdr.epri)
      if hdr.epri(idx) <= cur_epri
        epri_mask(idx) = 0;
      else
        cur_epri = hdr.epri(idx);
      end
    end
    file_mask = ones(size(records.relative_filename{board_idx}));
    for idx = 1:length(records.relative_rec_num{board_idx})-1
      recs = records.relative_rec_num{board_idx}(idx:idx+1);
      if all(epri_mask(recs(1):recs(2)-1) == 0)
        % Remove file from listing
        file_mask(idx) = 0;
      else
        records.relative_rec_num{board_idx}(idx) = sum(epri_mask(1:records.relative_rec_num{board_idx}(idx)-1)) + 1;
      end
    end
    epri_mask = logical(epri_mask);
    file_mask = logical(file_mask);
    
    hdr.seconds = hdr.seconds(epri_mask);
    hdr.fraction = hdr.fraction(epri_mask);
    hdr.epri = hdr.epri(epri_mask);
    hdr.offset = hdr.offset(epri_mask);
    hdr.utc_time_sod = hdr.utc_time_sod(epri_mask);
    hdr.file_idx = hdr.file_idx(epri_mask);
    
    if param.records.file_version ~= 101
      hdr.wfs = hdr.wfs(file_mask);
    end
    records.relative_filename{board_idx} = records.relative_filename{board_idx}(file_mask);
    records.relative_rec_num{board_idx} ...
      = [records.relative_rec_num{board_idx}(file_mask) length(hdr.seconds)+1];
  end
  
  fprintf('Searching for and correcting jumps in UTC time (%s)\n', datestr(now));
  
  
  %% Check for seconds of day roll over and unwrap (assume jump backward
  % of more than 23 hours is a roll over)
  wrap_idxs = find(abs(diff(hdr.utc_time_sod) + 86400) < 3600);
  for wrap_idx = wrap_idxs
    hdr.utc_time_sod(wrap_idx+1:end) = hdr.utc_time_sod(wrap_idx+1:end) + 86400;
    radar_time_notes = cat(2,radar_time_notes, ...
      sprintf('%d: Applied day wrap correction\n', wrap_idx));
  end
  
  %% Correct timing errors
  
  TIME_JUMP_THRESHOLD = 1.25; % In units of EPRI (1 = EPRI, 2 = 2*EPRI, etc)
  time_jumps_relative = (diff(hdr.utc_time_sod) ./ diff(double(hdr.epri)) - EPRI)/EPRI;
  time_jump_idxs = find(abs(time_jumps_relative) > TIME_JUMP_THRESHOLD);
  
  if ~isempty(time_jump_idxs)
    radar_time_notes = cat(2,radar_time_notes, ...
      sprintf('Time error threshold exceeded %d times.\n', length(time_jump_idxs)));
    warning('Time error threshold exceeded %d times.', length(time_jump_idxs));
  end
  
  ignore_idxs = [];
  while ~isempty(time_jump_idxs)
    fprintf('%d jumps left to fix\n', length(time_jump_idxs));
    %length(time_jump_idxs) % Useful for debugging
    bad_idx = time_jump_idxs(1)+1;
    if hdr.seconds(bad_idx)-hdr.seconds(bad_idx-1) == 1 && abs(hdr.utc_time_sod(bad_idx)-hdr.utc_time_sod(bad_idx-1) - 1) < 4*EPRI
      %% 1 PPS timing jump error
      hdr.seconds(bad_idx:end) = hdr.seconds(bad_idx:end) - 1;
      next_jump = find(diff(hdr.seconds(bad_idx:end)),1);
      if param.vectors.gps.utc_time_halved
        % Need to complete this
        keyboard
      else
        delta_fraction = hdr.fraction(bad_idx-1) + EPRI*param.radar.fs - hdr.fraction(bad_idx);
        hdr.fraction(bad_idx + (0:next_jump-1)) = hdr.fraction(bad_idx + (0:next_jump-1)) + delta_fraction;
        hdr.utc_time_sod = double(hdr.seconds) + double(hdr.fraction)/param.radar.fs;
      end
      radar_time_notes = cat(2,radar_time_notes, ...
        sprintf('%d: PPS_error\n', bad_idx));
    else
      %% Other timing jump error
      if param.records.manual_time_correct
        fprintf('Manual time correction enabled: fix hdr.utc_time_sod manually and then run dbcont\n');
        keyboard
        if 0
          % Time jump is real (do nothing)
          ignore_idxs = unique([ignore_idxs time_jump_idxs(1)]);
        end
        if 0
          % Time after jump is believed to be correct
          predicted_jump = test_time(bad_idx)-test_time(bad_idx-1);
          measured_jump = hdr.utc_time_sod(bad_idx) - hdr.utc_time_sod(bad_idx-1);
          hdr.utc_time_sod(1:bad_idx-1) = hdr.utc_time_sod(1:bad_idx-1) + measured_jump - predicted_jump;
        end
        if 0
          % Time before jump is believed to be correct
          predicted_jump = test_time(bad_idx)-test_time(bad_idx-1);
          measured_jump = hdr.utc_time_sod(bad_idx) - hdr.utc_time_sod(bad_idx-1);
          hdr.utc_time_sod(bad_idx:end) = hdr.utc_time_sod(bad_idx:end) - measured_jump + predicted_jump;
        end
        radar_time_notes = cat(2,radar_time_notes, ...
          sprintf(' %d: NaN manual_correction\n', bad_idx));
      else
        MIN_EPRI_BEFORE_FIX = 10;
        if abs(time_jumps_relative(bad_idx-1)) < MIN_EPRI_BEFORE_FIX
          % This jump is less than MIN_EPRI_BEFORE_FIX EPRI wide and may
          % be just a GPS/radar time update correction that is real
          % (i.e. radar time drifted and was corrected... not sure if this
          % is really the case)
          radar_time_notes = cat(2,radar_time_notes, ...
            sprintf(' %d: %f ignored\n', bad_idx, time_jumps_relative(bad_idx-1)));
          ignore_idxs = unique([ignore_idxs time_jump_idxs(1)]);
        else
          % Default Correction: Time before jump is believed to be correct
          figure(1); clf;
          plot(hdr.utc_time_sod);
          ylabel('UTC time (sec)');
          hold on;
          plot(bad_idx, hdr.utc_time_sod(bad_idx),'ro');
          hold off;
          figure(2); clf;
          plot(time_jumps_relative);
          ylabel('UTC time jump (EPRIs)');
          hold on;
          plot(bad_idx, time_jumps_relative(bad_idx-1),'ro');
          hold off;
          drawnow;
          warning('Verify in the plots that this is a time jump error before dbcont\n');
          %         keyboard
          hdr.utc_time_sod(bad_idx:end) = hdr.utc_time_sod(bad_idx:end) - time_jumps_relative(bad_idx-1)*EPRI;
          radar_time_notes = cat(2,radar_time_notes, ...
            sprintf(' %d: %f applied \n', bad_idx, time_jumps_relative(bad_idx-1)));
        end
      end
    end
    
    TIME_JUMP_THRESHOLD = 1.25; % In units of EPRI (1 = EPRI, 2 = 2*EPRI, etc)
    time_jumps_relative = (diff(hdr.utc_time_sod) ./ diff(double(hdr.epri)) - EPRI)/EPRI;
    time_jump_idxs = find(abs(time_jumps_relative) > TIME_JUMP_THRESHOLD);
    time_jump_idxs = setdiff(time_jump_idxs,ignore_idxs);
  end
  
  fprintf(radar_time_notes);
  
  p = polyfit(double(hdr.epri),hdr.utc_time_sod,1);
  utc_time_sod_corrected = double(hdr.epri)*p(1) + p(2);
  cur_board_clock_notes  ...
    = sprintf('Clock error %.12f\n  fs %.2f Hz\n  PRF %.6f Hz\n  max_error %.2f ms\n', ...
    p(1)/EPRI, param.radar.fs / (p(1)/EPRI), param.radar.prf / (p(1)/EPRI), ...
    max(abs(utc_time_sod_corrected-double(hdr.utc_time_sod))) * 1000);
  clock_notes = cat(2,clock_notes,cur_board_clock_notes);
  fprintf(cur_board_clock_notes);
  figure(1); clf;
  plot(utc_time_sod_corrected, utc_time_sod_corrected-double(hdr.utc_time_sod));
  title('UTC time correction error (should be < few ms)','fontsize',10);
  pause(0.1);
  drawnow;
  if ~isfield(param.records,'debug_level') || isempty(param.records.debug_level) || param.records.debug_level > 1
    fprintf('Please review the time clock comparison and then run "dbcont"\n');
    keyboard;
  end
  hdr.utc_time_sod = utc_time_sod_corrected;
  
  % Update the current board header
  board_hdrs{board_idx} = hdr;
end

% =====================================================================
%% Make sure all channels are aligned. Remove records that do not occur
% in every channel
% =====================================================================

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
  board_hdrs{board_idx}.fraction  = board_hdrs{board_idx}.fraction(good_idxs);
  board_hdrs{board_idx}.utc_time_sod = board_hdrs{board_idx}.utc_time_sod(good_idxs);
  board_hdrs{board_idx}.file_idx = board_hdrs{board_idx}.file_idx(good_idxs);
  board_hdrs{board_idx}.offset = board_hdrs{board_idx}.offset(good_idxs);
end

records.relative_rec_num = [];
records.offset = [];
for board_idx = 1:length(board_hdrs)
  unique_file_idxs = unique(board_hdrs{board_idx}.file_idx);
  for file_idx = 1:length(unique_file_idxs)
    records.relative_rec_num{board_idx}(file_idx) = find(board_hdrs{board_idx}.file_idx == unique_file_idxs(file_idx),1);
  end
  records.offset(board_idx,:) = board_hdrs{board_idx}.offset;
end

hdr = board_hdrs{board_idx};

% ===================================================================
%% Correlate GPS with radar data
% ===================================================================
fprintf('Loading GPS data (%s)\n', datestr(now));

if param.records.gps.en
  records = sync_radar_to_gps(param,records,hdr.utc_time_sod);
  
else
  records.lat = NaN*zeros(size(hdr.utc_time_sod));
  records.lon = NaN*zeros(size(hdr.utc_time_sod));
  records.elev = NaN*zeros(size(hdr.utc_time_sod));
  records.gps_time = NaN*zeros(size(hdr.utc_time_sod));
  records.roll = NaN*zeros(size(hdr.utc_time_sod));
  records.pitch = NaN*zeros(size(hdr.utc_time_sod));
  records.heading = NaN*zeros(size(hdr.utc_time_sod));
  records.gps_source = 'NA';
end

% =====================================================================
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
%  .relative_filename{board_idx == 1}{1 .. Nf(1)}
%  .relative_rec_num{board_idx == 1}(1 .. Nf(1))
for board_idx=1:length(board_hdrs)
  records.offset(board_idx,:) = board_hdrs{board_idx}.offset;
end
records.radar_name = param.radar_name;
records.ver = 3;
if param.records.file_version ~= 101
  records.raw.epri = orig_hdr.epri;
end
records.raw.seconds = orig_hdr.seconds;
records.raw.fraction = orig_hdr.fraction;

records.settings = [];

% Create the first entry in the records.settings field
records.settings.wfs_records(1) = 1;
for wf = 1:length(hdr.wfs{1})
  records.settings.wfs(1).wfs(wf).presums = hdr.wfs{1}(wf).presums;
end
if param.records.file_version == 1
  records.settings.wfs(1).wfs(1).num_sam = hdr.wfs{1}.num_sam;
end

records.notes = cat(2,sprintf('\nEPRI NOTES\n%s',epri_notes), ...
  sprintf('\nCLOCK NOTES\n%s',clock_notes), ...
  sprintf('\nRADAR TIME NOTES\n%s',radar_time_notes));
records.param_records = param;

fprintf('Saving records file %s (%s)\n',records_fn,datestr(now));
save(records_fn,'-v6','-struct','records');

% =====================================================================
% Create record aux files for faster loading times
% =====================================================================

fprintf('Creating auxiliary records files %s (%s)\n',records_fn,datestr(now));
create_records_aux_files(records_fn);


fprintf('Done (%s)\n\n', datestr(now));

