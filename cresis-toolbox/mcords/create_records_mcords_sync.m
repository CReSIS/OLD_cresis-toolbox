% script create_records_mcords_sync
%
% REMEMBER:
% 1. Before using for debugging a failed record creation, MAKE A LOCAL COPY
% 2. Clear local workspace "clear"
%
% This is a script designed to ease debugging of create_records_mcords.
%
% It is run from create_records_mcords.  If create_records_mcords fails
% for some reason, this script can be run from the command line to debug
% the problem.  Note that two create_records_mcords functions
%
% The variable "hdrs" must be cleared from memory before you can debug
% with this script though... the script uses the absence of this variable
% to trigger the load of the workspace.
%
% Author: John Paden

dbstack_info = dbstack;
reprocessing_mode = false;
if ~exist('param','var') || isempty(param) || length(dbstack_info) == 1
  % =====================================================================
  % Debug Setup
  % =====================================================================
  new_param = read_param_xls('/users/paden/scripts/branch/params-cr1/rds_param_2009_Antarctica_DC8.xls','20091028_02');
  %   new_param = read_param_xls('/users/paden/scripts/branch/params-cr1/rds_param_2011_Greenland_TO.xls','20110419_01');
  
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
  
  reprocessing_mode = true;
end

fprintf('Running %s correction and gps sync (%s)\n', param.day_seg, datestr(now));

%% Create radar's UTC time SOD from seconds/fractions
for adc_idx = 1:length(hdrs.seconds)
  hdrs.time{adc_idx} = double(hdrs.seconds{adc_idx}) ...
    + double(hdrs.fractions{adc_idx})/(param.radar.fs/2);
end

if reprocessing_mode
  % We are reprocessing the data (this call was not made from
  % create_records_mcords, but rather from the command line)
  
  if 0
    % Use this to reselect file range after you have discovered files
    % need to be removed (i.e. especially useful when there are huge
    % gaps in the EPRI)
    
    % First specify range of file indexes to keep (inclusive/closed
    % range so the [2 4] includes files 2 and 4, use of +/- inf valid)
    good_file_idxs = [2 inf];
    
    % Automated
    for adc_idx = 1:length(hdrs.epri)
      good_idxs = find(hdrs.file_idx{adc_idx} >= good_file_idxs(1) ...
        & hdrs.file_idx{adc_idx} <= good_file_idxs(2));
      hdrs.ver{adc_idx} = hdrs.ver{adc_idx}(good_idxs);
      hdrs.time{adc_idx} = hdrs.time{adc_idx}(good_idxs);
      hdrs.epri{adc_idx} = hdrs.epri{adc_idx}(good_idxs);
      hdrs.offset{adc_idx} = hdrs.offset{adc_idx}(good_idxs);
      hdrs.file_idx{adc_idx} = hdrs.file_idx{adc_idx}(good_idxs);
    end
  end
  
  if 0
    % Reselect a record range after you have already read in all the files
    
    % First specify range of file indexes to keep (inclusive/closed
    % range so the [2 4] includes files 2 and 4, use of +/- inf valid)
    param.records.file.rec_range = [1000 inf];
  end
  
  if 0
    %% Incomplete method for finding good indices
    plot(hdrs.epri{1},hdrs.time{1},'.')
    xlim([0 10e6]);
    ylim([0 1.4*86400]);
    
    [x,y] = ginput(2);
    hold on;
    plot([x(1) x(1) x(2) x(2) x(1)], [y(1) y(2) y(2) y(1) y(1)],'r');
    hold off;
    
    
    for adc_idx = 1:length(param.records.file.adcs)
      
      B = double(hdrs.epri{adc_idx});
      
      figure(1); clf;
      plot(B,hdrs.time{1},'.')
      xlim([0 10e6]);
      ylim([0 1.4*86400]);
      hold on;
      plot([x(1) x(1) x(2) x(2) x(1)], [y(1) y(2) y(2) y(1) y(1)],'r');
      hold off;
      
      good_idxs = find(B >= x(1) & B <= x(2) ...
        & hdrs.time{adc_idx} >= y(1) & hdrs.time{adc_idx} <= y(2));
      p = polyfit(B(good_idxs),hdrs.time{adc_idx}(good_idxs),1);
      EPRIs = round(x(1)) : round(x(2));
      utc_time_sod = polyval(p,EPRIs);
      hold on;
      plot(EPRIs,utc_time_sod,'r');
      plot(EPRIs,utc_time_sod+0.1,'k');
      plot(EPRIs,utc_time_sod-0.1,'k');
      hold off;
      
      time_expected = polyval(p,B(good_idxs));
      final_good_idxs = abs(hdrs.time{adc_idx}(good_idxs) - time_expected) < 0.1;
      final_good_idxs = good_idxs(final_good_idxs);
      
      plot(B(final_good_idxs), hdrs.time{adc_idx}(final_good_idxs))
      
      first_good_idx = find(final_good_idxs,1);
      hdrs.epri{adc_idx}(1:first_good_idx-1) ...
        = B(first_good_idx) + (-first_good_idx+1:-1);
      final_good_idxs(1:first_good_idx-1) = 1;
      
      % Correct EPRIs
      
      p = polyfit(1:length(final_good_idxs), B(final_good_idxs), 1);
      figure(1); clf;
      plot(1:length(final_good_idxs), B(final_good_idxs),'.')
      utc_time_sod = polyval(p,EPRIs);
      hold on;
      plot(EPRIs,utc_time_sod,'r');
      plot(EPRIs,utc_time_sod+0.1,'k');
      plot(EPRIs,utc_time_sod-0.1,'k');
      hold off;
      
      
      good_mask = zeros(size(B));
      good_mask(final_good_idxs) = 1;
      bad_mask = ~good_mask;
      bad_idxs = find(bad_mask);
      for bad_idx = bad_idxs
        B(bad_idx) = B(bad_idx-1) + 1;
      end
      
      keyboard
      %     hdrs.epri{adc_idx} = B;
    end
  end
  
end

%%
fprintf('Fix all the EPRIs\n');
for adc_idx = 1:length(param.records.file.adcs)
  B = double(hdrs.epri{adc_idx});
  
  done = false;
  fprintf('Fixing %d\n', adc_idx);
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
  
  hdrs.epri{adc_idx} = uint32(B);
  [hdrs.epri{adc_idx} unique_idxs] = unique(hdrs.epri{adc_idx});
  hdrs.seconds{adc_idx} = hdrs.seconds{adc_idx}(unique_idxs);
  hdrs.fractions{adc_idx} = hdrs.fractions{adc_idx}(unique_idxs);
  hdrs.file_idx{adc_idx} = hdrs.file_idx{adc_idx}(unique_idxs);
  hdrs.offset{adc_idx} = hdrs.offset{adc_idx}(unique_idxs);
end

if 0
  % For debugging
  figure(1); clf;
  for adc_idx = 1:length(param.records.file.adcs)
    title(sprintf('%d: %d elements',adc_idx, length(hdrs.epri{adc_idx})))
    plot(hdrs.epri{adc_idx})
    pause;
    hold on
    plot(hdrs.epri{adc_idx},'r')
  end
  keyboard
end

%% Synchronize ADC channels
% ===================================================================

% ===================================================================
% Find the first and last good EPRI that is available in all channels
%   - This is a robust way of finding this even when the first few
%     records are corrupt
%   - It works well as long as there is not a dropped record, even
%     then it performs okay
% ===================================================================
first_EPRIs = [];
last_EPRIs = [];
for adc_idx = 1:length(param.records.file.adcs)
  N = 21;
  actual_epri = double(hdrs.epri{adc_idx}(1:N));
  for offset = 1:N
    ideal_EPRIs = actual_epri(offset)-(offset-1) + (0:N-1);
    metric(offset) = sum(actual_epri(1:N) == ideal_EPRIs);
  end
  [max_metric best_offset] = max(metric);

  first_EPRIs(adc_idx) = actual_epri(best_offset) - (best_offset-1);

  N = 21;
  actual_epri = double(hdrs.epri{adc_idx}(end-N+1:end));
  for offset = 1:N
    ideal_EPRIs = actual_epri(offset)-(offset-1) + (0:N-1);
    metric(offset) = sum(actual_epri(1:N) == ideal_EPRIs);
  end
  [max_metric best_offset] = max(metric);

  last_EPRIs(adc_idx) = actual_epri(best_offset) - (best_offset-1);
end

first_EPRI = max(first_EPRIs);
last_EPRI = min(last_EPRIs);

% ===================================================================
% Check each record to see if there is a good EPRI and good time.
%  - Store status in adc_val
%    1: EPRI bad
%    2: time bad
%  - This will aid in interpolation later
% ===================================================================
fprintf('Format/sync records/ADCs\n');
confirmed_good_first_EPRI = false;
while ~confirmed_good_first_EPRI
  confirmed_good_first_EPRI = true;
  for adc_idx = 1:length(param.records.file.adcs)
    possible_offset = find(double(hdrs.epri{adc_idx}) == first_EPRI,1) - 1;
    if ~isempty(possible_offset)
      offset(adc_idx) = possible_offset;
    else
      possible_offset = find(double(hdrs.epri{adc_idx}) >= first_EPRI,1) - 1;
      offset(adc_idx) = possible_offset;
      if offset(adc_idx) ~= first_EPRI - first_EPRIs(adc_idx)
        fprintf('Predicted offset does not make sense... there is probably bad header fields\n');
        if param.debug_level >= 2
          fprintf('Type: plot(hdrs.epri{adc_idx}(1:%d))\n', max(offset(adc_idx), first_EPRI - first_EPRIs(adc_idx)));
          keyboard
        end
      end
      fprintf('  Reassigning first_EPRI\n');
      first_EPRI = hdrs.epri{adc_idx}(find(double(hdrs.epri{adc_idx}) >= first_EPRI,1));
      confirmed_good_first_EPRI = false;
    end
    first_times(adc_idx) = hdrs.time{adc_idx}(1+offset(adc_idx));
  end
end

cur_time = median(first_times());
dt = num_presum/param.radar.prf;

EPRIs = first_EPRI:last_EPRI;

%%
fprintf('Fix all the UTC times\n');
for adc_idx = 1:length(param.records.file.adcs)
  B = double(hdrs.time{adc_idx});
  
  done = false;
  fprintf('Fixing %d\n', adc_idx);
  iter=1;
  while ~done
    fprintf('  Iteration %d\n', iter);
    A = medfilt1(B,11);
    A(1:5) = median(B(1:6));
    A(end-4:end) = median(B(end-5:end));
    bad_mask = abs(A-B) > 50*dt;
    if bad_mask(1)
      first_good_idx = find(~bad_mask,1);
      bad_mask(1:first_good_idx-1) = 0;
      B(1:first_good_idx-1) = B(first_good_idx) + (-first_good_idx+1:-1)*dt;
    end
    bad_idxs = find(bad_mask);
    good_idxs = find(bad_mask)-1;
    B(bad_idxs) = B(good_idxs) + dt;
    if isempty(bad_idxs)
      done = true;
    end
    iter = iter+1;
  end
  
  hdrs.time{adc_idx} = B;
end

if 0
  % For debugging
  figure(2); clf;
  for adc_idx = 1:length(param.records.file.adcs)
    title(sprintf('%d: %d elements',adc_idx, length(hdrs.time{adc_idx})))
    plot(hdrs.time{adc_idx})
    pause;
    hold on
    plot(hdrs.time{adc_idx},'r')
    ylim([5e4 8e4])
  end
  keyboard
end

%% Synchronize each channel
adc_valid = zeros(size(EPRIs));
good_mask = logical(zeros(length(param.records.file.adcs),length(EPRIs)));
final_times = NaN*zeros(length(param.records.file.adcs),length(EPRIs));
for adc_idx = 1:length(param.records.file.adcs)
  valid_idxs = find(hdrs.epri{adc_idx} >= EPRIs(1) & hdrs.epri{adc_idx} <= EPRIs(end));
  good_mask(adc_idx, hdrs.epri{adc_idx}(valid_idxs) - EPRIs(1)+1) = 1;
  final_times(adc_idx, good_mask(adc_idx,:)) = hdrs.time{adc_idx}(valid_idxs);
end
bad_mask = sum(good_mask) < 1;
best_time = nanmedian(final_times);
best_time = best_time(~bad_mask);

%% Guarantee monotonically increasing time
old_time_idx = 1;
time_idx = 2;
time_mask = ones(size(best_time));
while time_idx < length(best_time)
  while best_time(time_idx) <= best_time(old_time_idx)
    time_mask(time_idx) = 0;
    time_idx = time_idx + 1;
  end
  time_idx = time_idx + 1;
end
if any(time_mask == 0)
  warning('Time is not monotonically increasing');
  keyboard
end

% ===================================================================
%% Correlate GPS with radar data
% ===================================================================
records = [];
records.raw.epri = hdrs.epri;
records.raw.seconds = hdrs.seconds;
records.raw.fractions = hdrs.fractions;

% ===================================================================
% Arrange the output variables and synchronize with the GPS
adc_idx = 1;
if param.records.gps.en
  records = sync_radar_to_gps(param,records,best_time);
  
else
  records.lat = NaN*zeros(1,sum(~bad_mask));
  records.lon = NaN*zeros(1,sum(~bad_mask));
  records.elev = NaN*zeros(1,sum(~bad_mask));
  records.gps_time = NaN*zeros(1,sum(~bad_mask));
  records.roll = NaN*zeros(1,sum(~bad_mask));
  records.pitch = NaN*zeros(1,sum(~bad_mask));
  records.heading = NaN*zeros(1,sum(~bad_mask));
  records.gps_source = 'NA';
end

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
for adc_idx = 1:length(param.records.file.adcs)
  valid_idxs = find(hdrs.epri{adc_idx} >= EPRIs(1) & hdrs.epri{adc_idx} <= EPRIs(end));
  unique_file_idxs = unique(hdrs.file_idx{adc_idx}(valid_idxs));
  for file_idx = 1:length(unique_file_idxs)
    first_valid_idx = find(hdrs.file_idx{adc_idx}(valid_idxs) == unique_file_idxs(file_idx),1);
    final_idx = hdrs.epri{adc_idx}(valid_idxs(first_valid_idx)) - EPRIs(1)+1;
    final_idx = final_idx - sum(bad_mask(1:final_idx-1));
    records.relative_rec_num{adc_idx}(file_idx) = final_idx;
    [fn_dir fn_name fn_ext] = fileparts(filenames{adc_idx}{unique_file_idxs(file_idx)});
    records.relative_filename{adc_idx}{file_idx} = [fn_name fn_ext];
  end
  records.offset(adc_idx,1:length(EPRIs)) = -2^31;
  records.offset(adc_idx,good_mask(adc_idx,:)) = hdrs.offset{adc_idx}(valid_idxs);
end
records.offset = records.offset(:,~bad_mask);
records.radar_name = param.radar_name;
records.ver = 3;

records.settings = [];

records.settings.wfs_records = 1;
wfs{1}.bit_shifts = wfs{1}.which_bits;
wfs{1} = rmfield(wfs{1},'which_bits');
for wf = 1:length(wfs{1}.num_sam)
  records.settings.wfs(wf).num_sam = wfs{1}.num_sam(wf);
  records.settings.wfs(wf).t0 = wfs{1}.t0(wf);
  records.settings.wfs(wf).presums = wfs{1}.presums(wf);
  records.settings.wfs(wf).bit_shifts = wfs{1}.bit_shifts(wf);
end

records.notes = '';
records.param_records = param;

records_fn = ct_filename_support(param,'','records');
fprintf('Saving records file %s (%s)\n',records_fn,datestr(now));
save(records_fn,'-v7.3','-struct','records'); % Handle large file sizes, so use v7.3

% =====================================================================
% Create record aux files for faster loading times
% =====================================================================

fprintf('Creating auxiliary records files %s (%s)\n',records_fn,datestr(now));
create_records_aux_files(records_fn);

fprintf('Done (%s)\n\n', datestr(now));

return;
