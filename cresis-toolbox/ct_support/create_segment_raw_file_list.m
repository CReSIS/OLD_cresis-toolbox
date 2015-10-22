% script create_segment_raw_file_list
%
% Helps find the file_prefix, start_idx, and stop_idx fields for
% creating segments in the parameter spreadsheet. This function will
% be typically run at the end of a mission when the segments are being
% added to the spreadsheet for the first time.
%
% Author: John Paden

% =========================================================================
%% User Settings
% =========================================================================
base_path = '/cresis/snfs1/data/Accum_Data/2013_Antarctica_Ground/20140110/';
radar_name = 'accum2';

% base_path = '/process/fmcw/snow/'; % USE create_segment_raw_file_list_v2.m
% radar_name = 'snow3';

% base_path = '/process/fmcw/kuband/'; % USE create_segment_raw_file_list_v2.m
% radar_name = 'kuband3';

% base_path = '/process/mcords/';
% radar_name = 'mcords3';

% base_path = '/cresis/data1/MCRDS/2009_Greenland/20090401/dataGISMO/';
% radar_name = 'mcrds';

min_seg_size = 2;

% Breaking of segments based on file name times uses the median time difference
% between file name time stamps and then adds a guard... sometimes the 
% file writes can get behind schedule/have a lot of variance with a single
% segment and larger values should be used. 0.5 is good when this is not
% the case.
% TIME_SEGMENT_GUARD = 0.5;
TIME_SEGMENT_GUARD = 1.5;

% Choose one method for looking at segment breaks:
% allbig: All big files sizes, defined as 90% of the max file size and
%   occuring at least twice, are considered full files and not a segment
%   break. This is the preferred method for FMCW radars with DDC because
%   their file sizes change often. However, create_segment_raw_file_list_v2.m
%   is now the preferred function to use for these systems.
% file_size_checking = 'allbig';
% anydifference: Any difference in the file size will cause a segment
%    break to happen. This is the preferred method for all other systems.
file_size_checking = 'anydifference';

% =========================================================================
%% Automated Section
% =========================================================================

% Get list of filenames
if any(strcmpi(radar_name,{'mcrds'}))
  raw_file_prefix = 'data';
  raw_file_suffix = 'raw';
else
  raw_file_prefix = radar_name;
  raw_file_suffix = 'bin';
end
fns = get_filenames(base_path,raw_file_prefix,'',raw_file_suffix);

if isempty(fns)
  error('No files matching %s/%s*%s', base_path, raw_file_prefix, raw_file_suffix);
end

% Get basic information from each filename (date and file prefix)
file_prefix_list = {};
file_dates = zeros(size(fns));
fn_size = zeros(size(fns));
for fn_idx = 1:length(fns)
  fn = fns{fn_idx};
  if any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3'}))
    fname_info = fname_info_fmcw(fn);
  elseif any(strcmpi(radar_name,{'mcrds'}))
    fname_info = fname_info_mcrds(fn);
  elseif any(strcmpi(radar_name,{'mcords'}))
    fname_info = fname_info_mcords(fn);
  elseif any(strcmpi(radar_name,{'mcords2','mcords3','mcords4'}))
    fname_info = fname_info_mcords2(fn);
  elseif any(strcmpi(radar_name,{'accum'}))
    fname_info = fname_info_accum(fn);
  elseif any(strcmpi(radar_name,{'accum2'}))
    fname_info = fname_info_accum2(fn);
  end
  file_dates(fn_idx) = fname_info.datenum;
  [fn_dir,fn_name] = fileparts(fn);
  if any(strcmpi(radar_name,{'mcrds'}))
    file_prefix = fn_name(1:find(fn_name=='.',1)+2);
  else
    file_prefix = fn_name(1:find(fn_name=='_',1)+2);
  end
  if ~any(strcmpi(file_prefix,file_prefix_list))
    file_prefix_list = cat(1,file_prefix_list,file_prefix);
  end
  fn_finfo = dir(fn);
  fn_size(fn_idx) = fn_finfo.bytes;
end

% Sort the files by date (segments must be numbered according to the time
% collected
[file_dates sort_idxs] = sort(file_dates);
fns = fns(sort_idxs);
fn_size = fn_size(sort_idxs);

% Plot the file dates
figure(1); clf;
plot(file_dates)

% Break up the files according to gaps in time between each raw file
segment_idxs = [1];

if strcmpi(file_size_checking,'allbig')
  % One way that segment boundaries are found is to look at files that are
  % no "full".  Since there may be more than one "full" size (due to
  % different record lengths), we look through all the file sizes and
  % find the ones that are most plausible to be the full sizes.
  % First, we assume that full sizes happen often, so we remove any file
  % size that is unique
  [~,unique_idxs] = unique(fn_size);
  not_unique_mask = logical(ones(size(fn_size)));
  not_unique_mask(unique_idxs) = 0;
  full_file_sizes = unique(fn_size(not_unique_mask));
  % Then, we assume that all the full sizes are pretty close to the max.
  full_file_sizes = full_file_sizes(full_file_sizes > 0.9*max(full_file_sizes));
  full_file_sizes
end

if 0
  % Segment breaks happen when max file size is not met...
  segment_idxs = [1; 1+find(fn_size < max_file_size)];
  segment_idxs = segment_idxs(1:end-1);
elseif 0
  segment_idxs = [1; 1+find(abs(diff(file_dates)) > 1.35*median(diff(file_dates)))];
else
  file_dates = file_dates * 86400;
  segment_idxs = [1];
  med_diff = median(diff(file_dates(1:4)));
  for idx = 2:length(file_dates)-1
    if strcmpi(file_size_checking,'allbig') & all(fn_size(idx-1) ~= full_file_sizes)
      fprintf('File size segment break %d\n  %s %i\n  %s %i\n  %s %i\n', ...
        idx, fns{idx-1}, fn_size(idx-1), fns{idx}, fn_size(idx), fns{idx+1}, fn_size(idx+1));
      segment_idxs = cat(1,segment_idxs,idx);
      med_diff = median(diff(file_dates(idx:idx+min(3,length(file_dates)-idx))));
    elseif strcmpi(file_size_checking,'anydifference') ...
        && (idx >= 3 && fn_size(idx-1) ~= fn_size(idx-2) && fn_size(idx) ~= fn_size(idx-1) ...
        || idx == 2 && fn_size(idx) ~= fn_size(idx-1))
      segment_idxs = cat(1,segment_idxs,idx);
      med_diff = median(diff(file_dates(idx:idx+min(3,length(file_dates)-idx))));
    elseif abs(file_dates(idx)-file_dates(idx-1) - med_diff) > TIME_SEGMENT_GUARD*med_diff
      if ~(any(strcmpi(radar_name,{'mcrds'})) && ~isempty(segment_idxs) && segment_idxs(end) == idx-1)
        fprintf('Time segment break %d\n  %s %i\n  %s %i\n  %s %i\n', ...
          idx, fns{idx-1}, fn_size(idx-1), fns{idx}, fn_size(idx), fns{idx+1}, fn_size(idx+1));
        segment_idxs = cat(1,segment_idxs,idx);
        med_diff = median(diff(file_dates(idx:idx+min(3,length(file_dates)-idx))));
      end
    end
  end
end
fprintf('\n');

% Print out segment list for parameter spreadsheet
fprintf('FORMAT: segment Segment file_prefix start_idx to stop_idx\n\n');
file_prefix_list_idxs = zeros(size(file_prefix_list));
segments = [];
for idx = 1:length(segment_idxs)
  start_fn = fns{segment_idxs(idx)};
  [start_fn_dir start_fn_name] = fileparts(start_fn);
  
  % Find the last file index associated with this filename in the sorted
  % list
  if idx < length(segment_idxs)
    end_idx = segment_idxs(idx+1)-1;
  else
    end_idx = length(fns);
  end
  
  % Determine offset indices to plug into parameter spreadsheet
  if any(strcmpi(radar_name,{'mcrds'}))
    file_prefix = start_fn_name(1:find(start_fn_name=='.',1)+2);
  else
    file_prefix = start_fn_name(1:find(start_fn_name=='_',1)+2);
  end  
  
  file_prefix_idx = strmatch(file_prefix,file_prefix_list);
  start_idx = file_prefix_list_idxs(file_prefix_idx) + 1;
  num_files = end_idx - segment_idxs(idx) + 1;
  stop_idx = start_idx + num_files - 1;
  file_prefix_list_idxs(file_prefix_idx) ...
    = file_prefix_list_idxs(file_prefix_idx) + num_files;

  if stop_idx - start_idx < min_seg_size
    % Print out this segment information
    fprintf('Segment %i %s %4i to %4i: TOO SHORT\n  %s\n  %s\n', ...
      idx, file_prefix, start_idx, stop_idx, start_fn, fns{end_idx});
  else
    % Print out this segment information
    fprintf('Segment %i %s %4i to %4i:\n  %s\n  %s\n', ...
      idx, file_prefix, start_idx, stop_idx, start_fn, fns{end_idx});
  end
  
  segments(idx).start_idx = start_idx;
  segments(idx).stop_idx = stop_idx;
  segments(idx).file_prefix = file_prefix;

end

% Print out some results that can be copied and pasted easily
seg_id = 0;
for idx = 1:length(segments)
  if segments(idx).stop_idx - segments(idx).start_idx >= min_seg_size
    seg_id = seg_id + 1;
    fprintf('%s\t%02d\t%d\t%d\t%s\n', datestr(file_dates(1)/86400,'yyyymmdd'), idx, segments(idx).start_idx, segments(idx).stop_idx, segments(idx).file_prefix);
  end
end

return;


