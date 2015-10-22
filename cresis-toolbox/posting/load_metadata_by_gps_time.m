function [metadata] = load_metadata_by_gps_time(param,gps_times)
% [metadata] = load_metadata_by_gps_time(param,gps_times)
%
% Loads metadata from one of the radars (e.g. rds) and returns the metadata
% associated with the specified gps times.
%
% param = parameters structure
%  .radar_name = e.g. 'accum', 'snow', etc
%  .season_name = e.g. '2009_Greenland_P3'
%  .day_seg = optional 'YYYYMMDD_SS' string to restrict the search
% gps_times: double vector of gps times
%
% metadata =
%   .day_seg: cell vector of frame name strings 'YYYYMMDD_SS'
%   .frm: double vector of frames
%   .record: double vector of records
%   .along_track: vector of along track position within that frame
%
% Example:
%   load_metadata_by_gps_time(struct('radar_name','snow','season_name','2009_Greenland_P3'))
%
% Author: John Paden

physical_constants;

% Get list of records files
if isfield(param,'day_seg')
  % Grab a specific segment
  fns = {ct_filename_support(param, '', 'records')};
else
  % Grab all segments
  records_dir = ct_filename_support(param, '', 'records');
  fns = get_filenames(records_dir,'records','','.nc');
end

if isempty(fns)
  error('No records file found in %s', records_dir);
end

% Load in first and last record from each records file
first_gps_time = [];
for file_idx = 1:length(fns)
  records_fn = fns{file_idx};
  [path name] = fileparts(records_fn);
  cdf_fn = fullfile(path, sprintf('%s.nc', name));
  
  try
    ncid = netcdf.open(cdf_fn,'NOWRITE');
  catch ME
    warning('Exception during file opening');
    ME
    keyboard
  end
  var_idx = netcdf.inqVarID(ncid,'gps_time');
  first_gps_time(file_idx) = netcdf.getVar(ncid,var_idx,[0 0],[1 1]);
  netcdf.close(ncid);
end

metadata = [];
for tidx = 1:length(gps_times)
  gps_time = gps_times(tidx);
  
  start_file_idx = find(gps_time >= first_gps_time,1,'last');
  stop_file_idx = find(gps_time >= first_gps_time,1,'last');
  
  if isempty(start_file_idx)
    error('Start time is before any radar data exists');
  end
  if isempty(stop_file_idx)
    error('Stop time is after any radar data exists');
  end
  if start_file_idx ~= stop_file_idx
    error('Time range includes times when no radar data exists');
  end
  
  records_fn = [fns{start_file_idx}(1:end-3) 'mat'];
  [tmp records_fn_name] = fileparts(records_fn);
  param.day_seg = records_fn_name(9:end);
  
  frames_fn = ct_filename_support(param, '', 'frames');
  load(frames_fn);
  
  records_ver = load(records_fn,'ver');
  if isfield(records_ver,'ver')
    records = load(records_fn);
  else
    load(records_fn, 'records');
  end
  
  [~,metadata.record(tidx)] = min(abs(records.gps_time - gps_time));
  
  % Determine which frames contain these records
  metadata.frame(tidx) = find(frames.frame_idxs <= metadata.record(tidx),1,'last');
  
  records_list = frames.frame_idxs(metadata.frame(tidx)) : metadata.record(tidx);
  along_track = geodetic_to_along_track(records.lat(records_list),records.lon(records_list));
  metadata.along_track(tidx) = along_track(end);
  
  metadata.day_seg{tidx} = records.param_records.day_seg;

end

return;
