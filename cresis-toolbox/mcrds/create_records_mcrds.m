function create_records_mcrds(param, param_override)
% create_records_mcrds(param, param_override)
%
% Function for creating records file for MCORDS data. This function can
% be called as a script by:
%  1. Commenting out the function line
%  2. Setting the default param structure
%  3. Uncommenting param = [] line
% This is useful for debugging.
%
% This function should be run after the GPS file has been created.
% For example, cresis-toobox/gps/missions/make_gps_2009_antarctica_DC8.m
%
% This function's output is used by create_frames.m.
%
% Author: John Paden
%
% See also: check_records_mcrds, create_records_mcrds_task,create_records_mcrds_post_sync

% =====================================================================
% General Setup
% =====================================================================

dbstack_info = dbstack;
if ~exist('param','var') || isempty(param) || length(dbstack_info) == 1
  % =====================================================================
  % Debug Setup
  % =====================================================================
  param = read_param_xls(ct_filename_param('rds_param_2006_Greenland_TO.xls'),'20060526_01');
  
  clear('param_override');
  param_override.sched.type = 'no scheduler';
  param_override.sched.rerun_only = true;

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
  
elseif ~isstruct(param)
  % Functional form
  param();
end
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', dbstack_info(1).name, param.day_seg, datestr(now,'HH:MM:SS'));
fprintf('=====================================================================\n');

% =====================================================================
% Prep work
% =====================================================================

% Get WGS84 ellipsoid parameters
physical_constants

% =====================================================================
%% Get the headers from all the files
% =====================================================================

% Get the list of files to include in this records file
fprintf('Getting files (%s)\n', datestr(now));
fns = {};
[base_dir,adc_folder_name,fns{1},file_idxs] = get_segment_file_list(param);

% hdrs = struct with header information from each file (this is the
%  foundation of the "records" variable)
%  .wfs = a struct vector of headers (contains complete header information
%    loaded from MCRDS). If the settings do not change in this segment
%    then wfs is length one. Each time the settings change, the header from
%    the first file for which the settings have been changed is stored in
%    this list.
%  .wfs_file = a number vector corresponding to .wfs which tells the
%    specific file that each wfs was loaded from.
%  .filenames = cell vector of strings containing filenames, same length
%    as file_idxs
%  .file_rec_offset = number vector with equal length to .filenames which
%    tells which record this file starts with
%  .comp_time = all the computer times from all the files in chronological
%    order (one computer time per record)
%  .radar_time = all the radar times from all the files in chronological
%    order (one radar time per record)
for file_idxs_idx = 1:length(file_idxs)
  file_idx = file_idxs(file_idxs_idx);
  fn = fns{1}{file_idx};
  [tmp fn_name fn_ext] = fileparts(fn);
  fprintf('  File %s %d of %d (%s)\n', fn_name, file_idxs_idx, length(file_idxs), datestr(now));
  
  hdr = basic_load_mcrds_hdr(fn);
  [comp_time radar_time] = basic_load_mcrds_time(fn,hdr);

  if file_idxs_idx == 1
    hdrs.wfs_file = 1;
    hdrs.wfs = hdr;
    hdrs.filenames{1} = [fn_name fn_ext];
    hdrs.file_rec_offset = 1;
    hdrs.comp_time = comp_time;
    hdrs.radar_time = radar_time;
  else
    hdrs.filenames{file_idxs_idx} = [fn_name fn_ext];
    hdrs.file_rec_offset(file_idxs_idx) = length(hdrs.comp_time) + 1;
    hdrs.comp_time = cat(2,hdrs.comp_time,comp_time);
    hdrs.radar_time = cat(2,hdrs.radar_time,radar_time);
    store_header = false;
    if hdr.IndexHeader ~= hdrs.wfs(end).IndexHeader ...
        || hdr.SampleFrequency ~= hdrs.wfs(end).SampleFrequency ...
        || hdr.PRF ~= hdrs.wfs(end).PRF ...
        || hdr.PreSum ~= hdrs.wfs(end).PreSum ...
        || hdr.NumberWaveforms ~= hdrs.wfs(end).NumberWaveforms ...
        || any(hdr.IndexRecordStart(:) ~= hdrs.wfs(end).IndexRecordStart(:)) ...
        || any(hdr.IndexRecordStop(:) ~= hdrs.wfs(end).IndexRecordStop(:)) ...
        || hdr.IndexData ~= hdrs.wfs(end).IndexData
      store_header = true;
    end
    for wf = 1:length(hdr.Waveform)
      if hdr.Waveform(wf).StartFrequency ~= hdrs.wfs(end).Waveform(wf).StartFrequency ...
          || hdr.Waveform(wf).StopFrequency ~= hdrs.wfs(end).Waveform(wf).StopFrequency ...
          || hdr.Waveform(wf).PulseDuration ~= hdrs.wfs(end).Waveform(wf).PulseDuration ...
          || hdr.Waveform(wf).ZeroPiModulation ~= hdrs.wfs(end).Waveform(wf).ZeroPiModulation ...
          || hdr.Waveform(wf).TxMultiplexer ~= hdrs.wfs(end).Waveform(wf).TxMultiplexer ...
          || any(hdr.Waveform(wf).TxAmpEnable(:) ~= hdrs.wfs(end).Waveform(wf).TxAmpEnable(:)) ...
          || hdr.Waveform(wf).ModCount0 ~= hdrs.wfs(end).Waveform(wf).ModCount0 ...
          || hdr.Waveform(wf).ModCount1 ~= hdrs.wfs(end).Waveform(wf).ModCount1 ...
          || any(hdr.Waveform(wf).NumberSamples(:) ~= hdrs.wfs(end).Waveform(wf).NumberSamples(:)) ...
          || any(hdr.Waveform(wf).SampleDelay(:) ~= hdrs.wfs(end).Waveform(wf).SampleDelay(:)) ...
          || any(hdr.Waveform(wf).RecordEnable(:) ~= hdrs.wfs(end).Waveform(wf).RecordEnable(:)) ...
          || any(hdr.Waveform(wf).BlankDelay(:) ~= hdrs.wfs(end).Waveform(wf).BlankDelay(:)) ...
          || hdr.Waveform(wf).RxAttenuation ~= hdrs.wfs(end).Waveform(wf).RxAttenuation
        store_header = true;
      end
    end
    if store_header
      warning('Header change in %d/%d', file_idxs_idx, file_idx);
      hdrs.wfs_file(end+1) = file_idxs_idx;
      hdrs.wfs(end+1) = hdr;
    end
  end
      
end

% Count the presums in one EPRI (assumes number of waveforms
% and number of presums does not change)
num_presum = hdrs.wfs(1).PreSum * hdrs.wfs(1).NumberWaveforms;

% ===================================================================
% ===================================================================
% Synchronize radar data to GPS data
% ===================================================================
% ===================================================================
if isfield(param,'tmp_path') && ~isempty(param.tmp_path)
  fn = ct_filename_tmp(param,param.records.records_fn,'records','workspace');
  fprintf('Saving workspace %s (%s)\n', fn, datestr(now));
  try
    fn_dir = fileparts(fn);
    if ~exist(fn_dir,'dir')
      mkdir(fn_dir);
    end
    save(fn);
  catch
    fprintf('  Saving workspace failed\n');
  end
end

create_records_mcrds_sync;

return;
