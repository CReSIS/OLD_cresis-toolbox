function hdrs = create_records_accum2_task(filenames,file_idxs);
% hdrs = create_records_accum2_task(filenames,file_idxs);
%
% This function is called by create_records_accum2.m.
%
% filenames = cell vector containing all the accum2 files that CAN be
%   loaded
% file_idxs = number vector containing specific indices that will be
%   loaded.
%
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
%
% Author: Anthony Hoch, John Paden, Logan Smith

tstart_create_records_accum2_task = tic;

% =====================================================================
% Get the headers from all the files
% =====================================================================

for file_idxs_idx = 1:length(file_idxs)
  file_idx = file_idxs(file_idxs_idx);
  fn = filenames{file_idx};
  [tmp fn_name fn_ext] = fileparts(fn);
  fprintf('  File %s (%.1f sec)\n', fn_name, ...
    toc(tstart_create_records_accum2_task));
  
  hdr = basic_load_accum2(fn,struct('clk',1e9,'all_hdrs',true));

  if file_idxs_idx == 1
    hdrs.wfs_file = 1;
    hdrs.wfs{file_idxs_idx} = hdr.wfs;
    hdrs.filenames{1}{file_idxs_idx} = [fn_name fn_ext];
    hdrs.file_rec_offset = 1;
    hdrs.comp_time = hdr.comp_time;
    hdrs.radar_time = hdr.radar_time;
    hdrs.radar_time_1pps = hdr.radar_time_1pps;
    hdrs.range_gate_start = hdr.range_gate_start;
    hdrs.range_gate_duration = hdr.range_gate_duration;
    hdrs.trigger_delay = hdr.trigger_delay;
    hdrs.num_coh_ave = hdr.num_coh_ave;
    % File index for all of these records will be the same
    hdrs.file_idx{1} = file_idxs_idx*ones(size(hdr.radar_time));
    % Offset is the byte offset in the file for each record (zero-indexed)
    if isfield(hdr.finfo,'bad_mask')
      % File was corrupted and frame_sync_info had to be used to find
      % the start of each record, we will use that... Notify user since
      % this should never happen (should just be able to type "dbcont")
      keyboard
      hdrs.offset{1} = hdr.finfo.syncs;
    else
      % File format is good, so we assume regular record spacing
      hdrs.offset{1} = hdr.finfo.syncs(1) + hdr.finfo.rec_size*(0:length(hdr.frame_sync)-1);
    end
    
  else
    % File index for all of these records will be the same
    hdrs.file_idx{1} = cat(2,hdrs.file_idx{1},file_idxs_idx*ones(size(hdr.radar_time)));
    % Offset is the byte offset in the file for each record (zero-indexed)
    if isfield(hdr.finfo,'bad_mask')
      warning('File was corrupted and frame_sync_info had to be used to find the start of each record, we will use that... Notifying user since this should never happen unless there was a write failure (should just be able to type "dbcont" to continue)');
      keyboard
      hdrs.offset{1} = cat(2,hdrs.offset{1},hdr.finfo.syncs);
    else
      % File format is good, so we assume regular record spacing
      hdrs.offset{1} = cat(2,hdrs.offset{1},hdr.finfo.syncs(1) + hdr.finfo.rec_size*(0:length(hdr.frame_sync)-1));
    end
    hdrs.filenames{1}{file_idxs_idx} = [fn_name fn_ext];
    hdrs.file_rec_offset(file_idxs_idx) = length(hdrs.comp_time) + 1;
    hdrs.comp_time = cat(2,hdrs.comp_time,hdr.comp_time);
    hdrs.radar_time = cat(2,hdrs.radar_time,hdr.radar_time);
    hdrs.radar_time_1pps = cat(2,hdrs.radar_time_1pps,hdr.radar_time_1pps);
    hdrs.range_gate_start = cat(2,hdrs.range_gate_start,hdr.range_gate_start);
    hdrs.range_gate_duration = cat(2,hdrs.range_gate_duration,hdr.range_gate_duration);
    hdrs.num_coh_ave = cat(2,hdrs.num_coh_ave,hdr.num_coh_ave);
%     store_header = false;
%     if hdr.IndexHeader ~= hdrs.wfs(end).IndexHeader ...
%         || hdr.SampleFrequency ~= hdrs.wfs(end).SampleFrequency ...
%         || hdr.PRF ~= hdrs.wfs(end).PRF ...
%         || hdr.PreSum ~= hdrs.wfs(end).PreSum ...
%         || hdr.NumberWaveforms ~= hdrs.wfs(end).NumberWaveforms ...
%         || any(hdr.IndexRecordStart(:) ~= hdrs.wfs(end).IndexRecordStart(:)) ...
%         || any(hdr.IndexRecordStop(:) ~= hdrs.wfs(end).IndexRecordStop(:)) ...
%         || hdr.IndexData ~= hdrs.wfs(end).IndexData
%       store_header = true;
%     end
%     for wf = 1:length(hdr.Waveform)
%       if hdr.Waveform(wf).StartFrequency ~= hdrs.wfs(end).Waveform(wf).StartFrequency ...
%           || hdr.Waveform(wf).StopFrequency ~= hdrs.wfs(end).Waveform(wf).StopFrequency ...
%           || hdr.Waveform(wf).PulseDuration ~= hdrs.wfs(end).Waveform(wf).PulseDuration ...
%           || hdr.Waveform(wf).ZeroPiModulation ~= hdrs.wfs(end).Waveform(wf).ZeroPiModulation ...
%           || hdr.Waveform(wf).TxMultiplexer ~= hdrs.wfs(end).Waveform(wf).TxMultiplexer ...
%           || any(hdr.Waveform(wf).TxAmpEnable(:) ~= hdrs.wfs(end).Waveform(wf).TxAmpEnable(:)) ...
%           || hdr.Waveform(wf).ModCount0 ~= hdrs.wfs(end).Waveform(wf).ModCount0 ...
%           || hdr.Waveform(wf).ModCount1 ~= hdrs.wfs(end).Waveform(wf).ModCount1 ...
%           || any(hdr.Waveform(wf).NumberSamples(:) ~= hdrs.wfs(end).Waveform(wf).NumberSamples(:)) ...
%           || any(hdr.Waveform(wf).SampleDelay(:) ~= hdrs.wfs(end).Waveform(wf).SampleDelay(:)) ...
%           || any(hdr.Waveform(wf).RecordEnable(:) ~= hdrs.wfs(end).Waveform(wf).RecordEnable(:)) ...
%           || any(hdr.Waveform(wf).BlankDelay(:) ~= hdrs.wfs(end).Waveform(wf).BlankDelay(:)) ...
%           || hdr.Waveform(wf).RxAttenuation ~= hdrs.wfs(end).Waveform(wf).RxAttenuation
%         store_header = true;
%       end
%     end
%     if store_header
%       warning('Header change');
%       hdrs.wfs_file(end+1) = file_idxs_idx;
%       hdrs.wfs(end+1) = hdr;
%     end
  end
  
%   tmp_hdr.offset      = Header.IndexData + 4+4+4+8+4 + Header.IndexRecordStart(adc_idx,1) - 1 + ...
%                         (0:(Header.NumberRecords-1)) .* [4+4+4+8+4+2*Header.IndexRecordStop(end,end)];                  
%   tmp_hdr.offset = tmp_hdr.offset + Header.IndexRecordStart(adc_idx,1) - 1;
     
%   for idx_wfs = 1:Header.NumberWaveforms
%       wfs.num_sam(idx_wfs)      = Header.Waveform(idx_wfs).NumberSamples(1);
%       wfs.which_bits(idx_wfs)   = 1;
%       wfs.samp_delay(idx_wfs)   = Header.Waveform(idx_wfs).SampleDelay(1);
%       wfs.presums(idx_wfs)      = Header.PreSum;
%   end
                    
end

return;
