function [success hdrs wfs] = create_records_acords_task(filenames,param,file_idxs,adc_idx)
% [success hdrs wfs] = create_records_acords_task(filenames,param,file_idxs,adc_idx)
%
% This function is called by create_records_mcrds.m.
% It has been created to allow the cluster to be used.

% success = boolean which is true when the function executes properly
%   if a task fails before it can return success, then success will be
%   empty
% hdrs = struct with header information from each file
%  FIELDS?
% wfs = struct with waveform information from the first header
%  FIELDS?

% When called from scheduler, we need to run 
try; toc; catch; tic; end


% =====================================================================
fn = filenames{file_idxs(1)};   % HAVE THIS COLLECT ALL RX INDEXES, THEN ONLY KEEP THOSE THAT MATTER OUTSIDE OF THE TASK (WASTEFUL, BUT MAINTAINS METHOD)
fname.name           = 'acords';
fname.group         = '';
fname.file_idx      = file_idxs(1);
fname.radar_num     = 1;
fname.rx            = 0;    % file contains all, this task returns all
time_stamp          = fn(end-22:end-9);
fname.datenum       = datenum(  str2double(time_stamp(1:4)), ...
                                str2double(time_stamp(5:6)), str2double(time_stamp(7:8)), ...
                                str2double(time_stamp(9:10)), str2double(time_stamp(11:12)), ...
                                str2double(time_stamp(13:14)));
fname.datevec       = datevec(fname.datenum);

% =====================================================================
% Prepare log file
% =====================================================================
filedir = fileparts(param.log_fn);
if ~exist(filedir,'dir')
  mkdir(filedir);
end
[log_fid,msg] = fopen(param.log_fn,'w');
if log_fid <= 1
  fprintf('Opening of log file failed, using stdout.\n');
  fprintf('  msg\n');
  log_fid = 1;
else
  fprintf('Opened log file %s\n', param.log_fn);
end
[status,result] = system('hostname');
fprintf(log_fid, 'hostname: %s', result);


% =====================================================================
% Get the headers from all the files
% =====================================================================

for file_idx_idx = 1:length(file_idxs)
  file_idx          = file_idxs(file_idx_idx);
  fn                = filenames{file_idx};
  [path name]       = fileparts(fn);
  fprintf(log_fid, '  File %s (%.1f sec)\n', name, toc);
  
  % currently using old scripts, will write in as possible
  Header                    = load_header_acords(fn);
  [TimeComputer TimeRadar]  = load_time_acords(fn,Header);
  
  % convert to newer CSARP variable structure  
  %[tmp_hdr wfs load_param] = load_mcords_hdr(fn, load_param);
  %%load_mcords_hdr is a MEX file, need to find structure of output?
  tmp_hdr.ver         = 1;
  tmp_hdr.time        = TimeComputer;
  tmp_hdr.epri        = floor(TimeRadar./(Header.PreSum*Header.NumberWaveforms/Header.PRF)); % creates an integer step mapping identical to EPRI concept
  
  tmp_hdr.offset      = Header.IndexData + 4+4+4+8+4 + Header.IndexRecordStart(adc_idx,1) - 1 + ...
                        (0:(Header.NumberRecords-1)) .* [4+4+4+8+4+2*Header.IndexRecordStop(end,end)];                  

  tmp_hdr.offset = tmp_hdr.offset + Header.IndexRecordStart(adc_idx,1) - 1;
     
  for idx_wfs = 1:Header.NumberWaveforms
      wfs.num_sam(idx_wfs)      = Header.Waveform(idx_wfs).NumberSamples(1);
      wfs.which_bits(idx_wfs)   = 1;
      wfs.samp_delay(idx_wfs)   = Header.Waveform(idx_wfs).SampleDelay(1);
      wfs.presums(idx_wfs)      = Header.PreSum;
  end
                    
  if file_idx_idx == 1
    hdrs.ver            = tmp_hdr.ver.*ones(1,length(tmp_hdr.offset),'uint16');
    hdrs.time           = tmp_hdr.time;
    hdrs.epri           = tmp_hdr.epri;
    hdrs.file_idx       = file_idx*ones(1,length(tmp_hdr.offset),'uint16');
    hdrs.offset         = tmp_hdr.offset;
  else
    hdrs.ver            = [hdrs.ver tmp_hdr.ver.*ones(1,length(tmp_hdr.offset),'uint16')];
    hdrs.time           = [hdrs.time tmp_hdr.time];
    hdrs.epri           = [hdrs.epri tmp_hdr.epri];
    hdrs.offset         = [hdrs.offset tmp_hdr.offset];
    hdrs.file_idx       = [hdrs.file_idx file_idx*ones(1,length(tmp_hdr.offset),'uint16')];
  end
end

success=1;

return;
