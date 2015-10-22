function [success hdrs wfs] = create_records_mcords_task(filenames,param,file_idxs)
% [success hdrs wfs] = create_records_mcords_task(filenames,param,file_idxs)
%
% This function is called by create_records_mcords.m.
% It has been created to allow the cluster to be used.
%
% success = boolean which is true when the function executes properly
%   if a task fails before it can return success, then success will be
%   empty
% hdrs = struct with header information from each file
%  FIELDS?
% wfs = struct with waveform information from the first header
%  FIELDS?
%
% See also: create_records_mcords.m, load_mcords_hdr.cpp, load_mcords_hdr_mfile.m
%
% Author: John Paden

% When called from scheduler, we need to run 
try; toc; catch; tic; end

% =====================================================================
% Find the first 10 sync markers
% =====================================================================
fn = filenames{file_idxs(1)};
fname = fname_info_mcords(fn);
if fname.file_idx == 0
  % Skip bad data at the beginning of the first file
  first_byte = 2^26;
else
  first_byte = 0;
end
% Mex Function:
%syncs = get_first10_sync(fn,first_byte);
% M-file Function:
syncs = get_first10_sync_mfile(fn,first_byte);

% Find the record size in a robust way
rec_size = median(diff(syncs));
% Find the first good header
first_byte = syncs(find(diff(syncs) == rec_size, 1));

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
load_param.first_byte = first_byte;
load_param.rec_size = rec_size;
load_param.sync_lock = 1;
load_param.last_bytes = uint8([]);
load_param.fs = param.radar.fs;
for file_idx_idx = 1:length(file_idxs)
  file_idx = file_idxs(file_idx_idx);
  fn = filenames{file_idx};
  [path name] = fileparts(fn);
  fprintf(log_fid, '  File %s (%.1f sec)\n', name, toc);
  
  % Mex Function:
  [tmp_hdr wfs load_param] = load_mcords_hdr(fn, load_param);
  % M-file Function:
  %[tmp_hdr wfs load_param] = load_mcords_hdr_mfile(fn, load_param);

  % load_mcords_hdr overwrites load_param, so reset param variables
  load_param.first_byte = 0;
  load_param.rec_size = rec_size;
  load_param.fs = param.radar.fs;
  
  if file_idx_idx == 1
    hdrs.ver = tmp_hdr.ver.';
    hdrs.seconds = tmp_hdr.seconds.';
    hdrs.fractions = tmp_hdr.fractions.';
    hdrs.epri = tmp_hdr.epri.';
    hdrs.file_idx = file_idx*ones(1,length(tmp_hdr.offset),'uint16');
    hdrs.offset = tmp_hdr.offset.';
  else
    hdrs.ver = [hdrs.ver tmp_hdr.ver.'];
    hdrs.seconds = [hdrs.seconds tmp_hdr.seconds.'];
    hdrs.fractions = [hdrs.fractions tmp_hdr.fractions.'];
    hdrs.epri = [hdrs.epri tmp_hdr.epri.'];
    hdrs.offset = [hdrs.offset tmp_hdr.offset.'];
    hdrs.file_idx = [hdrs.file_idx file_idx*ones(1,length(tmp_hdr.offset),'uint16')];
  end
end
if log_fid > 1
  fclose(log_fid)
end

success = 1;

return;
