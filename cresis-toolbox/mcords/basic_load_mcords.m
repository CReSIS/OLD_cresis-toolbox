function [hdr,data] = basic_load_mcords(fn, param)
% [hdr,data] = basic_load_mcords(fn, param)
%
% This is the only function which loads raw data directly from MCoRDS.
%
% Loads a single mcords radar file. This is primarily for debugging.
% NOTE: 64-bit computer is required to load a 512 MB file since it will
% consume 1 GB of memory after loading.
%
% fn = filename of MCoRDS data
% param = struct controlling loading of data
%   .clk = clock (Hz), default one, used to interpret
%     counts in the header fields
%     MCoRDS uses sampling frequency 1e9/9
%   .wfs = vector of waveforms to load (one indexed)
%   .recs = 2 element vector for records to load [start_rec num_rec]
%     start_rec uses zero-indexing
%   .first_byte = First byte to use. Default is zero. Primarily
%      used for skipping over "garbage" data at the front of the file
%      (scalar integer)
%
% hdr = file header for each record
% data = cell vector of single matrices of radar data where each entry
%   in the cell vector is a waveform
%
% Examples: run_basic_load_mcords
%
% Authors: John Paden

SAMP_SIZE = 2;

% ===================================================================
% Check input arguments
% ===================================================================
if ~exist('param','var') || isempty(param)
  param.clk = 1e9/9;
  param.wfs = [];
  param.recs = [];
  param.first_byte = 0;
end
if ~isfield(param,'clk');
  param.clk = 1e9/9;
end
if ~isfield(param,'wfs');
  param.wfs = [];
end
if ~isfield(param,'recs');
  param.recs = [];
end
if ~isfield(param,'first_byte');
  param.first_byte = 0;
end

% ===================================================================
% Get the locations of the first sync markers
% ===================================================================
syncs = get_first10_sync_mfile(fn, param.first_byte);
% Find the record size in a robust way
hdr.finfo.rec_size = median(diff(syncs));
% Find the first good header
hdr.finfo.syncs(1) = syncs(find(diff(syncs) == hdr.finfo.rec_size, 1));

% ===================================================================
% Determine the size of the output matrix
% ===================================================================
[fid,msg] = fopen(fn,'r','ieee-be');
if fid < 1
  fprintf('Could not open file %s\n', fn);
  error(msg);
end

fseek(fid,0,'eof');
hdr.finfo.file_size = ftell(fid);

hdr.finfo.num_rec = floor((hdr.finfo.file_size-hdr.finfo.syncs(1))/hdr.finfo.rec_size);

% Number of slow time records
if isempty(param.recs)
  recs = [0 hdr.finfo.num_rec];
else
  if param.recs(1) + param.recs(2) > hdr.finfo.num_rec
    error('Only %d records\n', hdr.finfo.num_rec);
  end
  recs = param.recs;
end

% Seek to first record
fseek(fid, hdr.finfo.syncs(1) + hdr.finfo.rec_size*recs(1), 'bof');

% Load in header data all at once
hdr_data = fread(fid, [40 recs(2)],'40*uint32',hdr.finfo.rec_size-160);

bad_syncs = find(hdr_data(1,:) ~= hex2dec('DEADBEEF'));
if ~isempty(bad_syncs)
  fprintf('  Bad sync marks, rereading data file the slow way\n');
  fprintf('  Searching for individual sync marks\n');
  hdr.finfo = frame_sync_info(fn,struct('cont_mode',false));
  fprintf('  Reading in data file\n');

  % Number of slow time records (recalculate with new information)
  if isempty(param.recs)
    recs = [0 hdr.finfo.num_rec];
  else
    if param.recs(1) + param.recs(2) > hdr.finfo.num_rec
      error('Only %d records\n', hdr.finfo.num_rec);
    end
    recs = param.recs;
  end

  % --------------------------------------------------------------
  % Read in first header so we can preallocate arrays
  fseek(fid,hdr.finfo.syncs(1),-1);
  hdr.frame_sync = fread(fid,1,'uint32');
  hdr.version = fread(fid,1,'uint32');
  hdr.seconds = fread(fid,1,'uint32');
  hdr.fraction = fread(fid,1,'uint32');
  hdr.epri = fread(fid,1,'uint32');
  hdr.num_wf = fread(fid,1,'uint32');
  hdr.bit_shifts = fread(fid,1,'uint32');
  hdr.dec_cfg = fread(fid,1,'uint32');
  wfs_raw = fread(fid,[2 16],'uint32');
  for wf = 1:hdr.num_wf(1)
    hdr.wfs(wf).num_sam = mod(wfs_raw(1,wf), 2^14);
    hdr.wfs(wf).bit_shifts = mod(floor(wfs_raw(2,wf) / 2^24), 2^4);
    hdr.wfs(wf).t0 = mod(floor(wfs_raw(2,wf) / 2^10), 2^14) / param.clk - 10.8e-6;
    hdr.wfs(wf).presums = 1 + mod(wfs_raw(2,wf), 2^10);
  end
  % -------------------------------------------------------------
  % Determine which waveforms are being loaded
  if isempty(param.wfs)
    if hdr.num_wf(1) == 0
      error('Zero waveforms header error: consider not using this file.');
    end
    param.wfs = 1:hdr.num_wf(1);
  end
  % -------------------------------------------------------------
  % Determine waveform offsets into each record
  offsets(1) = 160;
  for wf = 2:hdr.num_wf(1)
    offsets(wf) = offsets(wf-1) + hdr.wfs(wf-1).num_sam*SAMP_SIZE;
  end
  % --------------------------------------------------------------
  % Preallocate arrays
  hdr.frame_sync = zeros(1,recs(2));
  hdr.version = zeros(1,recs(2));
  hdr.seconds = zeros(1,recs(2));
  hdr.fraction = zeros(1,recs(2));
  hdr.epri = zeros(1,recs(2));
  hdr.num_wf = zeros(1,recs(2));
  hdr.bit_shifts = zeros(1,recs(2));
  hdr.dec_cfg = zeros(1,recs(2));
  for wf_idx = 1:length(param.wfs)
    wf = param.wfs(wf_idx);
    % Create the outputs matrix
    data{wf_idx} = zeros(hdr.wfs(wf).num_sam, recs(2), 'single');
  end
  % --------------------------------------------------------------
  % Read in data
  for rec = 1:recs(2)
    fseek(fid,hdr.finfo.syncs(rec),-1);
    hdr.frame_sync(rec) = fread(fid,1,'uint32');
    hdr.version(rec) = fread(fid,1,'uint32');
    hdr.seconds(rec) = fread(fid,1,'uint32');
    hdr.fraction(rec) = fread(fid,1,'uint32');
    hdr.epri(rec) = fread(fid,1,'uint32');
    hdr.num_wf(rec) = fread(fid,1,'uint32');
    hdr.bit_shifts(rec) = fread(fid,1,'uint32');
    hdr.dec_cfg(rec) = fread(fid,1,'uint32');
    % -------------------------------------------------------------
    % Load data matrices
    for wf_idx = 1:length(param.wfs)
      wf = param.wfs(wf_idx);
      % Seek to the waveform
      fseek(fid, hdr.finfo.syncs(rec) + offsets(wf), -1);
      data{wf_idx}(:,rec) = fread(fid, hdr.wfs(wf).num_sam, 'uint16');
    end
  end
  fprintf('  Done reading data file\n');
  
else
  % -------------------------------------------------------------
  % Parse header data matrix into header fields
  hdr.frame_sync = hdr_data(1,:);
  hdr.version = hdr_data(2,:);
  hdr.seconds = hdr_data(3,:);
  hdr.fraction = hdr_data(4,:);
  hdr.epri = hdr_data(5,:);
  hdr.num_wf = hdr_data(6,:);
  hdr.bit_shifts = hdr_data(7,:);
  hdr.dec_cfg = hdr_data(8,:);
  wfs_raw = reshape(hdr_data(9:40,1), [2 16]);
  for wf = 1:hdr.num_wf(1)
    hdr.wfs(wf).num_sam = mod(wfs_raw(1,wf), 2^14);
    hdr.wfs(wf).bit_shifts = mod(floor(wfs_raw(2,wf) / 2^24), 2^4);
    hdr.wfs(wf).t0 = mod(floor(wfs_raw(2,wf) / 2^10), 2^14) / param.clk - 10.8e-6;
    hdr.wfs(wf).presums = 1 + mod(wfs_raw(2,wf), 2^10);
  end
  % -------------------------------------------------------------
  % Determine waveform offsets into each record
  offsets(1) = 160;
  for wf = 2:hdr.num_wf(1)
    offsets(wf) = offsets(wf-1) + hdr.wfs(wf-1).num_sam*SAMP_SIZE;
  end
  % -------------------------------------------------------------
  % Determine which waveforms are being loaded
  if isempty(param.wfs)
    param.wfs = 1:hdr.num_wf(1);
  end
  % -------------------------------------------------------------
  % Load data matrices
  for wf_idx = 1:length(param.wfs)
    wf = param.wfs(wf_idx);
    % Seek to waveform in first record
    fseek(fid, hdr.finfo.syncs(1) + hdr.finfo.rec_size*recs(1) + offsets(wf), 'bof');
    
    prec_str = sprintf('%d*uint16=>single', hdr.wfs(wf).num_sam);
    data{wf_idx} = fread(fid, [hdr.wfs(wf).num_sam recs(2)], prec_str, ...
      hdr.finfo.rec_size-hdr.wfs(wf).num_sam*SAMP_SIZE);
  end
end

fclose(fid);

hdr.time_sod = double(hdr.seconds) + double(hdr.fraction)/(param.clk/2);
hdr.utc_time_sod = double(hdr.seconds) + double(hdr.fraction)/(param.clk/2);

return;
