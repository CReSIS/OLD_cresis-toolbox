function [finfo] = frame_sync_info(filename,param)
% [finfo] = frame_sync_info(filename,param)
%
% Returns basic information about frame syncs. Used with 
% MCoRDS and 1U DAQ raw files.
%
% param = struct
%   .fast = only looks in the beginning of the file for the
%      correct starting point, some finfo fields not filled.
%      Only looks at the first 20 syncs. Default is false.
%      (scalar binary)
%   .first_byte = First byte to use. Default is zero. Primarily
%      used for skipping over "garbage" data at the front of the file
%      (scalar integer)
%   .cont_mode = Continuous load mode (default is true)
%   .sync = hexidecimal format string like 'DEADBEEF', 'BADA55E5',
%     '1ACFFC1D'
%   .open_mode = open mode ('ieee-be' is default), see Matlab fopen
%
% Examples of usage:
%   filename = 'mcords.rec002.0001.r4-1.20100510065757.dat';
%   finfo = frame_sync_info(filename);
%   finfo = frame_sync_info(filename,struct('fast',true))
%
%   filename = 'mcords.rec002.0000.r4-1.20100510065757.dat';
%   finfo = frame_sync_info(filename,struct('first_byte',2^26))
%   finfo = frame_sync_info(filename,struct('fast',0,'first_byte',2^26))
%
% Author: John Paden

% ==================================================================
% Check user arguments
% ==================================================================
if ~exist('param','var') || isempty(param)
  param = [];
end
if ~isfield(param,'fast') || isempty(param.fast)
  param.fast = 0;
end
if ~isfield(param,'first_byte') || isempty(param.first_byte)
  param.first_byte = 0;
end
if ~isfield(param,'cont_mode') || isempty(param.cont_mode)
  param.cont_mode = true;
end
if ~isfield(param,'sync') || isempty(param.sync)
  param.sync = 'DEADBEEF';
end
if ~isfield(param,'rec_size') || isempty(param.rec_size)
  param.rec_size = [];
end
if ~isfield(param,'open_mode') || isempty(param.open_mode)
  param.open_mode = 'ieee-be';
end

% ==================================================================
% Open file
% ==================================================================
finfo.filename = filename;

[fid,msg] = fopen(filename,'r',param.open_mode);
if fid < 1
  fprintf('Could not open file %s\n', filename);
  error(msg);
end
fseek(fid,param.first_byte,'bof');

% =====================================================================
% Find frame sync (aka magic) 0xDEADBEEF
% =====================================================================
if strcmpi(param.open_mode,'ieee-be')
  frame_sync(1) = hex2dec(param.sync(1:2));
  frame_sync(2) = hex2dec(param.sync(3:4));
  frame_sync(3) = hex2dec(param.sync(5:6));
  frame_sync(4) = hex2dec(param.sync(7:8));
else
  frame_sync(1) = hex2dec(param.sync(7:8));
  frame_sync(2) = hex2dec(param.sync(5:6));
  frame_sync(3) = hex2dec(param.sync(3:4));
  frame_sync(4) = hex2dec(param.sync(1:2));
end

done = false;
finfo.syncs = [];
while ~done && ~feof(fid)
  % Read in overlapping blocks in case frame sync lies on a
  % block boundary
  fseek(fid,-length(frame_sync),'cof');
  block_ind = ftell(fid);
  % Read in a 2^16 byte block at a time and search of syncs
  % Byte alignment is unknown so we have to read in bytes
  % so we can sweep 4 byte frame_sync across all byte boundaries
  data_block = fread(fid,2^16,'uint8');
  % Look for matches on the first frame sync
  inds = find(data_block(1:end-(length(frame_sync-1)))==frame_sync(1));
  for idx = inds.'
    % Match the rest of the frame sync
    if length(frame_sync) == 1 || all(data_block(idx+(1:length(frame_sync)-1)).' == frame_sync(2:end))
      % Found a frame sync (aka magic)
      finfo.syncs(end+1) = block_ind + idx - 1;
      if param.fast && length(finfo.syncs) > 100
        done = true;
      end
    end
  end
end

fseek(fid,0,'eof');
finfo.file_size = ftell(fid);
fclose(fid);

if length(finfo.syncs) <= 1
  error('File too short, <= 1 syncs found');
end

if isempty(param.rec_size)
  % Estimate record size from file
  finfo.rec_size = median(diff(finfo.syncs));
  finfo.first_good = find(diff(finfo.syncs) == finfo.rec_size,1);
  if isempty(finfo.first_good)
    warning('Record size is ambiguous and median returned a number half way between two possible record sizes');
    finfo.first_good = 1;
  end
else
  finfo.rec_size = param.rec_size;
  finfo.first_good = find(diff(finfo.syncs) == finfo.rec_size,1);
  if isempty(finfo.first_good)
    warning('Record size is ambiguous and specified value returned no records');
    finfo.first_good = 1;
  end
end
  
if param.cont_mode
  finfo.num_rec = ceil((finfo.file_size-finfo.syncs(finfo.first_good))/finfo.rec_size);
else
  finfo.num_rec = length(finfo.syncs);
  while finfo.file_size - finfo.syncs(finfo.num_rec) < finfo.rec_size
    % Last record must be a full record, we are removing this last sync
    finfo.num_rec = finfo.num_rec - 1;
  end
end

finfo.ideal_sync_pos = finfo.syncs(finfo.first_good) ...
  + finfo.rec_size * (0:finfo.num_rec-1);

if param.fast
  return;
end

finfo.bad_mask = zeros(1,finfo.num_rec);
for ideal_rec = 1:length(finfo.ideal_sync_pos)
  rec_match = find(finfo.ideal_sync_pos(ideal_rec) == finfo.syncs,1);
  if isempty(rec_match)
    % No sync mark at this location
    finfo.bad_mask(ideal_rec) = 1;
  else
    % Check if this is a short record
    if rec_match == length(finfo.syncs)
      % Last record special case
      if finfo.file_size - finfo.syncs(rec_match) < finfo.rec_size
        finfo.bad_mask(ideal_rec) = 1;
      end
    else
      % Not last record
      if mod(finfo.syncs(rec_match+1) - finfo.syncs(rec_match),finfo.rec_size)
        finfo.bad_mask(ideal_rec) = 1;
      end
    end
  end
  
end

return;

