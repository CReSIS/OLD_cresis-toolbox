function sync_offsets = get_first10_sync_mfile(data_fn, first_byte, param)
% sync_offsets = get_first10_sync_mfile(data_fn, first_byte, param)
%
% Gets the location of the first 10 sync markers (0xDEADBEEF)
%
% data_fn = string, data file to open
% first_byte = first byte to read in (typically zero)
% param = optional argument
%  .sync = hexidecimal format string like 'DEADBEEF' or 'BADA55E5'
%  .num_sync = default is 10 (leave empty/undefined for default), number
%    of syncs returned will be at least this many (could be more)
%  .last = logical that when set to true reads syncs starting at the
%    end of the file (default is false and reads syncs starting at the
%    beginning of the file)
%
% Authors: John Paden
%
% See also: frame_sync_info.m, get_first10_sync.cpp,
%   create_records_mcords_task.m

% ==================================================================
% Check user arguments
% ==================================================================
if ~exist('param','var') || isempty(param)
  param = struct();
end
if ~isfield(param,'sync') || isempty(param.sync)
  param.sync = 'DEADBEEF';
end
if ~isfield(param,'num_sync') || isempty(param.num_sync)
  param.num_sync = 10;
end
if ~isfield(param,'last') || isempty(param.last)
  param.last = false;
end

% ==================================================================
% Open file
% ==================================================================
[fid,msg] = fopen(data_fn, 'r');
if fid < 1
  fprintf('Could not open file %s\n', data_fn);
  error(msg);
end
if param.last
  fseek(fid,first_byte,'eof');
else
  fseek(fid,first_byte,'bof');
end

% =====================================================================
% Find frame sync (aka magic) 0xDEADBEEF
% =====================================================================
frame_sync(1) = hex2dec(param.sync(1:2));
frame_sync(2) = hex2dec(param.sync(3:4));
frame_sync(3) = hex2dec(param.sync(5:6));
frame_sync(4) = hex2dec(param.sync(7:8));
if param.last
  frame_sync = frame_sync(end:-1:1);
end

done = false;
sync_offsets = [];
while ~done && ~feof(fid)
  if ~param.last
    % Read in overlapping blocks in case frame sync lies on a
    % block boundary
    fseek(fid,-length(frame_sync),'cof');
    block_ind = ftell(fid);
    % Read in a 2^16 byte block at a time and search of syncs
    % Byte alignment is unknown so we have to read in bytes
    % so we can sweep 4 byte frame_sync across all byte boundaries
    data_block = fread(fid,2^16,'uint8');
  else
    fseek(fid,length(frame_sync) - 2^16,'cof');
    block_ind = ftell(fid);
    data_block = fread(fid,2^16,'uint8');
    fseek(fid,-2^16,'cof');
    data_block = data_block(end:-1:1);
  end
  % Look for matches on the first frame sync
  inds = find(data_block(1:end-(length(frame_sync)-1))==frame_sync(1));
  if isempty(inds)
    inds = [];
  end
  for idx = inds.'
    % Match the rest of the frame sync
    if length(frame_sync) == 1 || all(data_block(idx+(1:length(frame_sync)-1)).' == frame_sync(2:end))
      % Found a frame sync (aka magic), record its offset
      if ~param.last
        sync_offsets(end+1) = block_ind + (idx - 1);
      else
        sync_offsets(end+1) = block_ind + length(data_block) - idx - length(frame_sync) + 1;
      end
      if length(sync_offsets) >= param.num_sync
        done = true;
      end
    end
  end
end

fclose(fid);

return;
