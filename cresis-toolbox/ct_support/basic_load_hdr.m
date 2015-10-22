function [file_size offset varargout] = basic_load_hdr(fn, param);
% [file_size offset varargout] = basic_load_hdr(fn, param);
%
% Loads header fields from radar data files that use frame sync. Function
% used by create_segment_raw_file_list_v2.m.
%
% fn = string containing filename to load
% param = structure controlling loading
%  .field_offsets = row vector containing the offsets to each 32 bit sized header
%    field that will be returned. Size of varargout is equal to the length of
%    this vector.  For example, [1 2 3], would return the 3 32 bit values
%    preceding each frame sync.
%  .frame_sync = One of the following frame sync values:
%    hex2dec('DEADBEEF'), hex2dec('1ACFFC1D'), hex2dec('BADA55E5')
%  .byte_offsets = optional field that specifies which byte boundaries the
%    frame sync will be searched on relative to 32 bit byte boundaries. The 
%    default is [0].
%    Setting .byte_offsets to [0 1 2 3] will check all byte offsets in that
%    order.
%   .byte_offsets_threshold = scalar, default is 100
%     If more than this many matches are found with a particular byte offset,
%     the program assumes this is the correct byte offset and does not try
%     any more byte offsets.
%
% [file_size offset fields1 ... fieldsN] = basic_load_hdr(fn, param);
%
% See also: create_segment_raw_file_list_v2.m

if ~isfield(param,'byte_offsets')
  param.byte_offsets = [0];
end

if ~isfield(param,'byte_offsets_threshold')
  param.byte_offsets_threshold = 100;
end

% FAILED EXPERIMENT:
%   Loading one byte at a time is slow
%
% fid = fopen(fn,'r','ieee-be');
% 
% data = fread(fid,'uint8').';
% 
% file_size = ftell(fid);
% 
% fclose(fid);
% 
% best_byte_offset = 0;
% best_score = -1;
% for byte_offset = param.byte_offsets
%   frame_sync_1 = mod(param.frame_sync, 2^8);
%   frame_sync_2 = mod( floor(param.frame_sync/2^8), 2^8);
%   frame_sync_3 = mod( floor(param.frame_sync/2^16), 2^8);
%   frame_sync_4 = floor(param.frame_sync/2^24);
%   
%   offset2 = find(data(byte_offset+1:end-3) == frame_sync_1 ...
%     & data(byte_offset+2:end-2) == frame_sync_2 ...
%     & data(byte_offset+3:end-1) == frame_sync_3 ...
%     & data(byte_offset+4:end) == frame_sync_4);
%   if length(offset2) > best_score
%     best_score = length(offset2);
%     best_byte_offset = byte_offset;
%   end
% end

cur_byte_offset = param.byte_offsets(1);

fid = fopen(fn,'r','ieee-be');

fseek(fid, cur_byte_offset, -1);

data = fread(fid,'uint32').';

file_size = ftell(fid);

fclose(fid);

offset = find(data == param.frame_sync);

% Check other offsets if very few frame syncs were found
best_score = length(offset);
best_byte_offset = cur_byte_offset;
if length(param.byte_offsets) > 1 && best_score <= param.byte_offsets_threshold
  
  if cur_byte_offset ~= 0
    fid = fopen(fn,'r','ieee-be');
    data = fread(fid,'uint32').';
    file_size = ftell(fid);
    fclose(fid);
  end

  for byte_offset = param.byte_offsets(2:end)
    frame_sync_1 = floor(param.frame_sync/2^(8*byte_offset));
    frame_sync_2 = mod(param.frame_sync,2^(8*byte_offset));
    offset2 = find(mod(data(1:end-1),2^(8*(4-byte_offset))) == frame_sync_1 ...
      & floor(data(2:end)/2^(8*(4-byte_offset))) == frame_sync_2);
    if length(offset2) > best_score
      best_score = length(offset2);
      best_byte_offset = byte_offset;
    end
    if best_score > param.byte_offsets_threshold
      % With this many matches, we assume that we have the right byte offset
      % and don't try the others. If the threshold is too low, the wrong
      % byte offset could be chosen and lead to header reader errors.
      break;
    end
  end
  if best_byte_offset ~= cur_byte_offset
    fprintf('  best byte offset %d with %d score\n', best_byte_offset, best_score);
    
    fid = fopen(fn,'r','ieee-be');
    
    fseek(fid, best_byte_offset, -1);
    
    data = fread(fid,'uint32').';
    
    file_size = ftell(fid);
    
    fclose(fid);
    
    offset = find(data == param.frame_sync);
  end
end

offset = offset(1:end-1);

for field_idx = 1:length(param.field_offsets)
  varargout{field_idx} = data(offset + param.field_offsets(field_idx));
end

%     % Open the file and get the header fields out of it using the frame
%     % sync locations from finfo.syncs
%     fid = fopen(fn,'r','ieee-be');
%     unknown = zeros(size(offset));
%     seconds = zeros(size(offset));
%     fraction = zeros(size(offset));
%     for offset_idx = 1:length(offset)
%       fseek(fid,offset(offset_idx)+4,-1);
%       unknown(offset_idx) = fread(fid,1,'uint32');
%       seconds(offset_idx) = fread(fid,1,'uint32');
%       fraction(offset_idx) = fread(fid,1,'uint32');
%     end
%     fclose(fid);

offset = (offset-1) *4 + best_byte_offset;


return;