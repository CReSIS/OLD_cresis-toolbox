function [hdr wfs status] = load_mcords_hdr_mfile(fn, param)
% [hdr wfs status] = load_mcords_hdr_mfile(fn, param)
%
% Used by create_records_mcords_task to load in the headers.
% It is robust to errors in the digital system and is efficient.
%
%  INPUTS:
%    fn = filename string
%    param = parameter structure
%      first_byte = first byte to read in (usually 0, except for the
%        first file index, then see example below)
%      fs = sampling frequency of radar, used to interpret counts/tics
%        in header field variables
%      rec_size = record size (see example for how to compute)
%      sync_lock = is the sync locked (i.e. does first_byte point to
%        a valid record)
%      last_bytes = bytes to prepend to the data, first_byte is relative
%        to these
%
%  OUTPUTS:
%    hdr.ver = num_recx1 vector of uint8
%    hdr.time = num_recx1 vector of double
%    hdr.epri = num_recx1 vector of uint32
%    hdr.offset = num_recx1 vector of int32
%      The first element may be negative which means that the start of
%      this record is in the previous file.
%
%    wfs.num_sam = hdr.num_wf x 1 vector of double
%    wfs.which_bits = hdr.num_wf x 1 vector of double
%    wfs.t0 = hdr.num_wf x 1 vector of double
%    wfs.presums = hdr.num_wf x 1 vector of double
%
%    status.last_bytes = last bytes in the file
%    status.sync_lock = current status of the sync

% HEADER FILE FORMAT
% Header Size: 4*8 + 16*8 = 160 bytes
% struct header
% {
%   uint32 frame_sync; // aka magic
%   uint32 version;    // file type version (e.g. May 2010 is 2)
%   uint32 seconds;    // # of seconds
%   uint32 fraction;   // # of fs/2 clock cycles
%   uint32 epri;       // effective PRI number
%   uint32 num_wf;     // number of waveforms
%   uint32 which_bits; // ???
%   uint32 dec_cfg;    // decimation config settings
%   struct wf wfs[16];
% }
% struct wf
% {
%   uint32(31:14) tmp;
%   uint32(13:0) num_sam;    // number of samples
%   uint32(31:28) tmp;
%   uint32(27:24) which_bits // whichbits for waveform
%   uint32(23:10) t0 // sample delay in fs clock cycles + 11 us
%   uint32(9:0) presums      // number of presums - 1
% }

rec_size = param.rec_size;
fid = fopen(fn, 'r', 'ieee-be');

% =======================================================================
% Determine the maximum number of complete records in the file
fseek(fid, 0, 'eof');
file_size = ftell(fid);
num_rec = floor((length(param.last_bytes) + file_size - param.first_byte)/rec_size);
% fprintf('Reading %d records\n', num_rec); % DEBUG

% Preallocate arrays
hdr.ver = zeros(num_rec,1,'uint8');
hdr.time = zeros(num_rec,1,'double');
hdr.epri = zeros(num_rec,1,'uint32');
hdr.offset = zeros(num_rec,1,'int32');

% =======================================================================
% Go to first record
fseek(fid, param.first_byte, 'bof');

% =======================================================================
% If sync is lost, find it before reading in records
status.sync_lock = param.sync_lock;
if ~status.sync_lock
  status.sync_lock = find_next_sync(fid, file_size - param.first_byte - rec_size);
  if (status.sync_lock)
    fprintf('No good records in this file.');
    return;
  end
end

% =======================================================================
% =======================================================================
% Read in headers
% =======================================================================
% =======================================================================

% =======================================================================
% Read in NUM_REC_BLOCK records at a time to speed up loading
NUM_REC_BLOCK = 2048;

rec = 0;
frame_sync = [hex2dec('DE'); hex2dec('AD'); hex2dec('BE'); hex2dec('EF')];
while ftell(fid) <= file_size-rec_size && status.sync_lock
  
  % Determine the number of records left in the file
  num_rec_left = floor((file_size-ftell(fid))/rec_size);
  if (num_rec_left > NUM_REC_BLOCK)
    num_rec_left = NUM_REC_BLOCK;
  end
  
  % Handle the beginning case when the user may have passed in some
  % bytes from the previous file and we need to prepend those to our
  % result.
  if (rec == 0 && length(param.last_bytes) > 0)
    [data_buffer bytes_read] = fread(fid,rec_size*num_rec_left - length(param.last_bytes),'uint8');
    data_buffer = [param.last_bytes; data_buffer];
  else
    [data_buffer bytes_read] = fread(fid,rec_size*num_rec_left,'uint8');
    bytes_read;
  end
  
  for sub_rec = 0:num_rec_left-1
    % ==============================================================
    % Check for bad sync marker
    if data_buffer(sub_rec*rec_size+0*4 + (1:4)) ~= frame_sync
      fprintf('  ==> Lost sync lock (0xDEADBEEF frame sync not found, %d)\n', ...
        ftell(fid)-bytes_read+sub_rec*rec_size);
      fseek(fid,-bytes_read+sub_rec*rec_size, 0);
      status.sync_lock = find_next_sync(fid, file_size - param.first_byte - rec_size);
      if ~status.sync_lock
        fprintf('  No more locks in this file.');
        return;
      end
      break;
    end
    
    % ==============================================================
    % Read, parse, and interpret header fields
    %  uint32 frame_sync; // aka magic
    %  uint32 version;    // file type version (e.g. May 2010 is 2)
    %  uint32 seconds;    // # of seconds
    %  uint32 fraction;   // # of fs/2 clock cycles
    %  uint32 epri;       // effective PRI number
    %  uint32 num_wf;     // number of waveforms
    %  uint32 which_bits; // ???
    %  uint32 dec_cfg;    // decimation config settings
    % ==============================================================
    % fprintf('0x%s\n", dec2hex(big2littleint32(&data_buffer[sub_rec*rec_size+0*4])); // DEBUG
    
    rec = rec + 1;
    hdr.ver(rec) = uint8(data_buffer(sub_rec*rec_size+1*4 + 4));
    seconds = sum(data_buffer(sub_rec*rec_size+2*4 + (1:4)).*[2^24; 2^16; 2^8; 1]);
    fractions = sum(data_buffer(sub_rec*rec_size+3*4 + (1:4)).*[2^24; 2^16; 2^8; 1]);
    hdr.time(rec) = seconds + 2*(fractions)/param.fs;
    hdr.epri(rec) = uint32(sum(data_buffer(sub_rec*rec_size+4*4 + (1:4)).*[2^24; 2^16; 2^8; 1]));
    hdr.offset(rec) = int32(ftell(fid)-(bytes_read+length(param.last_bytes))+sub_rec*rec_size);
    
    % ==============================================================
    % Read in waveform header
    %  uint32(31:14) tmp;
    %  uint32(13:0) num_sam;    // number of samples
    %  uint32(31:28) tmp;
    %  uint32(27:24) which_bits // whichbits for waveform
    %  uint32(23:10) t0 // sample delay in fs clock cycles
    %  uint32(9:0) presums      // number of presums - 1
    if rec == 1
      num_wfs = data_buffer(sub_rec*rec_size+5*4 + 4);
      wf = 0;
      while wf < num_wfs
        dword = sum(data_buffer(sub_rec*rec_size+8*4 + 8*wf + (1:4)).*[2^24; 2^16; 2^8; 1]);
        wfs.num_sam(wf+1) = mod(dword, 2^14);
        dword = sum(data_buffer(sub_rec*rec_size+9*4 + 8*wf + (1:4)).*[2^24; 2^16; 2^8; 1]);
        wfs.which_bits(wf+1) = mod(floor(dword/2^24), 2^4);
        wfs.t0(wf+1) = mod(floor(dword/2^10), 2^14) / param.fs - 11e-6;
        wfs.presums(wf+1) = mod(dword, 2^10) + 1;
        wf = wf + 1;
        % fprintf('%d of %d %f\n', wf, num_wfs, wfs(wf).presums); % DEBUG
      end
    end
  end
  param.last_bytes = [];
end


% fprintf('Read %d records of %d\n', rec, num_rec); % DEBUG
if (rec ~= num_rec)
  hdr.ver = hdr.ver(1:rec);
  hdr.time = hdr.time(1:rec);
  hdr.epri = hdr.epri(1:rec);
  hdr.offset = hdr.offset(1:rec);
end

% =======================================================================
% Go to last incomplete record and copy it
last_bytes_out_len = (file_size - ftell(fid));

status.last_bytes = fread(fid, last_bytes_out_len, 'uint8');

fclose(fid);

return;

% =====================================================================
% Find next frame sync (aka magic) 0xDEADBEEF
% =====================================================================
function sync_lock = find_next_sync(fid, bytes_left)

DATA_BLOCK_SIZE = 2^16;

frame_sync(1) = hex2dec('DE');
frame_sync(2) = hex2dec('AD');
frame_sync(3) = hex2dec('BE');
frame_sync(4) = hex2dec('EF');

total_read = 0;
while (total_read < bytes_left)
  
  if (bytes_left-total_read > DATA_BLOCK_SIZE)
    num_to_read = DATA_BLOCK_SIZE;
  else
    num_to_read = bytes_left-total_read;
  end
  
  % Read in overlapping blocks in case frame sync lies on a
  % block boundary
  fseek(fid,-length(frame_sync),0);
  block_ind = ftell(fid);
  
  % Read in a 2^16 byte block at a time and search of syncs
  % Byte alignment is unknown so we have to read in bytes
  % so we can sweep 4 byte frame_sync across all byte boundaries
  [data_block num_read] = fread(fid,num_to_read,'uint8');
  total_read = total_read + num_read;
  
  % Look for matches on the first frame sync
  inds = find(data_block(1:end-(length(frame_sync-1)))==frame_sync(1));
  for idx = inds.'
    % Match the rest of the frame sync
    if length(frame_sync) == 1 || all(data_block(idx+(1:length(frame_sync)-1)).' == frame_sync(2:end))
      % Found a frame sync (aka magic)
      fseek(fid, -num_read + idx - 1, 0);
      fprintf(' ==> Relocked (%d)\n', ftell(fid));
      sync_lock = 1;
      return;
    end
  end
  
end

return;


% =========================================================================
% =========================================================================
% Examples
% =========================================================================
% =========================================================================

fn = get_filenames('/cresis/data2/MCoRDS/2010_Antarctica/20100103A/chan1/','r3-1','','0085.dat');
fn = get_filenames('/cresis/data2/MCoRDS/2009_Chile/20091031/chan1/flight1011_LA_1_and_16presum_daqch1/flight1011_daqch1_rec1/','r1-1','','0010.dat');
fn = get_filenames('/cresis/data2/MCoRDS/2010_Greenland_DC8/20100324A/chan1/','','0010.r4','');
fn = get_filenames('/cresis/data2/MCoRDS/2009_Chile/20091020/chan1/seg5/','','','0012.dat');
fn = get_filenames('/cresis/data4/MCoRDS/2011_Greenland_TO/20110502/chan1/seg1/','','','0012.dat');
fn = fn{1};
syncs = get_first10_sync(fn, 0);
param.rec_size = median(diff(syncs));
first_byte = find(diff(syncs) == param.rec_size, 1);
param.first_byte = syncs(first_byte);
param.sync_lock = 1;
param.last_bytes = '';
param.fs = 1e9/8;
tic; [hdr wfs status] = load_mcords_hdr_mfile(fn, param); toc;
tic; [hdr2 wfs2 status2] = load_mcords_hdr(fn, param); toc;


