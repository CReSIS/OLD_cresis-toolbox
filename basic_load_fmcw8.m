function [hdr,data] = basic_load_fmcw8(fn,param)
% [hdr,data] = basic_load_fmcw8(fn, param)
%
% This is the only function which loads raw data directly.
%
% Loads a single fmcw8 radar file. This is primarily for debugging.
% NOTE: 64-bit computer may be essential to load a 256 MB file since it will
% consume 512 MB of memory after loading.
%
% If data is not specified as an output argument, only the header is returned
%
% fn = filename of FMCW8 data
% param = struct controlling loading of data
%   .clk = clock (Hz), default one, used to interpret
%     counts in the header fields
%     first system used 250e6
%   .recs = 2 element vector for records to load [start_rec num_rec]
%     start_rec uses zero-indexing (negative start_recs read from the
%     end of the file and only work with single header loading)
%   .file_version = default is 3
%     1: no loopback_mode or nyquist_zone setting + bug where 2 waveforms
%        were used because the DAQ could not support 1 waveform
%     2: loopback_mode or nyquist_zone settings added in place of utc_time2
%        + bug where 2 waveforms were used because the DAQ could not
%        support 1 waveform
%     3: no loopback_mode or nyquist_zone setting, 1 waveform
%     4: limited use
%     5: no loopback_mode or nyquist_zone setting, 1 waveform, 2013 Ant Basler and later used?
%   .records
%     .en = Special field for create_records_fmcw_accum.m. If true, all
%       headers will be loaded and the output arguments become:
%       hdr --> success flag (returns 1)
%       data --> hdr (all headers)
%   .debug_level = 1 is default, 2 generates plots/print-outs
%
% hdr = file header for each record (if data is not defined, behavior
%   depends on param.records.en variable, if it is false only the first hdr
%   is returned, otherwise all the headers are returned)
% data = Depends on param.records.en. When false, it is an optional output
%   array of radar data where dimensions are
%   1: fast-time/range-bin
%   2: slow-time/range-line/records
%
% Examples: See bottom of file
%   fn = '/process-archive/20170226/fmcw/snow/snow4_00_20170226_182432_0103.bin';
%   [hdr,data] = basic_load_fmcw8(fn,struct('clk',250e6));
%
% Authors: John Paden
%
% See also basic_load_*.m

% ===================================================================
% Check input arguments
% ===================================================================
if ~exist('param','var') || isempty(param)
  param.clk = 1;
  param.recs = [];
  param.debug_level = 1;
end
if ~isfield(param,'clk');
  param.clk = 1;
end
if ~isfield(param,'recs');
  param.recs = [];
end
if ~isfield(param,'file_version');
  param.file_version = 8;
end
if ~isfield(param,'records');
  param.records.en = false;
end

% Reset/clear hdr struct
hdr = [];

% ===================================================================
%% Data Format
% ===================================================================
% 32-bit header
% 32-bit EPRI
% 32-bit seconds in DCB (0 0 H H M M S S)
% 32-bit fraction
% 64-bit counter
% 64-bit computer time
% 8-bit zero/reserved
% 8-bit Nyquist zone
% 8-bit presums
% 8-bit bit shift
% 16-bit start
% 16-bit stop
% 64-bit Keysight waveform generator waveform ID
% DATA: depending on settings
% If raw data then 16 bit real with size stop-start
%
HEADER_SIZE = 48;
SAMPLE_SIZE = 2;

% ===============================================================
% Get first record position
% ===============================================================
if ~isempty(param.recs) && param.recs(1) < 0
  hdr.finfo.syncs = get_first10_sync_mfile(fn,0,struct('sync','BADA55E5','last',true));
else
  hdr.finfo.syncs = get_first10_sync_mfile(fn,0,struct('sync','BADA55E5'));
end

% ===============================================================
% Open file big-endian for reading
% ===============================================================
[fid,msg] = fopen(fn,'r','ieee-be');
if fid < 1
  fprintf('Could not open file %s\n', fn);
  error(msg);
end

fseek(fid, 0, 1);
eof_pos = ftell(fid);

% ===============================================================
% Read in waveform information + record size
% ===============================================================

if nargout < 2 && ~param.records.en
  % Seek to first record
  fseek(fid, hdr.finfo.syncs(1), -1);
  hdr.frame_sync = fread(fid,1,'uint32');
  hdr.epri = fread(fid,1,'uint32');
  hdr.seconds = fread(fid,1,'uint32').'; % From NMEA string converted to DCB
  hdr.seconds = BCD_to_seconds(hdr.seconds);
  hdr.fraction = fread(fid,1,'uint32');
  hdr.utc_time_sod = hdr.seconds + hdr.fraction / param.clk*2;
  hdr.counter = fread(fid,1,'uint64');
  hdr.comp_time_sod = double(fread(fid,1,'uint64'));
  fseek(fid,1,0);
  hdr.nyquist_zone = fread(fid,1,'uint8');
  hdr.presums = fread(fid, 1, 'uint8')+1; % presums are 0-indexed (+1)
  hdr.bit_shifts = -fread(fid, 1, 'int8');
  hdr.start_idx = fread(fid, 1, 'uint16');
  hdr.Tadc = hdr.start_idx / param.clk*2 - 10.8e-6;
  hdr.stop_idx = fread(fid, 1, 'uint16');
  hdr.waveform_ID = fread(fid,1,'uint64');
  
  fclose(fid);
  return;
elseif param.records.en
  error('Not supported');
  % Seek to first record
  fseek(fid, hdr.finfo.syncs(1) + 36, -1);
  hdr.start_idx = fread(fid, 1, 'uint16');
  hdr.stop_idx = fread(fid, 1, 'uint16');
  
  rline = 1;
  % Raw data
  hdr.num_sam(rline) = 2*(hdr.stop_idx(rline) - hdr.start_idx(rline));
  
  fseek(fid, hdr.finfo.syncs(1), -1);
  rec_size = HEADER_SIZE+SAMPLE_SIZE*hdr.num_sam;
  num_rec = floor( (eof_pos-ftell(fid)) / rec_size );
  hdr_data = fread(fid, [12 num_rec], '12*uint32', rec_size-12*4);
  data.epri = hdr_data(2,:);
  data.seconds = BCD_to_seconds(hdr_data(3,:));
  data.fraction = hdr_data(4,:);
  data.offset = hdr.finfo.syncs(1) + rec_size*(0:num_rec-1);
  fseek(fid, hdr.finfo.syncs(1) + 34, -1);
  data.wfs.presums = fread(fid, 1, 'uint8')+1;
  hdr = 1;
  fclose(fid);
  return;
end

% Seek to first record
fseek(fid, hdr.finfo.syncs(1), -1);

rline = 0;
hdr.finfo.rec_size = [];
FRAME_SYNC = hex2dec('BADA55E5');
data = zeros(0,0,'single');
while ftell(fid) <= eof_pos-HEADER_SIZE
  rline = rline + 1;
  hdr.frame_sync(rline) = fread(fid,1,'uint32');
  if hdr.frame_sync(rline) == FRAME_SYNC
    hdr.finfo.syncs(rline) = ftell(fid)-4;
  else
    % Search for next frame sync
%     keyboard
    found = false;
    while ~feof(fid)
      test = fread(fid,1,'uint32');
      if test == FRAME_SYNC
        found = true;
        break;
      end
    end
    if ~found
      rline = rline - 1;
      break;
    end
    hdr.finfo.syncs(rline) = ftell(fid)-4;
  end
  if ftell(fid) > eof_pos-HEADER_SIZE
    rline = rline - 1;
    break;
  end
  hdr.epri(rline) = fread(fid,1,'uint32');
  hdr.seconds(rline) = fread(fid,1,'uint32').'; % From NMEA string converted to DCB
  hdr.fraction(rline) = fread(fid,1,'uint32');
  hdr.counter(rline) = fread(fid,1,'uint64');
  hdr.comp_time_sod(rline) = double(fread(fid,1,'uint64'));
  fseek(fid,1,0);
  hdr.nyquist_zone(rline) = fread(fid,1,'uint8');
  hdr.presums(rline) = fread(fid, 1, 'uint8')+1; % presums are 0-indexed (+1)
  hdr.bit_shifts(rline) = -fread(fid, 1, 'int8');
  hdr.start_idx(rline) = fread(fid, 1, 'uint16');
  hdr.stop_idx(rline) = fread(fid, 1, 'uint16');
  hdr.waveform_ID = fread(fid,1,'uint64');
  
  % Raw data
  hdr.num_sam(rline) = 2*(hdr.stop_idx(rline) - hdr.start_idx(rline));
  
  if ftell(fid) > eof_pos - hdr.num_sam(rline)*SAMPLE_SIZE
    rline = rline - 1;
    break;
  end
  
  if param.records.en
    fseek(fid,hdr.num_sam(rline)*SAMPLE_SIZE,0);
  else
    % Real data
    data(1:hdr.num_sam(rline),rline) = fread(fid,hdr.num_sam(rline),'int16=>single');
  end
end
fclose(fid);

hdr.finfo.syncs = hdr.finfo.syncs(1:rline);
hdr.frame_sync = hdr.frame_sync(1:rline);
hdr.epri = hdr.epri(1:rline);
hdr.seconds = hdr.seconds(1:rline);
hdr.fraction = hdr.fraction(1:rline);
hdr.counter = hdr.counter(1:rline);
hdr.comp_time_sod = hdr.comp_time_sod(1:rline);
hdr.nyquist_zone = hdr.nyquist_zone(1:rline);
hdr.presums = hdr.presums(1:rline);
hdr.bit_shifts = hdr.bit_shifts(1:rline);
hdr.start_idx = hdr.start_idx(1:rline);
hdr.stop_idx = hdr.stop_idx(1:rline);
hdr.num_sam = hdr.num_sam(1:rline);

hdr.seconds = BCD_to_seconds(hdr.seconds);
hdr.utc_time_sod = hdr.seconds + hdr.fraction / param.clk*2;
hdr.Tadc = hdr.start_idx / param.clk*2 - 10.8e-6;

if param.records.en
  %% Remap outputs to match create_task.m output standard:
  %   [hdr] --> [success]
  %   [data] --> hdr
  hdr.offset = hdr.finfo.syncs;
  hdr.wfs.presums = hdr.presums(1);
  data = hdr;
  hdr = 1;
else
  data = data(:,1:rline);
end

return;

% ===============================================================
% ===============================================================
% Example
% ===============================================================
% ===============================================================

fn = '/process-archive/20170226/fmcw/snow/snow4_00_20170226_182432_0103.bin';
[hdr,data] = basic_load_fmcw8(fn,struct('clk',250e6));

