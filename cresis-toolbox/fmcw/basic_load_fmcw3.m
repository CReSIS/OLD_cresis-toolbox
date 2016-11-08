function [hdr,data] = basic_load_fmcw3(fn,param)
% [hdr,data] = basic_load_fmcw3(fn, param)
%
% This is the only function which loads raw data directly.
%
% Loads a single fmcw3 radar file. This is primarily for debugging.
% NOTE: 64-bit computer may be essential to load a 256 MB file since it will
% consume 512 MB of memory after loading.
%
% If data is not specified as an output argument, only the header is returned
%
% fn = filename of FMCW3 data
% param = struct controlling loading of data
%   .clk = clock (Hz), default one, used to interpret
%     counts in the header fields
%     snow/kuband use sampling frequency 1e9/8
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
%
%   fn = '/process/20130321/fmcw/kuband/kuband3_02_20130321_190648_0273.bin';
%   [hdr,data] = basic_load_fmcw3(fn,struct('clk',1e9/8));
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
  param.file_version = 3;
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
% 16-bit zeros
% 8-bit presums
% 8-bit bit shift
% 16-bit start
% 16-bit stop
% 16-bit DC offset
% 16-bit NCO value of the DDC
% 8-bit ADC external filter select (Nyquist zone select)
% 8-bit DDC filter select (FIRDEC rate select)
%  File Version 5:
%  0: One half rate filter applied (i.e. IQ DDC followed by x2 decimation)
%  1: Two half rate filters applied (i.e. IQ DDC followed by x4 decimation)
%  2: Three half rate filters applied (i.e. IQ DDC followed by x8 decimation)
%  3: Four half rate filters applied (i.e. IQ DDC followed by x16 decimation)
%  File Version 3:  --> WHEN LOADED IT IS CONVERTED TO FILE VERSION 5 VALUES
%  0: Two half rate filters applied (i.e. IQ DDC followed by x4 decimation)
%  1: Three half rate filters applied (i.e. IQ DDC followed by x8 decimation)
%  2: Four half rate filters applied (i.e. IQ DDC followed by x16 decimation)
%  3: Five half rate filters applied (i.e. IQ DDC followed by x32 decimation)
% 8-bit (Boolean 1 or 0) constant input or ADC input
% 8-bit (Boolean 1 or 0) raw data or ddc data
% DATA: depending on settings
% If raw data then 16 bit real with size stop-start
% If DDC data then 32 bit complex (16-bit I and Q) size (stop-start)/(2^DDC_Filter_Select)
%
% The actual size of the DDC data could vary by 1 depending on how the value is rounded.
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
  hdr.utc_time_sod = hdr.seconds + hdr.fraction / param.clk;
  hdr.counter = fread(fid,1,'uint64');
  hdr.comp_time_sod = double(fread(fid,1,'uint64'));
  fseek(fid,2,0);
  hdr.presums = fread(fid, 1, 'uint8')+1; % presums are 0-indexed (+1)
  hdr.bit_shifts = -fread(fid, 1, 'int8');
  hdr.start_idx = fread(fid, 1, 'uint16');
  hdr.Tadc = hdr.start_idx / param.clk - 10.8e-6;
  hdr.stop_idx = fread(fid, 1, 'uint16');
  hdr.DC_offset = fread(fid,1,'int16');
  hdr.NCO_freq = fread(fid,1,'uint16');
  hdr.nyquist_zone = fread(fid,1,'uint8');
  if param.file_version == 3
    hdr.DDC_filter_select = fread(fid,1,'uint8') + 1;
  else
    hdr.DDC_filter_select = fread(fid,1,'uint8');
  end
  hdr.input_selection = fread(fid,1,'uint8');
  if param.file_version == 3
    hdr.DDC_or_raw_select = fread(fid,1,'uint8');
  else
    hdr.DDC_or_raw_select = fread(fid,1,'uint8');
    if hdr.DDC_or_raw_select == 1
      hdr.DDC_or_raw_select = 0;
      hdr.DDC_filter_select = -1;
    end
  end
  
  fclose(fid);
  return;
elseif param.records.en
  % Seek to first record
  fseek(fid, hdr.finfo.syncs(1) + 36, -1);
  hdr.start_idx = fread(fid, 1, 'uint16');
  hdr.stop_idx = fread(fid, 1, 'uint16');
  fseek(fid,5,0);
  if param.file_version == 3
    hdr.DDC_filter_select = fread(fid,1,'uint8') + 1;
  else
    hdr.DDC_filter_select = fread(fid,1,'uint8');
  end
  fseek(fid,1,0);
  if param.file_version == 3
    hdr.DDC_or_raw_select = fread(fid,1,'uint8');
  else
    hdr.DDC_or_raw_select = fread(fid,1,'uint8');
    if hdr.DDC_or_raw_select == 1
      hdr.DDC_or_raw_select = 0;
      hdr.DDC_filter_select = -1;
    end
  end
  
  rline = 1;
  if hdr.DDC_or_raw_select(rline)
    % Raw data
    hdr.num_sam(rline) = hdr.stop_idx(rline) - hdr.start_idx(rline);
  else
    % DDC data
    hdr.num_sam(rline) = floor((hdr.stop_idx(rline) - hdr.start_idx(rline)) ...
      ./ 2.^(hdr.DDC_filter_select(rline) + 1));
  end
  
  fseek(fid, hdr.finfo.syncs(1), -1);
  rec_size = HEADER_SIZE+SAMPLE_SIZE*hdr.num_sam*(1 + ~hdr.DDC_or_raw_select(rline));
  num_rec = floor( (eof_pos-ftell(fid)) / rec_size );
  hdr_data = fread(fid, [12 num_rec], '12*uint32', rec_size-12*4);
  if all(hdr_data(1,:) == hdr_data(1,1)) ...
      && all(hdr_data(10,:) == hdr_data(10,1)) ...
      && all(mod(floor(hdr_data(12,:)/2^16),2^8) == mod(floor(hdr_data(12,1)/2^16),2^8)) ...
      && all(mod(hdr_data(12,:),2^8) == mod(hdr_data(12,1),2^8))
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
end

if 0

  fid = fopen(fn,'r','ieee-be');
  
  data = fread(fid,'uint32').';
  
  file_size = ftell(fid);
  
  fclose(fid);
  
  offset = find(data == param.frame_sync);
  offset = offset(1:end-1);
  
  % 
  
  for field_idx = 1:length(param.field_offsets)
    varargout{field_idx} = data(offset + param.field_offsets(field_idx));
  end
  
  offset = (offset-1) *4;
  
      [file_size offset epri seconds fraction hdr9 hdr10 hdr11] = basic_load_hdr(fn, hdr_param);
    start_idx = floor(hdr9/2^16);
    stop_idx = mod(hdr9,2^16);
    NCO_freq_step = mod(hdr10,2^16);
    nyquist_zone = mod(floor(hdr11/2^24),2^8);
    if param.file_version == 3
      DDC_filter_select = mod(floor(hdr11/2^16),2^8) + 1;
    else
      DDC_filter_select = mod(floor(hdr11/2^16),2^8);
    end
    if param.file_version == 3
      DDC_or_raw_select = mod(hdr11,2^8);
    else
      DDC_or_raw_select = mod(hdr11,2^8);
      DDC_filter_select(DDC_or_raw_select == 1) = DDC_filter_select(DDC_or_raw_select == 1) - 1;
      DDC_or_raw_select(DDC_or_raw_select == 1) = 0;
    end
    if DDC_or_raw_select
      % Raw data
      num_sam = stop_idx - start_idx;
    else
      % DDC data
      num_sam = floor(((stop_idx - start_idx) ./ 2.^(1+DDC_filter_select)));
    end
    HEADER_SIZE = 48;
    SAMPLE_SIZE = 2;
    expected_rec_size = HEADER_SIZE + SAMPLE_SIZE*num_sam.*(1+(DDC_or_raw_select==0));
    meas_rec_size = diff(offset);
    bad_mask = meas_rec_size ~= expected_rec_size(1:end-1);
    bad_mask(end+1) = file_size < offset(end) + expected_rec_size(end);
    if any(bad_mask)
      warning('Found %d record size errors', sum(bad_mask));
    end
    offset = offset(~bad_mask);
    epri = epri(~bad_mask);
    seconds = seconds(~bad_mask);
    fraction = fraction(~bad_mask);
    
    seconds = BCD_to_seconds(seconds);
  
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
  fseek(fid,2,0);
  hdr.presums(rline) = fread(fid, 1, 'uint8')+1; % presums are 0-indexed (+1)
  hdr.bit_shifts(rline) = -fread(fid, 1, 'int8');
  hdr.start_idx(rline) = fread(fid, 1, 'uint16');
  hdr.stop_idx(rline) = fread(fid, 1, 'uint16');
  hdr.DC_offset(rline) = fread(fid,1,'int16');
  hdr.NCO_freq(rline) = fread(fid,1,'uint16');
  hdr.nyquist_zone(rline) = fread(fid,1,'uint8');
  if param.file_version == 3
    hdr.DDC_filter_select(rline) = fread(fid,1,'uint8') + 1;
  else
    hdr.DDC_filter_select(rline) = fread(fid,1,'uint8');
  end
  hdr.input_selection(rline) = fread(fid,1,'uint8');
  if param.file_version == 3
    hdr.DDC_or_raw_select(rline) = fread(fid,1,'uint8');
  else
    hdr.DDC_or_raw_select(rline) = fread(fid,1,'uint8');
    if hdr.DDC_or_raw_select(rline) == 1
      hdr.DDC_or_raw_select(rline) = 0;
      hdr.DDC_filter_select(rline) = -1;
    end
  end
  
  if hdr.DDC_or_raw_select(rline)
    % Raw data
    hdr.num_sam(rline) = hdr.stop_idx(rline) - hdr.start_idx(rline);
  else
    % DDC data
    hdr.num_sam(rline) = floor((hdr.stop_idx(rline) - hdr.start_idx(rline)) ...
      ./ 2.^(hdr.DDC_filter_select(rline) + 1));
  end
  
  if rline < 2 || hdr.num_sam(rline) ~= hdr.num_sam(rline-1)
    % Update record size with new waveform
    hdr.finfo.rec_size(end+1) = HEADER_SIZE + ...
      (1 + ~hdr.DDC_or_raw_select(rline)) * SAMPLE_SIZE * hdr.num_sam(rline);
    % Preallocated records
    num_rec = floor((eof_pos - (ftell(fid)+hdr.num_sam(rline)*SAMPLE_SIZE*(1 + ~hdr.DDC_or_raw_select(rline)))) / (HEADER_SIZE + SAMPLE_SIZE*(1 + ~hdr.DDC_or_raw_select(rline))*hdr.num_sam(rline)));
    data(1,rline+num_rec) = 0;
  end
  
  if ftell(fid) > eof_pos - hdr.num_sam(rline)*SAMPLE_SIZE*(1 + ~hdr.DDC_or_raw_select(rline))
    rline = rline - 1;
    break;
  end
  
  if param.records.en
    fseek(fid,hdr.num_sam(rline)*SAMPLE_SIZE*(1 + ~hdr.DDC_or_raw_select(rline)),0);
  else
    if hdr.DDC_or_raw_select(rline)
      % Real data
      data(1:hdr.num_sam(rline),rline) = fread(fid,hdr.num_sam(rline),'int16=>single');
      % The following line is not necessary with some data files (e.g. 2015 Alaska TO_NRL)
      data(:,rline) = data(reshape([2:2:hdr.num_sam(rline);1:2:hdr.num_sam(rline)-1],[hdr.num_sam(rline) 1]),rline);
    else
      % Complex data
      tmp = fread(fid,2*hdr.num_sam(rline),'int16=>single');
      data(1:hdr.num_sam(rline),rline) = tmp(1:2:end) + 1i*tmp(2:2:end);
    end
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
hdr.presums = hdr.presums(1:rline);
hdr.bit_shifts = hdr.bit_shifts(1:rline);
hdr.start_idx = hdr.start_idx(1:rline);
hdr.stop_idx = hdr.stop_idx(1:rline);
hdr.num_sam = hdr.num_sam(1:rline);
hdr.DC_offset = hdr.DC_offset(1:rline);
hdr.NCO_freq = hdr.NCO_freq(1:rline);
hdr.nyquist_zone = hdr.nyquist_zone(1:rline);
hdr.DDC_filter_select = hdr.DDC_filter_select(1:rline);
hdr.input_selection = hdr.input_selection(1:rline);
hdr.DDC_or_raw_select = hdr.DDC_or_raw_select(1:rline);

hdr.seconds = BCD_to_seconds(hdr.seconds);
hdr.utc_time_sod = hdr.seconds + hdr.fraction / param.clk;
hdr.Tadc = hdr.start_idx / param.clk - 10.8e-6;

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

fn = '/process/20130321/fmcw/kuband/kuband3_02_20130321_190648_0273.bin';
[hdr,data] = basic_load_fmcw3(fn,struct('clk',1e9/8));
[hdr,data] = basic_load_fmcw3(fn,struct('clk',125e6,'recs',[0 500]));

fn = '/process/20130321/fmcw/kuband/kuband3_02_20130321_185201_0257.bin';
[hdr,data] = basic_load_fmcw3(fn,struct('clk',1e9/8));

