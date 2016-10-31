function [hdr,data] = basic_load_mcords4(fn,param)
% [hdr,data] = basic_load_mcords4(fn, param)
%
% This is the only function which loads raw data directly.
% This is for the MCORDS4 (MRI UWB Basler RDS system, 150-600 MHz)
%
% Loads a single mcords4 radar file. This is use for debugging and for
% records file generation.
% NOTE: 64-bit computer may be essential to load a 256 MB file since it will
% consume 512 MB of memory after loading.
%
% fn = filename of MCoRDS-4 data
% param = struct controlling loading of data
%   .clk = clock (Hz), default one, used to interpret
%     counts in the header fields. This defaults to the sampling clock
%     divided by 4.
%   .recs = 2 element vector for records to load [start_rec num_rec]
%     start_rec uses zero-indexing, default is [0 inf]
%   .debug_level = 1 is default, 2 generates plots/print-outs
%   .first_byte = first byte to start reading at (default is zero0)
%   .start_index_time_offset = time offset between transmit waveform start
%     and the start index (default is -10.8e-6 implying that the first
%     sample collected is 10.8 us before the transmit event starts)
%
% hdr = file header for each record
% data = cell vector of single matrices of radar data where each entry
%   in the cell vector is a 3-D array for that waveform. Dimensions
%   1: fast-time/range-bin
%   2: slow-time/range-line/records
%
% Examples: See bottom of file
%
%   fn = 'UWB_Ch0_20130829_171826_00_0001.bin';
%   param.clk = 500e6;
%   [hdr,data] = basic_load_mcords4(fn,param);
%
% Authors: John Paden

% ===================================================================
%% Check input arguments
% ===================================================================
if ~exist('param','var')
  param = struct();
end
if ~isfield(param,'clk');
  param.clk = 1e9/2 / 4; % Default to sampling clock divided by 4
end
if ~isfield(param,'recs');
  param.recs = [0 inf];
end
if ~isfield(param,'debug_level');
  param.debug_level = 1;
end
if ~isfield(param,'first_byte');
  param.first_byte = 0;
end
if ~isfield(param,'start_index_time_offset')
  param.start_index_time_offset = -1.0665e-05;
end
if ~isfield(param,'header_only')
  if nargout <= 1
    param.header_only = 1;
  else
    param.header_only = 0;
  end
end
if ~isfield(param,'record_mode')
  param.record_mode = 1;
end

% Reset/clear hdr struct
hdr = [];

% ===============================================================
%% Get first record position
% ===============================================================
hdr.sync_offsets = get_first10_sync_mfile(fn,0,struct('sync','1ACFFC1D','first_byte',param.first_byte));

% ===============================================================
%% Open file big-endian for reading
% ===============================================================
[fid,msg] = fopen(fn,'r','ieee-be');
if fid < 1
  fprintf('Could not open file %s\n', fn);
  error(msg);
end

% ===============================================================
%% Read in waveform information + record size
% ===============================================================
HEADER_SIZE = 40;
WF_HEADER_SIZE = 8;

fseek(fid, HEADER_SIZE+1+hdr.sync_offsets(1), -1);
num_waveforms = fread(fid, 1, 'uint8') + 1;

% Seek to beginning of first waveform
fseek(fid, HEADER_SIZE+hdr.sync_offsets(1), -1);
% rec_size = record size in 16-bit words
hdr.rec_size = HEADER_SIZE/2;
% wf_offset: helps keep track of the sample offset into each waveform
wf_offset = 4;

for wf = 1:num_waveforms
  % Read in waveform header
  fseek(fid,2,0);
  hdr.wfs(wf).presums = fread(fid, 1, 'uint8');
  hdr.wfs(wf).bit_shifts = -fread(fid, 1, 'int8');
  hdr.wfs(wf).start_idx = fread(fid, 1, 'uint16');
  hdr.wfs(wf).t0 = hdr.wfs(wf).start_idx / param.clk + param.start_index_time_offset;
  hdr.wfs(wf).num_sam = 4*fread(fid, 1, 'uint16');
  
  % Skip passed waveform data
  fseek(fid, 2*hdr.wfs(wf).num_sam, 0);
  
  % Update record size with new waveform
  hdr.rec_size = hdr.rec_size + WF_HEADER_SIZE/2 + hdr.wfs(wf).num_sam;
  
  % Keep track of the sample offset into each waveform
  hdr.wfs(wf).offset = wf_offset;
  wf_offset = wf_offset + WF_HEADER_SIZE/2 + hdr.wfs(wf).num_sam;
end

fseek(fid,0,1);
hdr.file_size = ftell(fid);

if hdr.rec_size ~= median(diff(hdr.sync_offsets))/2;
  warning('Estimated header size is wrong');
  keyboard;
  % For badly recorded files (e.g. DDS settings not matching ADC settings)
  % you can try uncommenting the following line:
  % hdr.rec_size = median(diff(hdr.sync_offsets))/2;
end

if nargout < 2
  % Seek to first record
  fseek(fid, param.recs(1) * hdr.rec_size*2 + hdr.sync_offsets(1), -1);

  % Parse header data
  hdr.frame_sync = fread(fid,1,'uint32');
  hdr.epri = fread(fid,1,'uint32');
  hdr.firmware_version = fread(fid,1,'uint32');
  hdr.measurement_type = fread(fid,1,'uint32');
  if 1
    hdr.seconds = fread(fid,1,'uint32').'; % From NMEA string converted to DCB
    hdr.seconds = BCD_to_seconds(hdr.seconds);
    hdr.fractions = fread(fid,1,'uint32');
    hdr.counter = fread(fid,1,'uint64');
  else
    ascii_HHMMSS = fread(fid,8,'char=>char').'; % From NMEA string
    hdr.seconds = str2double(ascii_HHMMSS(3:4))*3600 + str2double(ascii_HHMMSS(5:6))*60 + str2double(ascii_HHMMSS(7:8));
    hdr.fractions = fread(fid,1,'uint64');
  end
  hdr.utc_time_sod = hdr.seconds + hdr.fractions / param.clk;
  hdr.comp_time_sod = double(fread(fid,1,'uint64'));
  
  fclose(fid);
  
  return;
end

% ===============================================================
% Read in all file data and close file
% ===============================================================

% Seek to first record
fseek(fid, param.recs(1) * hdr.rec_size*2 + hdr.sync_offsets(1), -1);

% Read in all records
raw_file_data = fread(fid, [hdr.rec_size param.recs(2)], 'int16=>int16');

% Close file
fclose(fid);

% ===============================================================
% Parse frame sync header
% ===============================================================

% Convert header data from signed integers to unsigned integers
hdr_data = double(raw_file_data(1:HEADER_SIZE/2,:));
hdr_data(hdr_data<0) = 2^16+hdr_data(hdr_data<0);

% Parse header data
hdr.frame_sync = 2^16*hdr_data(1:hdr.rec_size:end,:) ...
  + hdr_data(2:hdr.rec_size:end,:);

if any(hdr.frame_sync ~= hdr.frame_sync(1))
  fprintf('  Loss of frame sync, loading file the slow way\n');
  hdr.finfo = frame_sync_info(fn,struct('sync','1ACFFC1D','cont_mode',0));
  [fid,msg] = fopen(fn,'r','ieee-be');
  if param.recs(2) > hdr.finfo.num_rec-param.recs(1)
    param.recs(2) = hdr.finfo.num_rec-param.recs(1);
  end
  raw_file_data = zeros(hdr.rec_size, param.recs(2),'int16');
  for record = 1:size(raw_file_data,2)
    fseek(fid,hdr.finfo.syncs(param.recs(1)+record),-1);
    [tmp,count] = fread(fid,hdr.rec_size,'int16=>int16');
    if count ~= hdr.rec_size
      fprintf('%d requested to read, only %d read\n', hdr.rec_size, count);
      keyboard
    end
    raw_file_data(:,record) = tmp;
  end
  fclose(fid);
  
  hdr_data = double(raw_file_data(1:HEADER_SIZE/2,:));
  hdr_data(hdr_data<0) = 2^16+hdr_data(hdr_data<0);

  % Parse header data
  hdr.frame_sync = 2^16*hdr_data(1:hdr.rec_size:end,:) ...
    + hdr_data(2:hdr.rec_size:end,:);
end

% ===============================================================
% Parse remaining header
% ===============================================================
hdr.epri = 2^16*hdr_data(3,:) ...
  + hdr_data(4,:);

% Convert seconds from NMEA ASCII string converted to Decimal Encoded
% Binary
%   32 bits: 0 0 H H M M S S
hdr.seconds = BCD_to_seconds(double(hdr_data(9,:))*2^16 + double(hdr_data(10,:)));

hdr.fractions = 2^16*hdr_data(11,:) ...
  + hdr_data(12,:);

hdr.counter = 2^48*hdr_data(13,:) ...
  + 2^32*hdr_data(14,:) ...
  + 2^16*hdr_data(15,:) ...
  + hdr_data(16,:);
  
hdr.utc_time_sod = hdr.seconds + hdr.fractions/param.clk;

hdr.comp_time_sod= (2^48*hdr_data(17,:) ...
  + 2^32*hdr_data(18,:) ...
  + 2^16*hdr_data(19,:) ...
  + hdr_data(20,:)) / 1000;

% ===============================================================
% Parse data into data matrix
% ===============================================================
param.wfs = 1:length(hdr.wfs);
for wf = 1:length(param.wfs)
  % Raw file data is 2-D and sorted in this order:
  %   1: fast-time, 2: slow-time
  % We want to sort in this order with each having its own dimension
  %   1: fast-time, 2: slow-time
  data{wf} = raw_file_data(HEADER_SIZE/2+hdr.wfs(wf).offset + (1:hdr.wfs(wf).num_sam),:);
  data{wf} = single(data{wf});
end

% ===============================================================
% Optional debug work
% ===============================================================
if param.debug_level >= 2
  for wf_idx=1:length(data)
    for adc_idx = 1:size(data{wf_idx},3)
      figure((wf_idx-1)*10 + adc_idx);
      imagesc(data{wf_idx}(:,:,adc_idx));
      colorbar;
      xlabel('Range line');
      ylabel('Range bin');
    end
  end

  figure(100);
  plot(hdr.frame_sync);
  title('Frame Sync');

  figure(101);
  plot(diff(hdr.epri));
  title('Diff of EPRI');
  fprintf('First EPRI = %.0f\n', hdr.epri(1));

  figure(102);
  plot(diff(hdr.utc_time_sod),'.');
  title('Diff of UTC time');
  PRI = mean(diff(hdr.utc_time_sod));
  fprintf('UTC Time = %.9f sec to %.9f sec\n', hdr.utc_time_sod(1), ...
    hdr.utc_time_sod(end));

  [tmp fn_name] = fileparts(fn);
  utc_time = datenum(str2double(fn_name(11:14)), str2double(fn_name(15:16)), ...
    str2double(fn_name(17:18)),0,0,hdr.utc_time_sod(1));
  fprintf('UTC Time = %s\n', datestr(utc_time));
  fprintf('PRF = %.3f Hz\n', 1/PRI);

  figure(103);
  plot(hdr.comp_time_sod);
  title('Computer time SOD');
  comp_time = datenum(str2double(fn_name(11:14)), str2double(fn_name(15:16)), ...
    str2double(fn_name(17:18)),0,0,hdr.comp_time_sod(1));
  fprintf('Comp Time = %s\n', datestr(comp_time));
  
  fprintf('Sync, first: %d\n', hdr.sync_offsets(1));
  fprintf('Sync, spacing: '); fprintf('%d,', diff(hdr.sync_offsets)); fprintf('\b\n');
end

return;

% ===============================================================
% ===============================================================
% Example
% ===============================================================
% ===============================================================

fn = '/cresis/scratch1/paden/mcords2/board_0/mcords_0_03012011_171839_01_0000.bin';
param.debug_level = 2;
param.clk = 125e6;
%param.recs = [0 inf];
[hdr,data] = basic_load_mcords2(fn,param);

fn = '/cresis/scratch1/paden/mcords2/board_1/mcords_1_03012011_171839_01_0000.bin';
param.debug_level = 2;
param.clk = 125e6;
%param.recs = [0 inf];
[hdr,data] = basic_load_mcords2(fn,param);

fn = '/cresis/scratch1/paden/mcords2/board_0/mcords_0_03012011_171946_01_0009.bin';
param.debug_level = 2;
param.clk = 125e6;
%param.recs = [0 inf];
[hdr,data] = basic_load_mcords2(fn,param);

fn = '/cresis/scratch1/paden/mcords2/board_1/mcords_1_03012011_171946_01_0009.bin';
param.debug_level = 2;
param.clk = 125e6;
%param.recs = [0 inf];
[hdr,data] = basic_load_mcords2(fn,param);

