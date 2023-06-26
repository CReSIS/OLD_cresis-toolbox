function [hdr,data] = basic_load_fmcw2(fn,param)
% [hdr,data] = basic_load_fmcw2(fn, param)
%
% This is the only function which loads raw data directly.
%
% THIS DATA FORMAT IS A HACK! TWO WAVEFORMS STORED IN THE FILE, BUT REALLY
% ONLY ONE WAVEFORM.  HAVE TO COMBINE TWO WAVEFORMS INTO ONE. CURRENTLY
% THIS IS DONE BY FORCING THE EFFECTIVE PRESUMS TO BE DOUBLED (ADDING
% TOGETHER WF1 and WF2).
%
% Loads a single fmcw2 radar file. This is primarily for debugging.
% NOTE: 64-bit computer may be essential to load a 256 MB file since it will
% consume 512 MB of memory after loading.
%
% If data is not specified as an output argument, only the header is returned
%
% fn = filename of FMCW2 data
% param = struct controlling loading of data
%   .clk = clock (Hz), default one, used to interpret
%     counts in the header fields
%     snow/kuband use sampling frequency 1e9/8
%   .recs = 2 element vector for records to load [start_rec num_rec]
%     start_rec uses zero-indexing
%   .file_version = default is 1
%     1: no loopback_mode or nyquist_zone setting + bug where 2 waveforms
%        were used because the DAQ could not support 1 waveform
%     2: loopback_mode or nyquist_zone settings added in place of utc_time2
%        + bug where 2 waveforms were used because the DAQ could not
%        support 1 waveform
%     3: no loopback_mode or nyquist_zone setting, 1 waveform
%   .records
%     .en = NOT USED. If true, all
%       headers will be loaded and the output arguments become:
%       hdr --> success flag (returns 1)
%       data --> hdr (all headers)
%     .epri = Expected EPRI
%     .force_all = boolean to force loading all record headers
%   .nohack = boolean that if set to true forces the file to load the file
%      as if there is only one waveform
%   .debug_level = 1 is default, 2 generates plots/print-outs
%
% hdr = file header for each record (if data is not defined, behavior
%   depends on param.records.en variable, if it is false only the first hdr
%   is returned, otherwise all the headers are returned)
% data = (optional output) array of radar data where dimensions are
%   1: fast-time/range-bin
%   2: slow-time/range-line/records
%
% Examples: See bottom of file
%
%   fn = 'snow_20120219_030552_00_0000.bin';
%   param.clk = 1e9/8;
%   [hdr,data] = basic_load_fmcw2(fn,param);
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
  param.nohack = false;
  param.records.en = false;
  param.file_version = 1;
end
if ~isfield(param,'clk');
  param.clk = 1e9/8;
end
if ~isfield(param,'recs');
  param.recs = [];
end
if ~isfield(param,'file_version');
  param.file_version = 1;
end
if ~isfield(param,'debug_level');
  param.debug_level = 1;
end
if ~isfield(param,'records');
  param.records.en = false;
end
if ~isfield(param,'nohack');
  param.nohack = false;
end

% Reset/clear hdr struct
hdr = [];

% ===============================================================
% Get first record position
% ===============================================================
hdr.finfo.syncs = get_first10_sync_mfile(fn,0,struct('sync','BADA55E5'));

% ===============================================================
% Open file big-endian for reading
% ===============================================================
[fid,msg] = fopen(fn,'r','ieee-be');
if fid < 1
  fprintf('Could not open file %s\n', fn);
  error(msg);
end

% ===============================================================
% Read in waveform information + record size
% ===============================================================
HEADER_SIZE = 8*4; % in bytes
NUM_ADC = 1;

% CRAPPY FILE FIX: Because of first record being bad sometimes, we have
% added this loop that tests the contents of the waveforms to see if
% they make sense before continuing
sync_offset = 1;
while sync_offset <= length(hdr.finfo.syncs)
  % Seek to number of waveforms field from first good record
  fseek(fid, HEADER_SIZE+1+hdr.finfo.syncs(sync_offset), -1);
  num_waveforms = fread(fid, 1, 'uint8')+1;
  
  % Seek to beginning of first waveform fields from first good record
  fseek(fid, HEADER_SIZE+hdr.finfo.syncs(sync_offset), -1);
  % rec_size = record size in bytes (initialize to header size and then
  %   add each waveform length in the following loop)
  hdr.finfo.rec_size = HEADER_SIZE;
  % wf_offset: helps keep track of the sample offset into each waveform
  wf_offset = 4;
  
  for wf = 1:num_waveforms
    % Read in waveform header
    hdr.wfs(wf).wf_idx = fread(fid, 1, 'uint8');
    hdr.wfs(wf).num_wfs = fread(fid, 1, 'uint8')+1;
    hdr.wfs(wf).presums = (fread(fid, 1, 'uint8')+1)*2; % presums are 0-indexed (+1) and the hack (*2)
    hdr.wfs(wf).bit_shifts = -fread(fid, 1, 'int8');
    hdr.wfs(wf).start_idx = fread(fid, 1, 'uint16');
    hdr.wfs(wf).Tadc = hdr.wfs(wf).start_idx / param.clk - 10.8e-6;
    hdr.wfs(wf).stop_idx = fread(fid, 1, 'uint16');
    hdr.wfs(wf).num_sam = hdr.wfs(wf).stop_idx - hdr.wfs(wf).start_idx;
    
    % Skip passed waveform data
    fseek(fid, 2*NUM_ADC*hdr.wfs(wf).num_sam, 0);
    
    % Update record size with new waveform
    hdr.finfo.rec_size = hdr.finfo.rec_size + (8 + 2*NUM_ADC*hdr.wfs(wf).num_sam);
    
    % Keep track of the sample offset into each waveform
    hdr.wfs(wf).offset = wf_offset;
    wf_offset = wf_offset + 4 + 1*hdr.wfs(wf).num_sam;
  end
  
  if num_waveforms == 1 || (hdr.wfs(1).num_wfs == hdr.wfs(2).num_wfs ...
      && hdr.wfs(1).presums == hdr.wfs(2).presums ...
      && hdr.wfs(1).bit_shifts == hdr.wfs(2).bit_shifts ...
      && hdr.wfs(1).start_idx == hdr.wfs(2).start_idx ...
      && hdr.wfs(1).Tadc == hdr.wfs(2).Tadc ...
      && hdr.wfs(1).stop_idx == hdr.wfs(2).stop_idx ...
      && hdr.wfs(1).num_sam == hdr.wfs(2).num_sam)
    break;
  end
  warning('  Record %i header was bad', sync_offset);
  sync_offset = sync_offset + 1;
end

% hdr.finfo.file_size: Seek to end of file to get file size in bytes
fseek(fid,0,1);
hdr.finfo.file_size = ftell(fid);

hdr.finfo.num_rec = floor((hdr.finfo.file_size - hdr.finfo.syncs(sync_offset)) / hdr.finfo.rec_size);
% Number of slow time records
if isempty(param.recs)
  param.recs = [0 hdr.finfo.num_rec];
  param_recs_set_by_function = true;
else
  if param.recs(1) + param.recs(2) > hdr.finfo.num_rec
    error('Only %i records in file',hdr.finfo.num_rec);
  end
end

if nargout < 2 && ~param.records.en
  % Seek to first record
  fseek(fid, param.recs(1) * hdr.finfo.rec_size + hdr.finfo.syncs(sync_offset), -1);
  
  % Parse header data
  hdr.frame_sync = fread(fid,1,'uint32');
  hdr.epri = fread(fid,1,'uint32');
  if param.file_version == 4
    % Convert seconds from NMEA ASCII string
    %   64 bits: 00 HH MM SS
    %   ASCII zero is "48"
    hdr_data(5:16,1) = fread(fid,12,'uint16');
    hdr.seconds = ((floor(hdr_data(6,:)/2^8)-48)*10 + mod(hdr_data(6,:),2^8)-48) * 3600 ...
      + ((floor(hdr_data(7,:)/2^8)-48)*10 + mod(hdr_data(7,:),2^8)-48) * 60 ...
      + ((floor(hdr_data(8,:)/2^8)-48)*10 + mod(hdr_data(8,:),2^8)-48);
    
    hdr.fraction = 2^48*hdr_data(9,:) ...
      + 2^32*hdr_data(10,:) ...
      + 2^16*hdr_data(11,:) ...
      + hdr_data(12,:);
    hdr.utc_time_sod = hdr.seconds + hdr.fraction/param.clk;
    
    hdr.comp_time_sod= (2^48*hdr_data(13,:) ...
      + 2^32*hdr_data(14,:) ...
      + 2^16*hdr_data(15,:) ...
      + hdr_data(16,:)) / 1000;
  else
    hdr.utc_time_sod = fread(fid,1,'uint32') + fread(fid,1,'uint32') / param.clk;
    hdr.comp_time_sod = double(fread(fid,1,'uint64'));
    tmp = fread(fid,1,'uint32');
    hdr.utc2_time_sod = tmp + fread(fid,1,'uint32') / param.clk;
  end
  if param.file_version == 2
    % Settings is 32 bit integer as hdr_data(7)
    %   bits 31 down to 0
    %   bits 17:16 --> loop back setting (1-2) --> we convert to 0-1
    %   bits  2:0 --> Nyquist zone (1-4) --> we convert to 0-3
    hdr.loopback = mod(floor(tmp/2^16),2^2) - 1;
    hdr.nyquist_zone = mod(tmp,2^3) - 1;
  end
  hdr.wfs = hdr.wfs(1);
  hdr = rmfield(hdr,'finfo');
  
  
  fclose(fid);
  
  hdr.wfs = hdr.wfs(1);
  return;
  
elseif param.records.en
  loading_failed = false;
  
  if ~param.records.force_all
    %% Quick Read (only looks at first and last record)
    
    % Read in first record
    fseek(fid, param.recs(1) * hdr.finfo.rec_size + hdr.finfo.syncs(sync_offset), -1);
    [hdr_data,samples_read] = fread(fid, HEADER_SIZE/4, 'uint32=>uint32');
    
    % Read in last record (this assumes an ideal file)
    fseek(fid, (param.recs(2)-1) * hdr.finfo.rec_size + hdr.finfo.syncs(sync_offset), -1);
    [last_hdr,samples_read] = fread(fid, HEADER_SIZE/4, 'uint32=>uint32');
    
    % Check to see if headers are both correct in terms of frame sync, EPRI
    % number and settings have not changed.
    if ~(last_hdr(1) == hdr_data(1) ...
        && last_hdr(2) == hdr_data(2) + param.recs(2)-1)
      loading_failed = true;
    end
    
    if ~(last_hdr(7) == hdr_data(7))
      % Settings changed in the file, load all records
      param.records.force_all = true;
    end
  end
  
  if ~loading_failed && ~param.records.force_all
    %% Quick Load using just first and last record
    fclose(fid);
    hdr.epri = hdr_data(2):last_hdr(2);
    if param.file_version == 4
      % Convert seconds from NMEA ASCII string
      %   64 bits: 00 HH MM SS
      %   ASCII zero is "48"
      hdr.seconds = ((mod(floor(hdr_data(3,:)/2^8),2^8)-48)*10 + mod(hdr_data(3,:),2^8)-48) * 3600 ...
        + ((floor(hdr_data(4,:)/2^24)-48)*10 + mod(floor(hdr_data(4,:)/2^16),2^8)-48) * 60 ...
        + ((mod(floor(hdr_data(4,:)/2^8),2^8)-48)*10 + mod(hdr_data(4,:),2^8)-48);
      hdr.seconds = hdr.seconds .* ones(size(hdr.epri),'uint32');
      
      hdr.fraction = 2^32*hdr_data(5,:) ...
        + hdr_data(6,:);
      hdr.fraction = double(hdr.fraction) + (0:param.recs(2)-1) * param.records.epri*param.clk;
    else
      hdr.seconds = hdr_data(3) .* ones(size(hdr.epri),'uint32');
      hdr.fraction = double(hdr_data(4)) + (0:param.recs(2)-1) * param.records.epri*param.clk;
    end
    need_second_jump = find(hdr.fraction/param.clk > 1,1);
    while ~isempty(need_second_jump)
      hdr.fraction(need_second_jump:end) = hdr.fraction(need_second_jump:end) - param.clk;
      hdr.seconds(need_second_jump:end) = hdr.seconds(need_second_jump:end) + 1;
      need_second_jump = find(hdr.fraction/param.clk > 1,1);
    end
    hdr.fraction = uint32(hdr.fraction);
    
    if param.file_version == 2
      % Settings is 32 bit integer as hdr_data(7)
      %   bits 31 down to 0
      %   bits 17:16 --> loop back setting (1-2) --> we convert to 0-1
      %   bits  2:0 --> Nyquist zone (1-4) --> we convert to 0-3
      hdr.loopback_mode = uint8(mod(floor(hdr_data(7,:)/2^16),2^2) - 1);
      hdr.nyquist_zone = uint8(mod(hdr_data(7,:),2^3) - 1);
      hdr.loopback_mode = hdr.loopback_mode .* ones(size(hdr.epri),'uint8');
      hdr.nyquist_zone = hdr.nyquist_zone .* ones(size(hdr.epri),'uint8');
    end
    
    hdr.offset = uint32(hdr.finfo.syncs(sync_offset) + hdr.finfo.rec_size*(0:param.recs(2)-1));
  else
    if ~loading_failed && param.records.force_all
      %% Load all records assuming no file errors
      
      % Read in all records
      fseek(fid, param.recs(1) * hdr.finfo.rec_size + hdr.finfo.syncs(sync_offset), -1);
      [hdr_data,samples_read] = fread(fid, [HEADER_SIZE/4 param.recs(2)], sprintf('%i*uint32=>uint32',HEADER_SIZE/4), hdr.finfo.rec_size-HEADER_SIZE);
      
      fclose(fid);
      
      % Parse header data
      hdr_data = double(hdr_data);
      hdr.frame_sync = hdr_data(1,:);
    end
    
    if loading_failed || any(hdr.frame_sync ~= hdr.frame_sync(1))
      %% Slow Header Load Case
      fprintf('  Loss of frame sync, loading file the slow way\n');
      hdr.finfo = frame_sync_info(fn,struct('sync','BADA55E5','cont_mode',0));
      [fid,msg] = fopen(fn,'r','ieee-be');
      if param.recs(2) > hdr.finfo.num_rec-param.recs(1)
        param.recs(2) = hdr.finfo.num_rec-param.recs(1);
      else
        if param_recs_set_by_function
          param.recs(2) = hdr.finfo.num_rec-param.recs(1);
        end
      end
      hdr_data = zeros(HEADER_SIZE/4, param.recs(2)-param.recs(1),'uint32');
      for record = 1:size(hdr_data,2)
        fseek(fid,hdr.finfo.syncs(param.recs(1)+record),-1);
        hdr_data(:,record) = fread(fid,HEADER_SIZE/4,'int32=>int32');
      end
      fclose(fid);
      
      % Parse header data
      hdr_data = double(hdr_data);
      
      hdr.offset = uint32(hdr.finfo.syncs(param.recs(1)+1:param.recs(2)));
    else
      hdr.offset = uint32(hdr.finfo.syncs(sync_offset) + hdr.finfo.rec_size*(0:param.recs(2)-1));
      hdr = rmfield(hdr,'frame_sync');
    end
    
    hdr.epri = uint32(hdr_data(2,:));
    if param.file_version == 4
      % Convert seconds from NMEA ASCII string
      %   64 bits: 00 HH MM SS
      %   ASCII zero is "48"
      hdr.seconds = ((floor(hdr_data(3,:)/2^8)-48)*10 + mod(hdr_data(3,:),2^8)-48) * 3600 ...
        + ((floor(hdr_data(4,:)/2^24)-48)*10 + mod(floor(hdr_data(4,:)/2^16),2^8)-48) * 60 ...
        + ((floor(hdr_data(4,:)/2^8)-48)*10 + mod(hdr_data(4,:),2^8)-48);
      
      hdr.fraction = 2^32*hdr_data(5,:) ...
        + hdr_data(6,:);
    else
      hdr.seconds = hdr_data(3,:);
      hdr.fraction = hdr_data(4,:);
    end
    if param.file_version == 2
      % Settings is 32 bit integer as hdr_data(7)
      %   bits 31 down to 0
      %   bits 17:16 --> loop back setting (1-2) --> we convert to 0-1
      %   bits  2:0 --> Nyquist zone (1-4) --> we convert to 0-3
      hdr.loopback_mode = uint8(mod(floor(hdr_data(7,:)/2^16),2^2) - 1);
      hdr.nyquist_zone = uint8(mod(hdr_data(7,:),2^3) - 1);
    end
  end
  
  hdr.wfs = hdr.wfs(1);
  hdr = rmfield(hdr,'finfo');
  
  %% Remap outputs to match create_task.m output standard:
  %   [hdr] --> [success]
  %   [data] --> hdr
  data = hdr;
  hdr = 1;
  return;
end

% ===============================================================
% Read in all file data and close file
% ===============================================================

% Seek to first record
fseek(fid, param.recs(1) * hdr.finfo.rec_size + hdr.finfo.syncs(sync_offset), -1);

% Read in all records
[raw_file_data,samples_read] = fread(fid, [hdr.finfo.rec_size/2 param.recs(2)], 'int16=>int16');

% Close file
fclose(fid);

% ===============================================================
% Parse frame sync header
% ===============================================================

% Convert header data from signed integers to unsigned integers
hdr_data = double(raw_file_data(1:HEADER_SIZE/2,:));
hdr_data(hdr_data<0) = 2^16+hdr_data(hdr_data<0);

% Parse header data
hdr.frame_sync = 2^16*hdr_data(1,:) + hdr_data(2,:);

if any(hdr.frame_sync ~= hdr.frame_sync(1))
  fprintf('  Loss of frame sync, loading file the slow way\n');
  hdr.finfo = frame_sync_info(fn,struct('sync','BADA55E5','cont_mode',0));
  [fid,msg] = fopen(fn,'r','ieee-be');
  if param.recs(2) > hdr.finfo.num_rec-param.recs(1)
    param.recs(2) = hdr.finfo.num_rec-param.recs(1);
  else
    if param_recs_set_by_function
      % If param.recs(2) was set by the function, then just readjust it
      % to match the correct file length from frame_sync_info
      param.recs(2) = hdr.finfo.num_rec-param.recs(1);
    end
  end
  raw_file_data = zeros(hdr.finfo.rec_size/2, param.recs(2),'int16');
  for record = 1:size(raw_file_data,2)
    fseek(fid,hdr.finfo.syncs(param.recs(1)+record),-1);
    raw_file_data(:,record) = fread(fid,hdr.finfo.rec_size/2,'int16=>int16');
  end
  fclose(fid);
  
  hdr_data = double(raw_file_data(1:HEADER_SIZE/2,:));
  hdr_data(hdr_data<0) = 2^16+hdr_data(hdr_data<0);
  
  % Parse header data
  hdr.frame_sync = 2^16*hdr_data(1,:) + hdr_data(2,:);
end

% ===============================================================
% Parse remaining header
% ===============================================================
hdr.epri = 2^16*hdr_data(3,:) + hdr_data(4,:);
if param.file_version == 4
  % Convert seconds from NMEA ASCII string
  %   64 bits: 00 HH MM SS
  %   ASCII zero is "48"
  hdr.seconds = ((floor(hdr_data(6,:)/2^8)-48)*10 + mod(hdr_data(6,:),2^8)-48) * 3600 ...
    + ((floor(hdr_data(7,:)/2^8)-48)*10 + mod(hdr_data(7,:),2^8)-48) * 60 ...
    + ((floor(hdr_data(8,:)/2^8)-48)*10 + mod(hdr_data(8,:),2^8)-48);
  
  hdr.fraction = 2^48*hdr_data(9,:) ...
    + 2^32*hdr_data(10,:) ...
    + 2^16*hdr_data(11,:) ...
    + hdr_data(12,:);
  hdr.utc_time_sod = hdr.seconds + hdr.fraction/param.clk;
  
  hdr.comp_time_sod= (2^48*hdr_data(13,:) ...
    + 2^32*hdr_data(14,:) ...
    + 2^16*hdr_data(15,:) ...
    + hdr_data(16,:)) / 1000;
else
  hdr.seconds = 2^16*hdr_data(5,:) + hdr_data(6,:);
  
  hdr.fraction = 2^16*hdr_data(7,:) + hdr_data(8,:);
  
  hdr.utc_time_sod = 2^16*hdr_data(5,:) + hdr_data(6,:) ...
    + (2^16*hdr_data(7,:) + hdr_data(8,:)) / param.clk;
  hdr.comp_time_sod = 2^16*hdr_data(9,:) + hdr_data(10,:) ...
    + (2^16*hdr_data(11,:) + hdr_data(12,:)) / 1000;
  hdr.utc2_time_sod = 2^16*hdr_data(13,:) + hdr_data(14,:) ...
    + (2^16*hdr_data(15,:) + hdr_data(16,:)) / param.clk;
end

if param.file_version == 2
  % Pack settings into 8 bit integer
  %   bits 7 down to 0
  %   bits 2 --> loopback
  %   bits 1:0 --> nyquist zone
  hdr.loopback_mode = uint8(mod(floor(hdr_data(13,:)/2^16),2^2) - 1);
  hdr.nyquist_zone = uint8(mod(hdr_data(14,:),2^3) - 1);
end

% ===============================================================
% Parse data into data matrix
% ===============================================================
param.wfs = 1:length(hdr.wfs);
% Raw file data is 2-D and sorted in this order:
%   1: fast-time (every other record in wf1/wf2), 2: slow-time
% We want to sort in this order
%   1: fast-time, 2: slow-time
if param.nohack
  data = single(raw_file_data(HEADER_SIZE/2+hdr.wfs(1).offset + (1:hdr.wfs(1).num_sam),:));
else
  % This file is the hack version and the 2 waveforms need to be added together
  data = single(raw_file_data(HEADER_SIZE/2+hdr.wfs(1).offset + (1:hdr.wfs(1).num_sam),:) ...
    + raw_file_data(HEADER_SIZE/2+hdr.wfs(2).offset + (1:hdr.wfs(2).num_sam),:));
end
hdr.wfs = hdr.wfs(1);

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
  
  fprintf('Sync, first: %d\n', hdr.finfo.syncs(1));
  fprintf('Sync, spacing: '); fprintf('%d,', diff(hdr.finfo.syncs)); fprintf('\b\n');
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

