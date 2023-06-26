function [hdr,data] = basic_load_accum2(fn, param)
% [hdr,data] = basic_load_accum2(fn, param)
%
% This is the only function which loads raw data directly.
%
% Loads a single accumulation radar file from the PXI-E Sundance Timothy
% Rink digital system. This function is primarily for debugging.
%
% fn = filename of accumulation radar data
% param = struct controlling loading of data
%   .clk = clock (Hz), default 1 GHz, used to interpret
%     counts in the header fields
%     trigger_delay count is clk/4 (250 MHz)
%     radar_time* count is clk/100 (10 MHz)
%   .recs = 2 element vector for records to load [start_rec num_rec]
%     start_rec uses one-indexing
%   .first_byte = First byte to use. Default is zero. Primarily
%      used for skipping over "garbage" data at the front of the file
%      (scalar integer)
%   .data_out_type = String containing a class type to make data. Typical
%      values are 'int16' and 'single'. Default is 'single'.
%   .all_hdrs = Boolean to load all hdr files and no data.
%
% hdr = file header information for each record unless data output 
%   argument is not specified
% data = cell vector of matrices. Each waveform is an entry in the cell
%   vector. The matrices are Nt by Nx.  First dimension is fast-time,
%   second dimension is slow-time.
%   This output argument is optional.
%   If not supplied, the function reads just the first record's header
%   information and then returns.
%
% Example:
%   hdr = basic_load_accum2('accum2_00_20120204_165304_0000.bin',struct('clk',1e9));
%
%   [hdr,data] = basic_load_accum2('accum2_00_20120204_165304_0000.bin',struct('clk',1e9));
%
%   fn = 'accum2_00_20120204_165304_0000.bin';
%   [hdr,data] = basic_load_accum2(fn,struct('clk',1e9));
%   [hdr,data] = basic_load_accum2(fn,struct('clk',1e9,'data_out_type','int16'));
%
% Authors: John Paden
%
% See also: basic_load_accum, run_basic_load_accum, proc_accum, run_proc_accum,
%  create_vectors_accum, run_create_vectors_accum

SAMP_SIZE = 2;
HEADER_SIZE = 4*16;
FOOTER_SIZE = 4*2;

% ===================================================================
% Check input arguments
% ===================================================================
if ~exist('param','var') || isempty(param)
  param.clk = 1e9;
  param.recs = [];
  param.first_byte = 0;
  param.data_out_type = 'single';
end
if ~isfield(param,'clk');
  param.clk = 1e9;
end
if ~isfield(param,'recs');
  param.recs = [];
end
if ~isfield(param,'first_byte');
  param.first_byte = 0;
end
if ~isfield(param,'data_out_type');
  param.data_out_type = 'single';
end
if ~isfield(param,'all_hdrs');
  param.all_hdrs = false;
end

data_type = str2func(param.data_out_type);

% ===================================================================
% Get the locations of the first sync markers
% ===================================================================
syncs = get_first10_sync_mfile(fn, param.first_byte,struct('sync','1ACFFC1D'));
if isempty(syncs)
  warning('No syncs found in this (empty? corrupt?) file so nothing to return\n');
  hdr = [];
  data = [];
  return;
end
% Find the record size in a robust way
hdr.finfo.rec_size = syncs(2)-syncs(1);%median(diff(syncs));
% Find the first good header
hdr.finfo.syncs(1) = syncs(1);%syncs(find(diff(syncs) == hdr.finfo.rec_size, 1));

% ===================================================================
% Open file
% ===================================================================
[fid,msg] = fopen(fn,'r','ieee-be');
if fid < 1
  fprintf('Could not open file %s\n', fn);
  error(msg);
end

% ===================================================================
% Determine the number of slow time records and which records
% will be read
% ===================================================================
fseek(fid,0,'eof');
hdr.finfo.file_size = ftell(fid);

hdr.finfo.num_rec = floor((hdr.finfo.file_size-hdr.finfo.syncs(1))/hdr.finfo.rec_size);

if isempty(param.recs)
  param.recs(1) = 1;
  param.recs(2) = hdr.finfo.num_rec;
end  
if param.recs(1) + param.recs(2)-1 > hdr.finfo.num_rec
  warning('Only %d records\n', hdr.finfo.num_rec);
end
recs = param.recs(1) + (0:param.recs(2)-1);
num_rec = recs(end)-recs(1)+1;

% Seek to first record
fseek(fid, hdr.finfo.syncs(1) + hdr.finfo.rec_size*(recs(1)-1), 'bof');

% ===================================================================
% Read in file
%   See MCORDS file format.docx
% ===================================================================

% Only header needs to be returned in output argument, so we read
% only the first record
hdr_data = fread(fid,[16 1],'16*uint32',hdr.finfo.rec_size-16*4);

hdr.frame_sync = hdr_data(1,:);
hdr.firmware_version = floor(hdr_data(2,:) / 2^24) ...
  + 0.1*mod(hdr_data(2,:)/2^16,2^8) ...
  + 0.01*mod(hdr_data(3,:)/2^8,2^8) ...
  + 0.001*mod(hdr_data(4,:),2^8);
hdr.radar_time = (hdr_data(3,:)*2^32 + hdr_data(4,:)) / (param.clk/100);
hdr.radar_time_1pps = (hdr_data(5,:)*2^32 + hdr_data(6,:)) / (param.clk/100);
hdr.range_gate_start = hdr_data(7,:);
hdr.range_gate_duration = floor(hdr_data(8,:)/2^16);
hdr.num_coh_ave = mod(floor(hdr_data(8,:)/2^8),2^8);
hdr.profile_num = mod(floor(hdr_data(8,:)/2^1),2^4);
hdr.zero_pi_en = mod(hdr_data(8,:),2);
hdr.trigger_delay = (hdr_data(9,:)*2^32 + hdr_data(10,:)) / (param.clk/4);
if 0
  % Original unofficial method
  hdr.comp_time = datenum_to_epoch(datenum( ...
    mod(floor(hdr_data(12,:)/2^24),2^8)*100 + mod(floor(hdr_data(12,:)/2^16),2^8), ...
    mod(floor(hdr_data(12,:)/2^8),2^8), ...
    mod(floor(hdr_data(12,:)/2^0),2^8), ...
    mod(floor(hdr_data(11,:)/2^16),2^8), ...
    mod(floor(hdr_data(11,:)/2^8),2^8), ...
    mod(floor(hdr_data(11,:)/2^0),2^8)));
else
  hdr.comp_time = hdr_data(12,:);
end
hdr.header_end = hdr_data(16,:);

wf_rec_size = HEADER_SIZE + FOOTER_SIZE + SAMP_SIZE*hdr.range_gate_duration;
hdr.num_wfs = hdr.finfo.rec_size / wf_rec_size;

if nargout == 1 && ~param.all_hdrs
  fclose(fid);
  
  hdr.wfs = [];
  for wf = 1:hdr.num_wfs
    hdr.wfs(wf).num_sam = hdr.range_gate_duration;
    hdr.wfs(wf).which_bits = 0;
    hdr.wfs(wf).t0 = hdr.trigger_delay(1) + hdr.range_gate_start(1)/param.clk - 9.68e-6;
    hdr.wfs(wf).presums = hdr.num_coh_ave(1);
  end
  
  return;
end

% Seek to first record again
fseek(fid, hdr.finfo.syncs(1) + hdr.finfo.rec_size*(recs(1)-1), 'bof');

% Read in all the header data
hdr_data = fread(fid,[16 num_rec],'16*uint32',hdr.finfo.rec_size-16*4);

% Check to make sure frame sync was never lost
hdr.frame_sync = hdr_data(1,:);
if any(hdr.frame_sync ~= hdr.frame_sync(1))
  fclose(fid);
  fprintf('  Loss of frame sync, loading file the slow way\n');
  slow_method = true;
  hdr.finfo = frame_sync_info(fn,struct('sync','1ACFFC1D','cont_mode',0,'rec_size',hdr.num_wfs*(18*4+2*hdr.range_gate_duration)));
  % Remove all short frames
  bad_mask = diff(hdr.finfo.syncs) < hdr.num_wfs*(18*4+2*hdr.range_gate_duration);
  hdr.finfo.syncs = hdr.finfo.syncs(~bad_mask);
  hdr.finfo.num_rec = length(hdr.finfo.syncs);
  [fid,msg] = fopen(fn,'r','ieee-be');
  if param.recs(2) > hdr.finfo.num_rec-param.recs(1)+1
    param.recs(2) = hdr.finfo.num_rec-param.recs(1)+1;
    recs = param.recs(1) + (0:param.recs(2)-1);
    num_rec = recs(end)-recs(1)+1;
  end
  tmp_data = zeros(hdr.finfo.rec_size/2, param.recs(2),'int16');
  for record = 1:size(tmp_data,2)
    fseek(fid,hdr.finfo.syncs(param.recs(1)+record-1),-1);
    tmp_data(:,record) = fread(fid,hdr.finfo.rec_size/2,'int16=>int16');
  end
  fclose(fid);
  
  hdr_data = double(tmp_data(1:HEADER_SIZE/2,:));
  hdr_data = 2^16*hdr_data(1:2:end,:) + hdr_data(2:2:end,:);
else
  slow_method = false;
end

hdr.frame_sync = hdr_data(1,:);
hdr.firmware_version = floor(hdr_data(2,:) / 2^24) ...
  + 0.1*mod(floor(hdr_data(2,:)/2^16),2^8) ...
  + 0.01*mod(floor(hdr_data(2,:)/2^8),2^8) ...
  + 0.001*mod(hdr_data(2,:),2^8);
hdr.radar_time = (hdr_data(3,:)*2^32 + hdr_data(4,:)) / (param.clk/100);
hdr.radar_time_1pps = (hdr_data(5,:)*2^32 + hdr_data(6,:)) / (param.clk/100);
hdr.range_gate_start = hdr_data(7,:);
hdr.range_gate_duration = floor(hdr_data(8,:)/2^16);
hdr.num_coh_ave = mod(floor(hdr_data(8,:)/2^8),2^8);
hdr.profile_num = mod(floor(hdr_data(8,:)/2^1),2^3);
hdr.zero_pi_en = mod(hdr_data(8,:),2);
hdr.trigger_delay = (hdr_data(9,:)*2^32 + hdr_data(10,:)) / (param.clk/4);
if 0
  % Original unofficial method
  hdr.comp_time = datenum_to_epoch(datenum( ...
    mod(floor(hdr_data(12,:)/2^24),2^8)*100 + mod(floor(hdr_data(12,:)/2^16),2^8), ...
    mod(floor(hdr_data(12,:)/2^8),2^8), ...
    mod(floor(hdr_data(12,:)/2^0),2^8), ...
    mod(floor(hdr_data(11,:)/2^16),2^8), ...
    mod(floor(hdr_data(11,:)/2^8),2^8), ...
    mod(floor(hdr_data(11,:)/2^0),2^8)));
else
  hdr.comp_time = hdr_data(12,:);
end
hdr.header_end = hdr_data(16,:);

% ===================================================================
% Read and format data from raw file
% ===================================================================

% DEBUG CHECK:
%num_samp = (hdr.finfo.rec_size-41*4)/hdr.num_wfs/SAMP_SIZE

num_sam = hdr.range_gate_duration(1);

hdr.wfs = [];
for wf = 1:hdr.num_wfs
  hdr.wfs(wf).num_sam = num_sam;
  hdr.wfs(wf).which_bits = 0;
  hdr.wfs(wf).t0 = hdr.trigger_delay(1) + hdr.range_gate_start(1)/param.clk - 9.68e-6;
  hdr.wfs(wf).presums = hdr.num_coh_ave(1);
end

if param.all_hdrs
  if ~slow_method
    fclose(fid);
  end
  return
end

if slow_method
  tmp_data = data_type(tmp_data(HEADER_SIZE/2+1:end-FOOTER_SIZE/2,:));
else
  precision = sprintf('%d*int16=>%s',num_sam,param.data_out_type);
  fseek(fid, hdr.finfo.syncs(1) + hdr.finfo.rec_size*(recs(1)-1) + 16*4, 'bof');
  tmp_data = fread(fid,[num_sam hdr.num_wfs*num_rec],precision,18*4);
  fclose(fid);
end

for wf = 1:hdr.num_wfs
  data{wf} = squeeze(tmp_data(:,wf:hdr.num_wfs:end));
end

return;

