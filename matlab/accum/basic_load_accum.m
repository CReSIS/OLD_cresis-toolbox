function [hdr,data] = basic_load_accum(fn, param)
% [hdr,data] = basic_load_accum(fn, param)
%
% This is the only function which loads raw data directly.
%
% Loads a single accumulation radar file. This is primarily for debugging.
% NOTE: 64-bit computer is required to load a 512 MB file since it will
% consume 1 GB of memory after loading.
%
% fn = filename of accumulation radar data
% param = struct controlling loading of data
%   .clk = clock (Hz), default one, used to interpret
%     counts in the header fields
%     2011 Greenland accum radar uses sampling frequency 1e9/8
%   .num_wfs = number of waveforms
%   .recs = 2 element vector for records to load [start_rec num_rec]
%     start_rec uses one-indexing
%   .first_byte = First byte to use. Default is zero. Primarily
%      used for skipping over "garbage" data at the front of the file
%      (scalar integer)
%   .data_out_type = String containing a class type to make data. Typical
%      values are 'uint16' and 'single'. Default is 'single'.
%
% hdr = file header information for each record unless data output 
%   argument is not specified
% data = single matrix of radar data. First dimension is fast-time
%   second dimension is slow-time, and third dimension is the 16 channels.
%   This output argument is optional.
%   If not supplied, the function reads just the first record's header
%   information and then returns.
%
% Example:
%   hdr = basic_load_accum('accum_20110418_13031573_0205.dat',struct('clk',1e9/8));
%
%   [hdr,data] = basic_load_accum('accum_20110418_13031573_0205.dat',struct('clk',1e9/8));
%
%   fn = '/cresis/data3/Accum_Data/2011_Greenland_P3/20110329/seg0/data.03292011.0235.dat';
%   [hdr,data] = basic_load_accum(fn,struct('clk',1e9/8));
%   [hdr,data] = basic_load_accum(fn,struct('clk',1e9/8,'data_out_type','uint16'));
%
% Authors: John Paden
%
% See also: run_basic_load_accum, proc_accum, run_proc_accum,
%  create_vectors_accum, run_create_vectors_accum

SAMP_SIZE = 2;
HEADER_SIZE = 160;

% ===================================================================
% Check input arguments
% ===================================================================
if ~exist('param','var') || isempty(param)
  param.clk = 1;
  param.wfs = [];
  param.recs = [];
  param.first_byte = 0;
  param.data_out_type = 'single';
end
if ~isfield(param,'clk');
  param.clk = 1;
end
if ~isfield(param,'num_wfs');
  param.num_wfs = 16;
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
if ~isfield(param,'hdr');
  param.hdr = false;
end

data_type = str2func(param.data_out_type);

% ===================================================================
% Get the locations of the first sync markers
% ===================================================================
syncs = get_first10_sync_mfile(fn, param.first_byte);
% Find the record size in a robust way
hdr.finfo.rec_size = median(diff(syncs));
% Find the first good header
hdr.finfo.syncs(1) = syncs(find(diff(syncs) == hdr.finfo.rec_size, 1));

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
  param.recs(2) = hdr.finfo.num_rec - param.recs(1) + 1;
end
recs = param.recs(1) + (0:param.recs(2)-1);
num_rec = recs(end)-recs(1)+1;

% Seek to first record
fseek(fid, hdr.finfo.syncs(1) + hdr.finfo.rec_size*(recs(1)-1), 'bof');

% ===================================================================
% Read in file
%   See MCORDS file format.docx
% ===================================================================

if nargout == 1
  % Only header needs to be returned in output argument, so we read
  % only the first record
  hdr_data = fread(fid,[40 1],'40*uint32',hdr.finfo.rec_size-40*4);
  fclose(fid);
  
  hdr.frame_sync = hdr_data(1,:);
  hdr.unknown1 = hdr_data(2,:);
  hdr.seconds = hdr_data(3,:);
  hdr.fraction = hdr_data(4,:);
  hdr.unknown2 = hdr_data(5,:);
  hdr.unknown3 = hdr_data(6,:);
  hdr.unknown4 = hdr_data(7,:);
  hdr.unknown5 = hdr_data(8,:);
  hdr.utc_time_sod = hdr.seconds + hdr.fraction/(param.clk/2);
  
  hdr.wfs = [];
  for wf = 1:param.num_wfs
    hdr.wfs(wf).num_samp = hdr_data(9+2*(wf-1),1);
    %hdr.wfs(wf).unknown = hdr_data(10+2*(wf-1),1);
    hdr.wfs(wf).bit_shifts = mod(floor(hdr_data(10+2*(wf-1),1) / 2^24), 2^4);
    hdr.wfs(wf).t0 = mod(floor(hdr_data(10+2*(wf-1),1) / 2^10), 2^14) / param.clk - 8.3450e-06;
    hdr.wfs(wf).presums = 1 + mod(hdr_data(10+2*(wf-1),1), 2^10);
  end
  
  return;
end

% Read in all the header data
hdr_data = fread(fid,[40 num_rec],'40*uint32',hdr.finfo.rec_size-40*4);

% Check to make sure frame sync was never lost
hdr.frame_sync = hdr_data(1,:);
if any(hdr.frame_sync ~= hdr.frame_sync(1))
  fclose(fid);
  fprintf('  Loss of frame sync, loading file the slow way\n');
  slow_method = true;
  hdr.finfo = frame_sync_info(fn,struct('sync','DEADBEEF','cont_mode',0));
  [fid,msg] = fopen(fn,'r','ieee-be');
  if param.recs(2) > hdr.finfo.num_rec-param.recs(1)
    param.recs(2) = hdr.finfo.num_rec-param.recs(1);
    recs = param.recs(1) + (0:param.recs(2)-1);
    num_rec = recs(end)-recs(1)+1;
  end
  data = zeros(hdr.finfo.rec_size/2, param.recs(2),'uint16');
  for record = 1:size(data,2)
    fseek(fid,hdr.finfo.syncs(param.recs(1)+record),-1);
    data(:,record) = fread(fid,hdr.finfo.rec_size/2,'uint16=>uint16');
  end
  fclose(fid);
  
  hdr_data = double(data(1:HEADER_SIZE/2,:));
  hdr_data = 2^16*hdr_data(1:2:end,:) + hdr_data(2:2:end,:);
else
  slow_method = false;
end

hdr.frame_sync = hdr_data(1,:);
hdr.version = hdr_data(2,:);
hdr.seconds = hdr_data(3,:);
hdr.fraction = hdr_data(4,:);
hdr.unknown2 = hdr_data(5,:);
hdr.unknown3 = hdr_data(6,:);
hdr.unknown4 = hdr_data(7,:);
hdr.unknown5 = hdr_data(8,:);
hdr.utc_time_sod = hdr.seconds + hdr.fraction/(param.clk/2);

hdr.wfs = [];
for wf = 1:param.num_wfs
  hdr.wfs(wf).num_samp = hdr_data(9+2*(wf-1),1);
  %hdr.wfs(wf).unknown = hdr_data(10+2*(wf-1),1);
  hdr.wfs(wf).bit_shifts = mod(floor(hdr_data(10+2*(wf-1),1) / 2^24), 2^4);
  hdr.wfs(wf).t0 = mod(floor(hdr_data(10+2*(wf-1),1) / 2^10), 2^14) / param.clk - 8.3450e-06;
  hdr.wfs(wf).presums = 1 + mod(hdr_data(10+2*(wf-1),1), 2^10);
end

% ===================================================================
% Read and format data from raw file
% ===================================================================

% DEBUG CHECK:
%num_samp = (hdr.finfo.rec_size-41*4)/param.num_wfs/SAMP_SIZE

num_samp = hdr.wfs(1).num_samp;
num_samp_per_rec = param.num_wfs*num_samp;

if slow_method
  data = data_type(data(HEADER_SIZE/2+3:end,:));
else
  precision = sprintf('%d*uint16=>%s',num_samp_per_rec,param.data_out_type);
  fseek(fid, hdr.finfo.syncs(1) + hdr.finfo.rec_size*(recs(1)-1) + 81*2, 'bof');
  data = fread(fid,[num_samp_per_rec num_rec],precision,hdr.finfo.rec_size-num_samp_per_rec*SAMP_SIZE);
  fclose(fid);
end

data = reshape(data,[num_samp param.num_wfs size(data,2)]);
% Trim off last 4 values (not valid data)
num_samp = num_samp - 4;
data = data(1:num_samp,:,:);
data = permute(data,[1 3 2]);

return;























