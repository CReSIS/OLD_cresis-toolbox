function [hdr,data] = basic_load_tidsor(fn, param)
% [hdr,data] = basic_load_tidsor(fn, param)
%
% This is the only function which loads raw data directly.
%
% Loads a single TIDSoR radar file. This is primarily for debugging.
% NOTE: 64-bit computer is required to load a 512 MB file since it will
% consume 1 GB of memory after loading. Use cohererent averaging to
% prevent memory allocation failure.
%
% fn = filename of accumulation radar data
% param = struct controlling loading of data
%   .clk = clock (Hz), default one, used to interpret
%     counts in the header fields
%     2011 Greenland accum radar uses sampling frequency 1e9/8
%   .wfs = which waveforms to load (ignored, currently all wf loaded)
%   .recs = 2 element vector for records to load [start_rec num_rec]
%     start_rec uses one-indexing
%   .coh_ave = slow-time coherent averages or presums (default is 1)
%     These are done as the data is loaded to reduce memory usage.
%   .first_byte = First byte to use. Default is zero. Primarily
%      used for skipping over "garbage" data at the front of the file
%      (scalar integer)
%
% hdr = file header information for each record unless data output 
%   argument is not specified
% data = single matrix of radar data. First dimension is fast-time
%   second dimension is slow-time, and third dimension is the waveforms.
%   This output argument is optional.
%   If not supplied, the function reads just the first record's header
%   information and then returns.
%
% Example:
%   fn = '/cresis/data2/TIDSoR/2011_Chile/Fligth_test_Sat_Jul23_2011/data.07232011.0001.dat';
%   [hdr,data] = basic_load_tidsor(fn,struct('clk',1e9/8,'coh_ave',100));
%
%   fn = '/cresis/data2/TIDSoR/2011_Chile/Fligth_test_Sat_Jul23_2011/data.07232011.0001.dat';
%   [hdr,data] = basic_load_tidsor(fn,struct('clk',1e9/8,'recs',[1 100]));
%
%   fn = '/cresis/data2/TIDSoR/2011_Chile/Fligth_test_Sat_Jul23_2011/data.07232011.0000.dat';
%   [hdr,data] = basic_load_tidsor(fn,struct('clk',1e9/8));
%
% Authors: John Paden
%
% See also: run_basic_load_tidsor

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
end
if ~isfield(param,'clk')
  param.clk = 1;
end
if ~isfield(param,'wfs')
  param.wfs = [];
end
if ~isfield(param,'recs')
  param.recs = [];
end
if ~isfield(param,'first_byte')
  param.first_byte = 0;
end
if ~isfield(param,'coh_ave') || isempty(param.coh_ave)
  param.coh_ave = 1;
end

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
end
% Round number of records to nearest multiple of the coherent averages
param.recs(2) = floor(param.recs(2)/param.coh_ave) * param.coh_ave;
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
  hdr.version = hdr_data(2,:);
  hdr.seconds = hdr_data(3,:);
  hdr.fraction = hdr_data(4,:);
  hdr.epri = hdr_data(5,:);
  hdr.num_wf = hdr_data(6,:);
  hdr.unknown1 = hdr_data(7,:);
  hdr.unknown2 = hdr_data(8,:);
  hdr.utc_time_sod = hdr.seconds + hdr.fraction/(param.clk/2);
  
  hdr.wfs = [];
  for wf = 1:hdr.num_wf(1)
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
  hdr.finfo = frame_sync_info(fn,struct('sync','DEADBEEF','cont_mode',0,'first_byte',param.first_byte));
  [fid,msg] = fopen(fn,'r','ieee-be');
  if param.recs(2) > hdr.finfo.num_rec-param.recs(1)
    param.recs(2) = hdr.finfo.num_rec-param.recs(1);
    recs = param.recs(1) + (0:param.recs(2)-1);
    num_rec = recs(end)-recs(1)+1;
  end
  % Round number of records to nearest multiple of the coherent averages
  param.recs(2) = floor(param.recs(2)/param.coh_ave) * param.coh_ave;

  % Read in data accounting for coherent averaging
  num_rlines = param.recs(2)/param.coh_ave;
  data = zeros(hdr.finfo.rec_size/2 - HEADER_SIZE/2 - 2, num_rlines, 'single');
  hdr_data = zeros(HEADER_SIZE/4, param.recs(2), 'double');
  rline = 1;
  for record = 1:param.recs(2)
    fseek(fid,hdr.finfo.syncs(param.recs(1)+record),-1);
    tmp_data = fread(fid,hdr.finfo.rec_size/2,'uint16=>float32');
    hdr_data(:,record) = 2^16*tmp_data(1:2:HEADER_SIZE/2) + tmp_data(2:2:HEADER_SIZE/2);
    data(:,rline) = data(:,rline) + tmp_data(HEADER_SIZE/2+3:end);
    if ~mod(record,param.coh_ave)
      rline = rline + 1;
    end
  end
  fclose(fid);
  data = data/param.coh_ave;
else
  slow_method = false;
end

hdr.frame_sync = hdr_data(1,:);
hdr.version = hdr_data(2,:);
hdr.seconds = hdr_data(3,:);
hdr.fraction = hdr_data(4,:);
hdr.epri = hdr_data(5,:);
hdr.num_wf = hdr_data(6,:);
hdr.unknown1 = hdr_data(7,:);
hdr.unknown2 = hdr_data(8,:);
hdr.utc_time_sod = hdr.seconds + hdr.fraction/(param.clk/2);

hdr.wfs = [];
for wf = 1:hdr.num_wf(1)
  hdr.wfs(wf).num_samp = hdr_data(9+2*(wf-1),1);
  %hdr.wfs(wf).unknown = hdr_data(10+2*(wf-1),1);
  hdr.wfs(wf).bit_shifts = mod(floor(hdr_data(10+2*(wf-1),1) / 2^24), 2^4);
  hdr.wfs(wf).t0 = mod(floor(hdr_data(10+2*(wf-1),1) / 2^10), 2^14) / param.clk - 8.3450e-06;
  hdr.wfs(wf).presums = 1 + mod(hdr_data(10+2*(wf-1),1), 2^10);
end

% Coherently average header data
hdr.frame_sync = fir_dec(hdr.frame_sync,param.coh_ave);
hdr.version = fir_dec(hdr.version,param.coh_ave);
hdr.seconds = fir_dec(hdr.seconds,param.coh_ave);
hdr.fraction = fir_dec(hdr.fraction,param.coh_ave);
hdr.epri = fir_dec(hdr.epri,param.coh_ave);
hdr.num_wf = fir_dec(hdr.num_wf,param.coh_ave);
hdr.unknown1 = fir_dec(hdr.unknown1,param.coh_ave);
hdr.unknown2 = fir_dec(hdr.unknown2,param.coh_ave);
hdr.utc_time_sod = fir_dec(hdr.utc_time_sod,param.coh_ave);

% ===================================================================
% Read and format data from raw file
% ===================================================================

% DEBUG CHECK:
%num_samp = (hdr.finfo.rec_size-41*4)/hdr.num_wf(1)/SAMP_SIZE

num_samp = hdr.wfs(1).num_samp;
num_samp_per_rec = hdr.num_wf(1)*num_samp;

if ~slow_method
  precision = sprintf('%d*uint16=>float32',num_samp_per_rec);
  fseek(fid, hdr.finfo.syncs(1) + hdr.finfo.rec_size*(recs(1)-1) + 81*2, 'bof');
  num_rlines = param.recs(2)/param.coh_ave;
  data = zeros(hdr.finfo.rec_size/2 - HEADER_SIZE/2 - 2, num_rlines, 'single');
  for rline = 1:num_rlines
    tmp_data = fread(fid,[num_samp_per_rec param.coh_ave],precision,hdr.finfo.rec_size-num_samp_per_rec*SAMP_SIZE);
    data(:,rline) = mean(tmp_data,2);
  end
  fclose(fid);
end

data = reshape(data,[num_samp hdr.num_wf(1) size(data,2)]);
% Trim off last 4 values (not valid data)
num_samp = num_samp - 4;
data = data(1:num_samp,:,:);
data = permute(data,[1 3 2]);

return;
