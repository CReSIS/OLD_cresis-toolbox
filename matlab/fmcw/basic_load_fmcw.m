function [hdr,data] = basic_load_fmcw(fn, param)
% [hdr,data] = basic_load_fmcw(fn, param)
%
% This is the only function which loads raw data directly.
%
% Loads a single FMCW (snow,kuband) radar file. This is primarily for debugging.
% NOTE: 64-bit computer is required to load a 512 MB file since it will
% consume 1 GB of memory after loading.
%
% fn = filename of FMCW data
% param = struct controlling loading of data
%   .clk = clock (Hz), default one, used to interpret
%     counts in the header fields
%     2010 Antarctica FMCW radars use sampling frequency 1e9/16
%   .utc_time_halved = boolean (default 0)
%   .recs = 2 element vector for records to load [start_rec num_rec]
%     start_rec uses zero-indexing. An error is thrown if records beyond
%     the end of file are requested.
%   .first_byte = First byte to use. Default is zero. Primarily
%      used for skipping over "garbage" data at the front of the file
%      (scalar integer)
%   .records
%     .en = NOT USED. If true, all
%       headers will be loaded and the output arguments become:
%       hdr --> success flag (returns 1)
%       data --> hdr (every header)
%     .epri = Expected EPRI
%     .force_all = boolean to force loading all record headers
%
% hdr = file header for each record (if data is not defined, behavior
%   depends on param.records.en variable, if it is false only the first hdr
%   is returned, otherwise all the headers are returned)
% data = (optional output) array of radar data where dimensions are
%   1: fast-time/range-bin
%   2: slow-time/range-line/records
%
% Examples: see bottom of file
%
% Authors: John Paden
%
% See also basic_load_*.m

% ===================================================================
% Check input arguments
% ===================================================================
if ~exist('param','var') || isempty(param)
  param.clk = 1;
  param.utc_time_halved = 0;
  param.recs = [];
  param.first_byte = 0;
end
if ~isfield(param,'clk');
  param.clk = 1;
end
  % We need to make sure the field is 0 or 1 for computations done below
if ~isfield(param,'utc_time_halved') || isempty(param.utc_time_halved) ...
    || ~param.utc_time_halved
  param.utc_time_halved = 0;
elseif param.utc_time_halved
  param.utc_time_halved = 1;
end
if ~isfield(param,'recs');
  param.recs = [];
end
if ~isfield(param,'first_byte');
  param.first_byte = 0;
end
if ~isfield(param,'records');
  param.records.en = false;
end

% ===================================================================
% Get the locations of the first sync markers
% ===================================================================
hdr.finfo.syncs = get_first10_sync_mfile(fn, param.first_byte);
% Find the record size in a robust way (rec_size in bytes)
hdr.finfo.rec_size = median(diff(hdr.finfo.syncs));
% Find the first good header
hdr.finfo.syncs = hdr.finfo.syncs(find(diff(hdr.finfo.syncs) == hdr.finfo.rec_size));

% ===================================================================
% Determine the size of the output matrix
% ===================================================================
[fid,msg] = fopen(fn,'r','ieee-be');
if fid < 1
  fprintf('Could not open file %s\n', fn);
  error(msg);
end

fseek(fid,0,'eof');
hdr.finfo.file_size = ftell(fid);

hdr.finfo.num_rec = floor((hdr.finfo.file_size-hdr.finfo.syncs(1))/hdr.finfo.rec_size);

% Number of slow time records
param_recs_set_by_function = false;
if isempty(param.recs)
  param.recs = [0 hdr.finfo.num_rec];
  param_recs_set_by_function = true;
else
  if param.recs(1) + param.recs(2) > hdr.finfo.num_rec
    error('Only %i records in file',hdr.finfo.num_rec);
  end
end

% Seek to first record
sync_offset = 1;
fseek(fid, hdr.finfo.syncs(sync_offset) + hdr.finfo.rec_size*param.recs(1), 'bof');

% ===================================================================
% Read in file
%   See MCORDS file format.docx
% ===================================================================
HEADER_SIZE = 7*4; % in bytes

if nargout < 2 && ~param.records.en
  % Only header needs to be returned in output argument, so we read
  % only the first record
  hdr_data = fread(fid,[HEADER_SIZE/4 1],sprintf('%i*uint32',HEADER_SIZE/4),hdr.finfo.rec_size-HEADER_SIZE);
  fclose(fid);
  hdr_data(5:6) = hdr_data(5:6) / (1+param.utc_time_halved);
  
  hdr.frame_sync = hdr_data(1,:);
  hdr.epri = hdr_data(2,:);
  hdr.unknown2 = hdr_data(3,:);
  hdr.unknown3 = hdr_data(4,:);
  hdr.seconds = hdr_data(5,:);
  hdr.fraction = hdr_data(6,:);
  hdr.unknown4 = hdr_data(7,:);
  
  hdr.wfs = [];
  wf = 1;
  hdr.wfs(wf).num_sam = (hdr.finfo.rec_size)/2 - 14;
  hdr.wfs(wf).bit_shifts = 0;
  hdr.wfs(wf).Tadc = 0;
  hdr.wfs(wf).presums = 4;
  
  hdr.utc_time_sod = hdr.seconds + hdr.fraction/param.clk;
  
  return;
  
elseif param.records.en
  loading_failed = false;
  
  if ~param.records.force_all
    %% Quick Read (only looks at first and last record)
    
    % Read in first record
    fseek(fid, param.recs(1) * hdr.finfo.rec_size + hdr.finfo.syncs(sync_offset), -1);
    [hdr_data,samples_read] = fread(fid, HEADER_SIZE/4, 'uint32=>uint32');
    hdr_data(5:6) = hdr_data(5:6) / (1+param.utc_time_halved);
    
    % Read in last record (this assumes an ideal file)
    fseek(fid, (param.recs(2)-1) * hdr.finfo.rec_size + hdr.finfo.syncs(sync_offset), -1);
    [last_hdr,samples_read] = fread(fid, HEADER_SIZE/4, 'uint32=>uint32');
    last_hdr(5:6) = last_hdr(5:6) / (1+param.utc_time_halved);
    
    % Check to see if headers are both correct in terms of frame sync, EPRI
    % number and settings have not changed.
    if ~(last_hdr(1) == hdr_data(1) ...
        && last_hdr(2) == hdr_data(2) + param.recs(2)-1)
      loading_failed = true;
    end
  end
  
  if ~loading_failed && ~param.records.force_all
    %% Quick Load using just first and last record
    fclose(fid);
    hdr.epri = hdr_data(2):last_hdr(2);
    hdr.seconds = hdr_data(5) .* ones(size(hdr.epri),'uint32');
    hdr.fraction = double(hdr_data(6)) + (0:param.recs(2)-1) * param.records.epri*param.clk;
    need_second_jump = find(hdr.fraction/param.clk > 1,1);
    while ~isempty(need_second_jump)
      hdr.fraction(need_second_jump:end) = hdr.fraction(need_second_jump:end) - param.clk;
      hdr.seconds(need_second_jump:end) = hdr.seconds(need_second_jump:end) + 1;
      need_second_jump = find(hdr.fraction/param.clk > 1,1);
    end
    hdr.fraction = uint32(hdr.fraction);
    if hdr.seconds(end) ~= last_hdr(5)
      warning('Potential error because hdr.seconds(end) (%d) and last_hdr(%d) should be the same', hdr.seconds(end), last_hdr(5));
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
      hdr_data(5:6,:) = hdr_data(5:6,:) / (1+param.utc_time_halved);
      hdr.frame_sync = hdr_data(1,:);
    end
    
    if loading_failed || any(hdr.frame_sync ~= hdr.frame_sync(1))
      %% Slow Header Load Case
      fprintf('  Loss of frame sync, loading file the slow way\n');
      hdr.finfo = frame_sync_info(fn,struct('sync','DEADBEEF','cont_mode',0));
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
        hdr_data(:,record) = fread(fid,HEADER_SIZE/4,'uint32=>uint32');
      end
      fclose(fid);
      
      % Parse header data
      hdr_data = double(hdr_data);
      hdr_data(5:6,:) = hdr_data(5:6,:) / (1+param.utc_time_halved);
      
      hdr.offset = uint32(hdr.finfo.syncs(param.recs(1)+1:param.recs(2)));
    else
      hdr.offset = uint32(hdr.finfo.syncs(sync_offset) + hdr.finfo.rec_size*(0:param.recs(2)-1));
      hdr = rmfield(hdr,'frame_sync');
    end
    
    hdr.epri = hdr_data(2,:);
    hdr.seconds = hdr_data(5,:);
    hdr.fraction = hdr_data(6,:);
  end
  
  hdr.wfs = [];
  wf = 1;
  hdr.wfs(wf).presums = NaN;
  hdr.wfs(wf).bit_shifts = 0;
  hdr.wfs(wf).Tadc = NaN;
  hdr.wfs(wf).start_idx = NaN;
  hdr.wfs(wf).num_sam = (hdr.finfo.rec_size)/2 - 14;
  hdr = rmfield(hdr,'finfo');
  
  %% Remap outputs to match create_task.m output standard:
  %   [hdr,data] --> [success,hdr]
  data = hdr;
  hdr = 1;
  return;
end

% Load in data matrix all at once
data = fread(fid,[hdr.finfo.rec_size/2 param.recs(2)],'uint16=>single');

bad_syncs = find(data(1,:) ~= hex2dec('DEAD') | data(2,:) ~= hex2dec('BEEF'));
if ~isempty(bad_syncs)
  fprintf('  Bad sync marks, rereading data file the slow way\n');
  fprintf('  Searching for individual sync marks\n');
  hdr.finfo = frame_sync_info(fn,struct('cont_mode',false));
  fprintf('  Reading in data file\n');
  
  % Number of slow time records
  if param.recs(2) > hdr.finfo.num_rec-param.recs(1)
    param.recs(2) = hdr.finfo.num_rec-param.recs(1);
  else
    if param_recs_set_by_function
      param.recs(2) = hdr.finfo.num_rec-param.recs(1);
    end
  end
  
  % Load in data
  fid = fopen(fn,'r','ieee-be');
  data = zeros(hdr.finfo.rec_size/2,param.recs(2),'single');
  for rec = 1:param.recs(2)
    fseek(fid,hdr.finfo.syncs(param.recs(1) + rec),-1);
    data(:,rec) = fread(fid,hdr.finfo.rec_size/2,'uint16');
  end
  fclose(fid);
  fprintf('  Done reading data file\n');
else
  fclose(fid);
end

hdr.frame_sync = 2^16*double(data(1,:))+double(data(2,:));
hdr.epri = 2^16*double(data(3,:))+double(data(4,:));
hdr.unknown2 = 2^16*double(data(5,:))+double(data(6,:));
hdr.unknown3 = 2^16*double(data(7,:))+double(data(8,:));
hdr.seconds = (2^16*double(data(9,:))+double(data(10,:))) / (1+param.utc_time_halved);
hdr.fraction = (2^16*double(data(11,:))+double(data(12,:))) / (1+param.utc_time_halved);
hdr.unknown4 = 2^16*double(data(13,:))+double(data(14,:));

hdr.wfs = [];
wf = 1;
hdr.wfs(wf).num_sam = size(data,1)-14;
hdr.wfs(wf).bit_shifts = 0;
hdr.wfs(wf).Tadc = 0;
hdr.wfs(wf).presums = 4;

% Each pair of samples are flipped around, so we fix that here
data = data(reshape([16:2:end;15:2:end-1],[size(data,1)-14 1]),:);

hdr.utc_time_sod = double(hdr.seconds) + double(hdr.fraction)/param.clk;

return;


% ===================================================================
% Examples
% ===================================================================

% Loading the whole file
fn = 'E:\snow\data00.10302010.0092.dat';
[hdr,data] = basic_load_fmcw(fn,struct('clk',1e9/16));
data = data(1:end-1000,:);
for rline = 1:size(data,2)
  data(:,rline) = data(:,rline) - mean(data(:,rline));
end
for rbin = 1:size(data,1)
  data(rbin,:) = data(rbin,:) - mean(data(rbin,:));
end
data = fft(data);
imagesc(lp(data));

% Loading records range 7000 to 8499
fn = 'E:\snow\data00.10302010.0092.dat';
[hdr,data] = basic_load_fmcw(fn,struct('clk',1e9/16,'recs',[7000 1500]));
data = data(1:end-1000,:);
for rline = 1:size(data,2)
  data(:,rline) = data(:,rline) - mean(data(:,rline));
end
for rbin = 1:size(data,1)
  data(rbin,:) = data(rbin,:) - mean(data(rbin,:));
end
data = fft(data);
imagesc(lp(data));
