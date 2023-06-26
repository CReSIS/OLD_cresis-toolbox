function [hdr,data] = basic_load_mcrds(filename,param)
% [hdr,data] = basic_load_mcrds(filename,param)
%
% filename: string containing the path to a raw MCRDS file
% param = struct controlling loading of data
%   .recs = 2 element vector for records to load [start_rec num_rec]
%     start_rec uses zero-indexing (default is [0 inf] or all records)
%     num_rec may be infinity, which means to end of file
%   .coh_ave = double scalar, number of coherent averages, (default
%     is 1)
%
% hdr = hdr fields from file
%   .comp_time: 1 by Nx double vector of computer time for each record
%   .radar_time: 1 by Nx double vector of radar time for each record
% data = cell vector of waveform data
%   data{wf}(Nt,Nx,Nc) where wf = waveform, Nt = fast-time, Nx = slow-time
%   Nc = adc channels
%
% Example
%   fn = '/cresis/data1/MCRDS/2008_Greenland/20080627A/dataGISMO/data.20080627100855.0001.raw';
%   [hdr,data] = basic_load_mcrds(fn);
%
% Author: Anthony Hoch, John Paden
%
% See also: basic_load_mcrds_hdr.m, basic_load_mcrds_time.m,
%   basic_load_mcords.m, basic_load_mcords2.m, basic_load_fmcw.m,
%   basic_load_accum.m

% ===================================================================
% Check input arguments
% ===================================================================
if ~exist('param','var') || isempty(param)
  param.recs = [0 inf];
  param.coh_ave = 1;
end
if ~isfield(param,'recs');
  param.recs = [0 inf];
end
if ~isfield(param,'coh_ave');
  param.coh_ave = 1;
end
coh_ave = param.coh_ave;

[fid,msg] = fopen(filename,'r');
if fid == -1
  fprintf('Could not open file (%s)',filename);
  warning(msg);
  return
end

% Determine the file size
fseek(fid,0,1);
file_size = ftell(fid);
fseek(fid,0,-1);

% Initialization for reading in headers
hdrs = [];
curNumDataRecords = 0;

numSam(1) = 0;
numSam(2) = 0;

% ------------------------------------------------------------
% ------------------------------------------------------------
% Read in Generic File Header
% ------------------------------------------------------------
% ------------------------------------------------------------
[fileType itemsreturned] = fread(fid,32,'char');
[year itemsreturned] = fread(fid,1,'int16');
[subYearVer itemsreturned] = fread(fid,1,'int8');
status = fseek(fid,5,'cof');

[datatype itemsreturned] = fread(fid,1,'int32');
if itemsreturned < 1
  fprintf('%s\n',ferror(fid));
  fclose(fid);
  error(sprintf('File format error (truncated hdr)'));
end

% ------------------------------------------------------------
% ------------------------------------------------------------
% Read the hdr
% ------------------------------------------------------------
% ------------------------------------------------------------
[hdr.sampFreq itemsreturned] ...
  = fread(fid,1,'float64');
if itemsreturned < 1
  fclose(fid);
  error(sprintf('File format error (truncated hdr)'));
end
[hdr.prfCount itemsreturned] = fread(fid,1,'int32');
if itemsreturned < 1
  fclose(fid);
  error(sprintf('File format error (truncated hdr)'));
end
[hdr.numAve itemsreturned] = fread(fid,1,'int32');
if itemsreturned < 1
  fclose(fid);
  error(sprintf('File format error (truncated hdr)'));
end
[hdr.rxAtten itemsreturned] = fread(fid,64,'int32');
if itemsreturned < 64
  fclose(fid);
  error(sprintf('File format error (truncated hdr)'));
end
hdr.rxAtten ...
  = permute(reshape(hdr.rxAtten,[4 2 8]),[3 2 1]);
[hdr.rxBlank itemsreturned] = fread(fid,32,'int32');
if itemsreturned < 32
  fclose(fid);
  error(sprintf('File format error (truncated hdr)'));
end
hdr.rxBlank ...
  = permute(reshape(hdr.rxBlank,[4 1 8]),[3 2 1]);
[hdr.calModeEn itemsreturned] = fread(fid,1,'int8');
if itemsreturned < 1
  fclose(fid);
  error(sprintf('File format error (truncated hdr)'));
end
[hdr.numWaveforms itemsreturned] ...
  = fread(fid,1,'int8');
if itemsreturned < 1
  fclose(fid);
  error(sprintf('File format error (truncated hdr)'));
end
status = fseek(fid,2,'cof');
if status == -1
  fprintf('%s\n',ferror(fid));
  fclose(fid);
  error(sprintf('File format error (truncated hdr)'));
end
[hdr.calNumOfPnts itemsreturned] ...
  = fread(fid,1,'int32');
if itemsreturned < 1
  fclose(fid);
  error(sprintf('File format error (truncated hdr)'));
end
[hdr.calStartFreq itemsreturned] ...
  = fread(fid,1,'float64');
if itemsreturned < 1
  fclose(fid);
  error(sprintf('File format error (truncated hdr)'));
end
[hdr.calStopFreq itemsreturned] ...
  = fread(fid,1,'float64');
if itemsreturned < 1
  fclose(fid);
  error(sprintf('File format error (truncated hdr)'));
end
[hdr.calDelay itemsreturned] ...
  = fread(fid,1,'float64');
if itemsreturned < 1
  fclose(fid);
  error(sprintf('File format error (truncated hdr)'));
end
[hdr.calDuration itemsreturned] ...
  = fread(fid,1,'float64');
if itemsreturned < 1
  fclose(fid);
  error(sprintf('File format error (truncated hdr)'));
end
% ------------------------------------------------------------
% -- Read the waveforms
% 1. Keep track of the number of samples per daq channel
numSam = zeros(8,1);
for ind = 1:hdr.numWaveforms
  [hdr.wf(ind).startFreq itemsreturned] ...
    = fread(fid,1,'float64');
  if itemsreturned < 1
    fclose(fid);
    error(sprintf('File format error (truncated hdr)'));
  end
  [hdr.wf(ind).stopFreq itemsreturned] ...
    = fread(fid,1,'float64');
  if itemsreturned < 1
    fclose(fid);
    error(sprintf('File format error (truncated hdr)'));
  end
  [hdr.wf(ind).pulseDuration itemsreturned] ...
    = fread(fid,1,'float64');
  if itemsreturned < 1
    fclose(fid);
    error(sprintf('File format error (truncated hdr)'));
  end
  [hdr.wf(ind).calFreq itemsreturned] ...
    = fread(fid,1,'float64');
  if itemsreturned < 1
    fclose(fid);
    error(sprintf('File format error (truncated hdr)'));
  end
  [hdr.wf(ind).calDelay itemsreturned] ...
    = fread(fid,1,'float64');
  if itemsreturned < 1
    fclose(fid);
    error(sprintf('File format error (truncated hdr)'));
  end
  [hdr.wf(ind).calDuration itemsreturned] ...
    = fread(fid,1,'float64');
  if itemsreturned < 1
    fclose(fid);
    error(sprintf('File format error (truncated hdr)'));
  end
  [hdr.wf(ind).bandSelect itemsreturned] ...
    = fread(fid,1,'int8');
  if itemsreturned < 1
    fclose(fid);
    error(sprintf('File format error (truncated hdr)'));
  end
  [hdr.wf(ind).zeroPiMod itemsreturned] ...
    = fread(fid,1,'int8');
  if itemsreturned < 1
    fclose(fid);
    error(sprintf('File format error (truncated hdr)'));
  end
  [hdr.wf(ind).txMult itemsreturned] ...
    = fread(fid,1,'int8');
  if itemsreturned < 1
    fclose(fid);
    error(sprintf('File format error (truncated hdr)'));
  end
  [hdr.wf(ind).rxMode itemsreturned] ...
    = fread(fid,1,'int8');
  if itemsreturned < 1
    fclose(fid);
    error(sprintf('File format error (truncated hdr)'));
  end
  [hdr.wf(ind).txAmpEn itemsreturned] ...
    = fread(fid,2,'int8');
  if itemsreturned < 2
    fclose(fid);
    error(sprintf('File format error (truncated hdr)'));
  end
  status = fseek(fid,2,'cof');
  if status == -1
    fprintf('%s\n',ferror(fid));
    fclose(fid);
    error(sprintf('File format error (truncated hdr)'));
  end
  [hdr.wf(ind).modCount0 itemsreturned] ...
    = fread(fid,1,'int32');
  if itemsreturned < 1
    fclose(fid);
    error(sprintf('File format error (truncated hdr)'));
  end
  [hdr.wf(ind).modCount1 itemsreturned] ...
    = fread(fid,1,'int32');
  if itemsreturned < 1
    fclose(fid);
    error(sprintf('File format error (truncated hdr)'));
  end
  [hdr.wf(ind).numSam itemsreturned] ...
    = fread(fid,4,'int32');
  if itemsreturned < 4
    fclose(fid);
    error(sprintf('File format error (truncated hdr)'));
  end
  if ind == 1
    hdr.wf(ind).numSam = hdr.wf(ind).numSam - 1;
  end
  hdr.wf(ind).numSam ...
    = reshape(repmat(hdr.wf(ind).numSam,[1 2]).',[2*length(hdr.wf(ind).numSam) 1]);
  [hdr.wf(ind).sampleDelayCount itemsreturned] ...
    = fread(fid,4,'int32');
  if itemsreturned < 4
    fclose(fid);
    error(sprintf('File format error (truncated hdr)'));
  end
  hdr.wf(ind).sampleDelayCount ...
    = reshape(repmat(hdr.wf(ind).sampleDelayCount,[1 2]).',[2*length(hdr.wf(ind).sampleDelayCount) 1]);
  [hdr.wf(ind).recordEn itemsreturned] ...
    = fread(fid,8,'int8');
  if itemsreturned < 8
    fclose(fid);
    error(sprintf('File format error (truncated hdr)'));
  end
  [hdr.wf(ind).recordStart itemsreturned] ...
    = fread(fid,8,'int32');
  if itemsreturned < 8
    fclose(fid);
    error(sprintf('File format error (truncated hdr)'));
  end
  [hdr.wf(ind).recordStop itemsreturned] ...
    = fread(fid,8,'int32');
  if itemsreturned < 8
    fclose(fid);
    error(sprintf('File format error (truncated hdr)'));
  end
  if hdr.wf(ind).recordEn(1)
    numSam(1) = numSam(1) + hdr.wf(ind).numSam(1);
  end
  if hdr.wf(ind).recordEn(2)
    numSam(2) = numSam(2) + hdr.wf(ind).numSam(1);
  end
  if hdr.wf(ind).recordEn(3)
    numSam(3) = numSam(3) + hdr.wf(ind).numSam(2);
  end
  if hdr.wf(ind).recordEn(4)
    numSam(4) = numSam(4) + hdr.wf(ind).numSam(2);
  end
  if hdr.wf(ind).recordEn(5)
    numSam(5) = numSam(5) + hdr.wf(ind).numSam(3);
  end
  if hdr.wf(ind).recordEn(6)
    numSam(6) = numSam(6) + hdr.wf(ind).numSam(3);
  end
  if hdr.wf(ind).recordEn(7)
    numSam(7) = numSam(7) + hdr.wf(ind).numSam(4);
  end
  if hdr.wf(ind).recordEn(8)
    numSam(8) = numSam(8) + hdr.wf(ind).numSam(4);
  end
  [hdr.wf(ind).blankDelayCount itemsreturned] ...
    = fread(fid,2,'int32');
  if itemsreturned < 2
    fclose(fid);
    error(sprintf('File format error (truncated hdr)'));
  end
end
hdr.numSam = numSam;

% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% Prepare to read in Data
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

% Determine offsets for waveform wfInd
%   Create index offsets into data block such that:
%     offsets(daq,wfInd) points to the first sample
%     atten(daq,wfInd) gives the attenuation settings
%   hdr.wf(wfInd).numSam(daq) already contains the number of samples
offInd = 1;
offset = 1;
for daq = 1:8
  for wfInd = 1:length(hdr.wf)
    offsets(daq,wfInd) = offset;
    atten(daq,wfInd) = hdr.rxAtten(daq,1,hdr.wf(wfInd).rxMode+1) ...
      + hdr.rxAtten(daq,2,hdr.wf(wfInd).rxMode+1);
    offset = offset + hdr.wf(wfInd).numSam(daq);
  end
end


% Determine the number of complete records in the file
rec_size = (4+8+8+4+2*sum(hdr.numSam));
num_rec = floor((file_size - ftell(fid)) / rec_size);
if param.recs(1) >= num_rec
  error('Start record is beyond end of file');
end
if param.recs(1) + param.recs(2) > num_rec
  param.recs(2) = num_rec - param.recs(1);
end

% Allocate data
hdr.comp_time = zeros(1,floor(param.recs(2)/coh_ave));
hdr.radar_time = zeros(1,floor(param.recs(2)/coh_ave));
for wfInd = 1:length(hdr.wf)
  % 1. Determine the number of DAQs recorded for this waveform
  numDaq = sum(hdr.wf(wfInd).recordEn);
  % 2. Determine the longest record (account for error sample in first waveform)
  numBins = max(hdr.wf(wfInd).numSam);
  % 3. Allocate memory
  data.wf{wfInd} = zeros(numBins,floor(param.recs(2)/param.coh_ave),numDaq);
end

% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% Read in data
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
offset_index = 0;
index = 1;
indexRaw = 1;
cohInd = 1;
% Skip to the first record to load, param.recs(1)
fseek(fid,param.recs(1)*rec_size,'cof');

% ----------------------------------------------------------------------
% Begin reading in data
while indexRaw <= param.recs(2)
  % Read datatype
  [dt itemsreturned] = fread(fid,1,'int32');
  if itemsreturned < 1
    warning(sprintf('Only %.0f records instead of requested %.0f records',indexRaw-1,param.recs(2)));
    if (cohInd == 1)
      for wfInd = 1:length(hdr.wf)
        data.wf{wfInd} = data.wf{wfInd}(:,1:index-1,:);
      end
      comp_time = comp_time(1:index-1);
      radar_time = radar_time(1:index-1);
    else
      for wfInd = 1:length(hdr.wf)
        data.wf{wfInd} = data.wf{wfInd}(:,1:index,:);
      end
      comp_time = comp_time(1:index);
      radar_time = radar_time(1:index);
    end
    break;
  end
  % Read computer timestamp
  [seconds itemsreturned] = fread(fid,1,'int32');
  if itemsreturned < 1
    fclose(fid);
    error(sprintf('File format error (truncated timestamp)'));
  end
  [useconds itemsreturned] = fread(fid,1,'int32');
  if itemsreturned < 1
    fclose(fid);
    error(sprintf('File format error (truncated timestamp)'));
  end
  if cohInd == 1
    comp_time(index) = seconds + (1e-6 * useconds);
  elseif cohInd == param.coh_ave
    comp_time(index) = 1/param.coh_ave * (seconds + (1e-6 * useconds) + comp_time(index));
  else
    comp_time(index) = seconds + (1e-6 * useconds) + comp_time(index);
  end
  % Read radar timestamp
  [curRadarTime itemsreturned] = fread(fid,1,'int64');
  if itemsreturned < 1
    fclose(fid);
    error(sprintf('File format error (truncated timestamp)'));
  end
  if index > 1 & abs(curRadarTime - oldRadarTime ...
      - 10e6*hdr.prfCount/10e6*length(hdr.wf)*hdr.numAve) > 1
    fprintf('Skip in radar time: %.0f, should be %.0f\n',curRadarTime - oldRadarTime,...
      10e6*hdr.prfCount/10e6*length(hdr.wf)*hdr.numAve);
  end
  oldRadarTime = curRadarTime;
  if cohInd == 1
    radar_time(index) = curRadarTime;
  elseif cohInd == param.coh_ave
    radar_time(index) = 1/param.coh_ave * (curRadarTime + radar_time(index));
  else
    radar_time(index) = curRadarTime + radar_time(index);
  end

  % Skip past DAQ error bytes
  status = fseek(fid,4,'cof');
  
  % Read raw data
  [tmpData itemsreturned] = ...
    fread(fid,sum(hdr.numSam),'uint16');
  if itemsreturned < sum(hdr.numSam)-length(find(hdr.numSam~=0))
    fprintf('%s\n',ferror(fid));
    fclose(fid);
    error(sprintf('File format error (truncated data)'));
  end
  % Parse/sort data (also coherently average)
  if cohInd == 1
    for daq = find(hdr.numSam~=0).'
      for wfInd = 1:length(hdr.wf)
        data.wf{wfInd}(1:hdr.wf(wfInd).numSam(daq),index,daq) ...
          = tmpData(offsets(daq,wfInd):offsets(daq,wfInd)+hdr.wf(wfInd).numSam(daq)-1);
      end
    end
    if (cohInd ~= param.coh_ave)
      cohInd = cohInd + 1;
    else
      index = index  + 1;
    end
  elseif cohInd == param.coh_ave
    for daq = find(hdr.numSam~=0).'
      for wfInd = 1:length(hdr.wf)
        data.wf{wfInd}(1:hdr.wf(wfInd).numSam(daq),index,daq) ...
          = 1/param.coh_ave*(data.wf{wfInd}(1:hdr.wf(wfInd).numSam(daq),index,daq) ...
          + tmpData(offsets(daq,wfInd):offsets(daq,wfInd)+hdr.wf(wfInd).numSam(daq)-1));
      end
    end
    index = index + 1;
    cohInd = 1;
  else
    for daq = find(hdr.numSam~=0).'
      for wfInd = 1:length(hdr.wf)
        data.wf{wfInd}(1:hdr.wf(wfInd).numSam(daq),index,daq) ...
          = data.wf{wfInd}(1:hdr.wf(wfInd).numSam(daq),index,daq) ...
          + tmpData(offsets(daq,wfInd):offsets(daq,wfInd)+hdr.wf(wfInd).numSam(daq)-1);
      end
    end
    cohInd = cohInd + 1;
  end

  indexRaw = indexRaw + 1;
end

fclose(fid);


if length(comp_time) > 1
  prfPeriod = mean(diff(comp_time));
  p = polyfit([0:1:length(comp_time)-1],comp_time,1);
  comp_time = polyval(p,[0:1:length(comp_time)-1]) - 0.5*prfPeriod ...
    + 0.5*prfPeriod*hdr.prfCount/10e6*8;
end

% Reformat radar time to seconds
radar_time = radar_time/10.0e6;

hdr.radar_time = radar_time;
hdr.comp_time = comp_time;

return;
