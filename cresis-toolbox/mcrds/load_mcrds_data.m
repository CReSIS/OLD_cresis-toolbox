function [param] = load_mcrds_data(param)
% [param] = load_mcrds_data(param)
%
% Function for loading MCRDS data. Supports loading over multiple files,
% pulse compression and fast decimation and presumming while loading
% to minimize memory usage. It supports loading arbitrary channels,
% waveform pairs and combining.  Multiple combination can be loaded
% simulataneously (each combination is referred to as an "image").
%
% Note: If pulse compression is not applied, then the param.radar.rx_path(rx).td
% correction is not applied and windowing is not applied.  These are both
% applied at the same time as the frequency domain matched filtering for
% pulse compression.
%
% Inputs:
% param: structure describing how the data is to be loaded
%  .load
%   .file_idx = file indexes for each record to load (from records)
%   .offset = offset into each record to load (from records)
%   .rec_data_size = size of all the data (not header) bytes in the record
%   .filenames = cell array of strings which file_idx references
%   .imgs = cell vector of wf_adc_list's. A wf_adc_list is an Nx2 array
%     of N (wf,adc) pairs from which the image will be formed. The first
%     column is the waveform and the second column is the receiver. These
%     wf/adc references are absolute (i.e. not relative indices)
%
%  .proc
%   .presums = number of presums to perform (may cause records on the
%     end to be dropped if not a factor of the number of records being loaded)
%   .combine_rx = true/false
%   .trim_vals = [2x1] integer vector, first entry is the number of
%     samples to drop at the beginning of the raw record and the second
%     entry is the number of samples to drop at the end of the raw record.
%   .pulse_rfi = struct for RFI removal
%    .en = true/false enable
%    .inc_ave = how many fast time samples to average
%    .tresh_scale = signal must exceed this threshold to remove
%    .pulse_comp = perform pulse compression
%    .ft_dec = fast time decimation
%
%  .radar = radar calibration information
%   .rx_path = struct array of receiver equalization coefficients
%    .chan_equal = scalar complex double (data DIVIDED by this)
%    .td = time delay correction (applied during pulse compression and
%       not applied if full compression or fast time)
%
%  .wfs = standard waveform processing structure (not all fields used)
%    This array is usually created with load_mcords_wfs.
%   .quantization_to_V = conversion factor from ADC quantization levels
%     to voltage
%   .Tpd = pulse duration
%   .f0 = start frequency of chirp
%   .f1 = stop frequency of chirp
%   .t0 = start time of recording relative to tx pulse
%   .blank = start and stop of receiver blank (leave empty if not blanked)
%   .Nt_ref = length of reference in fast-time
%   .Nt_raw = length of raw data in fast-time
%   .Nt_pc = length of pulse compressed data in fast-time
%   .offset = byte offset in the data part of the record
%   .ref = pulse compression reference waveform (conj. freq domain)
%   .time_raw = time axis of raw data
%   .freq_inds = frequency inds that will be kept after fast time decimation
%   .dc_shift = DC shift caused by fast time decimation
%   .fc = center frequency
%   .dt = sample spacing of output data product
%   .df = frequency spacing of output data product
%   .Nt = number of samples in output data product
%   .fs = sampling frequency of output data product
%   .freq = frequency axis of output data product
%   .time = time axis of output data product
%
% Outputs:
% param: updates to this structure are made to allow loading across
%   file boundaries
%
% Global Outputs:
% g_data: radar data is placed here
%
% Examples at bottom of file
%
% Authors: John Paden
%
% See also get_heights.m, csarp.m

if ~isfield(param.proc,'raw_data')
  param.proc.raw_data = false;
end

HEADER_SIZE = 160;

global g_data;
wfs = param.wfs;

physical_constants;

accum = build_img_load_struct(param.load.imgs, param.load.adcs);

% ===================================================================
% Preallocate data matrix
total_rec = param.load.recs(2)-param.load.recs(1)+1;
Nx = floor(total_rec/param.proc.presums);
if ~iscell(g_data)
  g_data = cell(size(param.load.imgs));
end
for img_idx = 1:length(param.load.imgs)
  wf = param.load.imgs{img_idx}(1,1);
  Nt = param.wfs(wf).Nt;
  if param.proc.combine_rx
    Nc = 1;
  else
    Nc = size(param.load.imgs{img_idx},1);
  end
  if any([size(g_data{img_idx},1) size(g_data{img_idx},2) size(g_data{img_idx},3)] ~= [Nt Nx Nc])
    g_data{img_idx} = zeros(Nt,Nx,Nc,'single');
  end
end

% ===================================================================
% Load data
% ===================================================================

num_accum = 0;
out_idx = 0;
rec = 1;
while rec < total_rec;
  % Get the filename and the number of records to load from this file
  % 1. We search for the file after the desired one
  fn_idx = find(param.load.file_rec_offset{1} > param.load.recs(1)+rec-1,1);
  num_rec = total_rec - rec + 1;
  if isempty(fn_idx)
    % 2. If no such file exists, then we want the last file
    fn_idx = length(param.load.filenames{1});
    fn = fullfile(param.load.filepath,param.load.filenames{1}{fn_idx});
  else
    % 3. If it does, grab the previous one and load the smaller of
    %    the remaining records to be loaded, block_size, or the number of records in
    %    this file. block_size is the max number of records to read at once.
    block_size = 500;
    fn_idx = fn_idx - 1;
    fn = fullfile(param.load.filepath,param.load.filenames{1}{fn_idx});
    num_rec = min([num_rec block_size param.load.file_rec_offset{1}(fn_idx+1) ...
      - (param.load.recs(1)+rec-1)]);
  end
  
  % ===============================================================
  % Load the records
  % ===============================================================
  
  fprintf('Loading data from %s\n', fn);
  fid = fopen(fn, 'r');
  % Seek to the offset of the current record to load
  %   Header bytes: 484 + 176*size(param.wfs.offsets,2)
  %   Records: 24 + 2*number of samples in record
  %   Current record header: 24 bytes
  fseek(fid,484 + 176*length(param.wfs) ...
    + (24+2*param.load.rec_data_size) * (param.load.recs(1)+rec-1-param.load.file_rec_offset{1}(fn_idx)) ...
    + 24,-1);
  % Load in records
  [rec_data num_items] = fread(fid,[param.load.rec_data_size num_rec], ...
    sprintf('%d*uint16=>single',param.load.rec_data_size),24);
  if num_items/param.load.rec_data_size ~= num_rec
    error('Read or calc error');
  end
  
  init_rec = rec;
  for rec = rec : rec+num_rec-1
    % ===============================================================
    % Presum in the accumulator
    % ===============================================================
    for accum_idx = 1:length(accum.wf)
      adc = accum.adc(accum_idx);
      wf = accum.wf(accum_idx);
      tmp = rec_data([1:wfs(wf).Nt_raw] + wfs(1).offset(adc,wf)-1,rec-init_rec+1);
      % Convert to volts
      tmp = (tmp-mean(tmp)) * wfs(wf).quantization_to_V;
      % Accumulate (presum)
      if num_accum == 0
        accum.data{accum_idx} = tmp;
      else
        accum.data{accum_idx} = accum.data{accum_idx} + tmp;
      end
    end
    
    % ===============================================================
    % Process record when presums are done
    % ===============================================================
    num_accum = num_accum + 1;
    if num_accum >= param.proc.presums
      num_accum = 0;
      out_idx = out_idx + 1;
      for accum_idx = 1:length(accum.wf)
        adc = accum.adc(accum_idx);
        wf = accum.wf(accum_idx);
        img_idx = accum.img_idx(accum_idx);
        wf_adc_idx = accum.wf_adc_idx(accum_idx);
        
        accum.data{accum_idx}([1:param.proc.trim_vals(1) end-param.proc.trim_vals(2)+1:end]) = 0;
        if param.proc.pulse_rfi.en
          pdata = abs(accum.data{accum_idx}).^2;
          inc_ave = param.proc.pulse_rfi.inc_ave;
          thresh_scale = param.proc.pulse_rfi.thresh_scale;
          thresh = filter(ones(inc_ave,1)/inc_ave,1,pdata);
          thresh = [thresh(1+(inc_ave-1)/2:end,:); repmat(thresh(end),[(inc_ave-1)/2 size(thresh,2)])];
          bad_idxs = find(pdata > thresh_scale*thresh);
          for bad_idx = bad_idxs.'
            if bad_idx >= 7 && bad_idx <= length(pdata)-6
              vals = pdata(bad_idx-6:bad_idx+6);
              thresh = median(vals) * 10^(17/10);
              bad_idxs2 = find(vals > thresh);
              accum.data{accum_idx}(bad_idx-6 + (bad_idxs2-1)) = 0;
            end
          end
        end
        % Apply channel compensation
        if ~param.proc.raw_data
          chan_equal = 10.^(param.radar.wfs(wf).chan_equal_dB(param.radar.wfs(wf).rx_paths(adc))/20) ...
            .* exp(1i*param.radar.wfs(wf).chan_equal_deg(param.radar.wfs(wf).rx_paths(adc)) ...
              + param.adc_phase_corr_deg(rec,adc))/180*pi);
          accum.data{accum_idx} = accum.data{accum_idx}/chan_equal;
          adc_gain = param.load.wfs(find(param.load.recs(1) + rec -1 >= param.load.wfs_records,1)).wfs(wf).adc_gains;
          accum.data{accum_idx} = accum.data{accum_idx}/adc_gain;
        end
        
        if param.proc.pulse_comp
          % ===========================================================
          % Do pulse compression
          % Apply blank (only should enable if sidelobe problems present)
          if ~isempty(wfs(wf).blank)
            % Blank is larger of two numbers passed in through radar worksheet blank parameter:
            %   Number 1 is added to surface time delay and is usually equal to pulse duration
            %   Number 2 is unmodified and is usually equal to hardware blank setting
            %   Set either number to -inf to disable
            blank_time = max(param.surface(rec) + wfs(wf).blank(1),wfs(wf).blank(2));
            %accum.data{accum_idx}(wfs(wf).time_raw>wfs(wf).blank(1) & wfs(wf).time_raw<wfs(wf).blank(2)) = 0;
            accum.data{accum_idx}(wfs(wf).time_raw-param.radar.wfs(wf).Tsys(adc) <= blank_time) = 0;
          end
          % Apply matched filter
          % Zero pad front: (the standard)
          accum.data{accum_idx} = fft([zeros(wfs(wf).pad_length,1); accum.data{accum_idx}]);
          % Zero pad end: (debug only)
          %accum.data{accum_idx} = fft(accum.data{accum_idx}, wfs(wf).Nt_pc);
          
          % Apply matched filter and transform back to time domain
          accum.data{accum_idx} = ifft(accum.data{accum_idx} .* wfs(wf).ref{adc});
          
          if param.proc.ft_dec
            % Digital down conversion and decimation
            accum.data{accum_idx} = accum.data{accum_idx}.*exp(-1i*2*pi*wfs(wf).fc*wfs(wf).time_raw);
            accum.data{accum_idx} = resample(double(accum.data{accum_idx}), param.wfs(1).ft_dec(1), param.wfs(1).ft_dec(2));
          end
          
        elseif param.proc.ft_dec
          accum.data{accum_idx} = fft(accum.data{accum_idx},wfs(wf).Nt_raw);
          accum.data{accum_idx} = ifft(accum.data{accum_idx}(wfs(wf).freq_inds));
          accum.data{accum_idx} = accum.data{accum_idx}.*exp(-1i*2*pi*wfs(wf).fc*wfs(wf).time_raw);
          accum.data{accum_idx} = resample(double(accum.data{accum_idx}), param.wfs(1).ft_dec(1), param.wfs(1).ft_dec(2));
        end
        if param.proc.combine_rx
          if wf_adc_idx == 1
            g_data{img_idx}(:,out_idx) = accum.data{accum_idx} / param.proc.presums / size(param.load.imgs{img_idx},1);
          else
            g_data{img_idx}(:,out_idx) = g_data{img_idx}(:,out_idx) + accum.data{accum_idx} / param.proc.presums / size(param.load.imgs{img_idx},1);
          end
        else
          g_data{img_idx}(:,out_idx,wf_adc_idx) = accum.data{accum_idx} / param.proc.presums;
        end
      end
    end
  end
  
  fclose(fid);
  rec = rec + 1;
  
  % ===============================================================
end

return;

% ===================================================================
% ===================================================================
% load_mcrds_data.m Examples
% ===================================================================
% ===================================================================

% None yet

% ===================================================================
% ===================================================================
% build_img_load_struct support function
% ===================================================================
% ===================================================================
function accum = build_img_load_struct(imgs, adcs)
% accum = build_img_load_struct(imgs, adcs)
%
% Builds an accumulator structure for load_mcrds_data. This is required
% to run load_mcrds_data.
%
% imgs = cell vector of imgs (each entry is a separate wf_adc_list)
%   Each wf_adc_list is an Nx2 array where N is the number of channels,
%   the first column is the waveform, and the second column is the receiver.
%   Absolute index of wf and adc are used in this array.
% adcs = vector of valid receivers (absolute index of each receiver)
%
% accum = structure that helps load_mcrds_data load data ("accumulator")
%   This structure contains lists of how to accumulate and store data 
%   when loading from load_mcrds_data.  load_mcrds_data first
%   accumulates the data (presums):
%     accum.data{accumulator-instance}
%   Once the presums are finished, the data is processed and stored in
%   the output variable:
%     g_data{img_idx}(fast-time,slow-time,wf_adc_idx)
%   The four fields in accum are all length Kx1 where K is the number
%   of accumulator instances.  The number of accumulator instances is
%   determined by the length of the fields.
%  .adc = K x 1 vector indicating which adc this instance is pulled from
%  .wf = K x 1 vector indicating which waveform this instance is pulled from
%  .wf_adc_idx = an index in the final output array unless receivers are
%    combined in which case size(g_data{:}, 3) == 1.
%  .img_idx = an index in the final output array
%
% Examples: At the bottom of this file
%
% Author: John Paden
%
% See also: load_mcrds_data

% Verify that all receiver entries are valid
for img_idx = 1:length(imgs)
  for wf_adc_idx = 1:size(imgs{img_idx},1)
    if ~ismember(imgs{img_idx}(wf_adc_idx,2), adcs)
      error('Rx %d is invalid', imgs{img_idx}(wf_adc_idx,2));
    end
  end
end

% Build accum structure
accum.adc = [];
accum.wf = [];
accum.wf_adc_idx = [];
accum.img_idx = [];
for img_idx = 1:length(imgs)
  for wf_adc_idx = 1:size(imgs{img_idx},1)
    accum.adc(end+1) = imgs{img_idx}(wf_adc_idx,2);
    accum.wf(end+1) = imgs{img_idx}(wf_adc_idx,1);
    accum.wf_adc_idx(end+1) = wf_adc_idx;
    accum.img_idx(end+1) = img_idx;
  end
end

return;

% ===================================================================
% ===================================================================
% build_img_load_struct Examples
% ===================================================================
% ===================================================================

% Example 1
imgs = {[1 6; 1 7; 1 8],[2 1; 2 2; 2 3; 2 4; 2 5]};
adcs = [1 2 3 4 5 7 8];
accum = build_img_load_struct(imgs, adcs);

% Example 2
imgs = {[1 6; 1 7; 1 8],[2 1; 2 2; 2 3; 2 4; 2 5]};
adcs = [1 2 3 4 5 6 7 8];
accum = build_img_load_struct(imgs, adcs);

% Example 3
imgs = {[1 1; 1 2; 1 3; 1 4; 1 5; 1 6; 1 7; 1 8], ...
  [2 1; 2 2; 2 3; 2 4; 2 5; 2 6; 2 7; 3 1; 3 2; 3 3; 3 4; 3 5; 3 6; 3 7; 3 8]};
adcs = [1 2 3 4 5 6 7 8];
accum = build_img_load_struct(imgs, adcs);



