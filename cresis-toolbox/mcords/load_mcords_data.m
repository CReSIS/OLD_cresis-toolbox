function [param] = load_mcords_data(param)
% [param] = load_mcords_data(param)
%
% Function for loading MCoRDS data. Supports loading over multiple files,
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
%      see .adcs for intrepretation
%   .offset = offset into each record to load (from records)
%      see .adcs for intrepretation
%   .rec_data_size = size of all the data (not header) bytes in the record
%   .filenames = cell array of strings which file_idx references
%   .adcs = list of unique adcs that are being loaded (offset and file_idx
%      are cell vectors of the same length as this vector and this vector
%      tells which entry in offset/file_idx goes to which adc)
%   .imgs = cell vector of wf_adc_list's. A wf_adc_list is an Nx2 array
%     of N (wf,adc) pairs from which the image will be formed. The first
%     column is the waveform and the second column is the adc. These
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
%    This array is usually created with load_mcords_wfs and is explained
%    in that function.
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
% See also get_heights.m, csarp.m, load_mcords.m, load_mcords_wfs.m

if ~isfield(param.load,'wf_adc_comb')
  param.load.wf_adc_comb.en = 0;
end

HEADER_SIZE = 160;

global g_data;
wfs = param.wfs;

physical_constants;

accum = build_img_load_struct(param.load.imgs, param.load.adcs);

% ===================================================================
% Preallocate data matrix
total_rec = length(param.load.file_idx{1});
Nx = floor(total_rec/param.proc.presums);
if ~iscell(g_data)
  g_data = cell(size(param.load.imgs));
end
for img_idx = 1:length(param.load.imgs)
  wf = param.load.imgs{img_idx}(1,1);
  if param.load.wf_adc_comb.en
    Nt = param.load.wf_adc_comb.Nt;
  else
    Nt = param.wfs(wf).Nt;
  end
  if param.proc.combine_rx
    Nc = 1;
  else
    Nc = size(param.load.imgs{img_idx},1);
  end
  if any([size(g_data{img_idx},1) size(g_data{img_idx},2) size(g_data{img_idx},3)] ~= [Nt Nx Nc])
    g_data{img_idx} = zeros(Nt,Nx,Nc,'single');
  elseif param.proc.combine_rx
    g_data{img_idx}(:) = 0;
  end
end

% ===================================================================
% Load data
% ===================================================================
for adc_idx = 1:length(param.load.adcs)
  adc = param.load.adcs(adc_idx);

  if param.load.offset{adc_idx}(1) < 0 && param.load.offset{adc_idx}(1) ~= -2^31
    % Get the last incomplete record from the previous file
    fn_idx = param.load.file_idx{adc_idx}(1) - 1;
    fn = param.load.filenames{adc_idx}{fn_idx};
    [fid,msg] = fopen(fn,'r');
    if fid<0
      error('File open failed (%s)\n%s',fn, msg);
    end
    fseek(fid,param.load.offset{adc_idx}(1),'eof');
    last_record = fread(fid,-double(param.load.offset{adc_idx}(1)),'uint8');
    fclose(fid);
    get_last_record = true;
  else
    % Do not need the last incomplete record from the previous file
    last_record = uint8([]);
    get_last_record = false;
  end

  num_accum = 0;
  out_idx = 0;
  rec = 1;
  while rec < total_rec;
    % Get the filename
    fn_idx = param.load.file_idx{adc_idx}(rec);
    fn = param.load.filenames{adc_idx}{fn_idx};

    % Get number of records to load
    num_rec = find(param.load.file_idx{adc_idx}(rec:total_rec ) > fn_idx, 1) - 1;

    % Check to see if we crossed a file boundary
    if isempty(num_rec)
      % Not crossing a file boundary
      num_rec = length(rec:total_rec);
      get_last_record = false;
    else
      % Crossing a file boundary
      get_last_record = true;
    end

    % ===============================================================
    % Load the records
    % ===============================================================
    
    fprintf('Loading data from %s\n', fn);
    [fid,msg] = fopen(fn, 'r');
    if fid<0
      error('File open failed (%s)\n%s',fn, msg);
    end
    
    dropped_record_first_rec_check = rec;
    for rec = rec : rec+num_rec-1
      if param.load.offset{adc_idx}(rec) == -2^31
        % This indicates that this record was dropped and no data exists
        % ideally we should interpolate and fill in... what actually
        % happens is that the data gets filled in with the last good data
        % except if it is the first record being loaded in which case
        % it is filled with zeros.
        dropped_record = true;
      else
        dropped_record = false;
      end
      
      if ~dropped_record
        offset = param.load.offset{adc_idx}(rec) + length(last_record);
        fseek(fid, offset, 'bof');
        if HEADER_SIZE+param.load.rec_data_size-length(last_record) > 0
          rec_data = [last_record; fread(fid, HEADER_SIZE+param.load.rec_data_size-length(last_record), 'uint8')];
        else
          rec_data = last_record(1 : HEADER_SIZE+param.load.rec_data_size);
        end
        last_record = uint8([]);
      else
        % Dropped record... no data exists
        if rec == dropped_record_first_rec_check
          rec_data = zeros(HEADER_SIZE+param.load.rec_data_size,1,'uint8');
        else
          % Just use the data from the last record
        end
      end
      
      % ===============================================================
      % Process record
      % ===============================================================
      for accum_idx = 1:length(accum(adc).wf)
        wf = accum(adc).wf(accum_idx);
        % Convert little endian load into big endian values
        tmp = single(rec_data(1 + 2*(0:wfs(wf).Nt_raw-1) + HEADER_SIZE + wfs(wf).offset)) * 256 ...
          + single(rec_data(2 + 2*(0:wfs(wf).Nt_raw-1) + HEADER_SIZE + wfs(wf).offset));
        % Convert to volts, remove DC-bias, and apply trim
        mean_tmp = mean(tmp(1+param.proc.trim_vals(1):end-param.proc.trim_vals(2)));
        tmp([1:param.proc.trim_vals(1) end-param.proc.trim_vals(2)+1:end]) = mean_tmp;
        tmp = (tmp-mean_tmp) * wfs(wf).quantization_to_V;
        % Accumulate (presum)
        if num_accum == 0
          accum(adc).data{accum_idx} = tmp;
        else
          accum(adc).data{accum_idx} = accum(adc).data{accum_idx} + tmp;
        end
      end

      num_accum = num_accum + 1;
      if num_accum >= param.proc.presums
        num_accum = 0;
        out_idx = out_idx + 1;    
        for accum_idx = 1:length(accum(adc).wf)
          wf = accum(adc).wf(accum_idx);
          img_idx = accum(adc).img_idx(accum_idx);
          wf_adc_idx = accum(adc).wf_adc_idx(accum_idx);
          
          if param.proc.pulse_rfi.en
            pdata = abs(accum(adc).data{accum_idx}).^2;
            inc_ave = param.proc.pulse_rfi.inc_ave;
            thresh_scale = param.proc.pulse_rfi.thresh_scale;
            thresh = filter(ones(inc_ave,1)/inc_ave,1,pdata);
            %thresh = medfilt1(double(thresh),inc_ave*4+1);
            thresh = [thresh(1+(inc_ave-1)/2:end,:); repmat(thresh(end),[(inc_ave-1)/2 size(thresh,2)])];
            bad_idxs = find(pdata > thresh_scale*thresh);
            for bad_idx = bad_idxs.'
              if bad_idx >= 7 && bad_idx <= length(pdata)-6
                vals = pdata(max(1,bad_idx-6):min(bad_idx+6,end));
                comp_vals = pdata(max(1,bad_idx-30):min(bad_idx+30,end));
                thresh2 = median(comp_vals) * 10^(10/10);
                bad_idxs2 = find(vals > thresh2);
                accum(adc).data{accum_idx}(bad_idx-6 + (bad_idxs2-1)) = 0;
              end
            end
            if 0
              figure(1); clf;
              plot(lp(pdata))
              hold on
              plot(lp(thresh_scale*thresh),'r');
              plot(lp(abs(accum(adc).data{accum_idx}).^2),'g');
              hold off;
              keyboard;              
            end
          end
          % Apply channel compensation
          chan_equal = 10.^(param.radar.wfs(wf).chan_equal_dB(param.radar.wfs(wf).rx_paths(adc))/20) ...
            .* exp(j*param.radar.wfs(wf).chan_equal_deg(param.radar.wfs(wf).rx_paths(adc))/180*pi);
          accum(adc).data{accum_idx} = accum(adc).data{accum_idx}/chan_equal;
            accum(adc).data{accum_idx} = accum(adc).data{accum_idx}/wfs(wf).adc_gains(adc);
          if param.proc.pulse_comp
            % ===========================================================
            % Do pulse compression
            % Apply blank (only should enable if sidelobe problems present)
            if ~isempty(wfs(wf).blank)
              % accum(adc).data{accum_idx}(wfs(wf).time_raw>wfs(wf).blank(1) & wfs(wf).time_raw<wfs(wf).blank(2)) = 0;
              accum(adc).data{accum_idx}(wfs(wf).time_raw-param.radar.wfs(wf).Tsys(adc) <= param.surface(rec) + wfs(wf).blank) = 0;
            end
            % Apply matched filter
            % Zero pad front: (the standard)
            accum(adc).data{accum_idx} = fft([zeros(wfs(wf).pad_length,1); accum(adc).data{accum_idx}], wfs(wf).Nt_pc);
            % Zero pad end: (debug only)
            %accum(adc).data{accum_idx} = fft(accum(adc).data{accum_idx}, wfs(wf).Nt_pc);
            
            % Apply matched filter and transform back to time domain
            accum(adc).data{accum_idx} = ifft(accum(adc).data{accum_idx} .* wfs(wf).ref{adc});
            
            if param.proc.ft_dec
              % Digital down conversion and decimation
              accum(adc).data{accum_idx} = accum(adc).data{accum_idx}.*exp(-1i*2*pi*wfs(wf).fc*wfs(wf).time_raw);
              accum(adc).data{accum_idx} = resample(double(accum(board+1).data{accum_idx}), param.wfs(1).ft_dec(1), param.wfs(1).ft_dec(2));
            end
            
          elseif param.proc.ft_dec
            accum(adc).data{accum_idx} = fft(accum(adc).data{accum_idx},wfs(wf).Nt_raw);
            accum(adc).data{accum_idx} = ifft(accum(adc).data{accum_idx}(wfs(wf).freq_inds));
            accum(adc).data{accum_idx} = accum(adc).data{accum_idx}.*exp(-1i*2*pi*wfs(wf).fc*wfs(wf).time_raw);
            accum(adc).data{accum_idx} = resample(double(accum(board+1).data{accum_idx}), param.wfs(1).ft_dec(1), param.wfs(1).ft_dec(2));
          end
          if ~param.load.wf_adc_comb.en
            % Regular wf-adc pair loading: no combining wf-adc pairs in fast-time
            if param.proc.combine_rx
              g_data{img_idx}(:,out_idx) = g_data{img_idx}(:,out_idx) + accum(adc).data{accum_idx} / param.proc.presums / size(param.load.imgs{img_idx},1);
            else
              g_data{img_idx}(:,out_idx,wf_adc_idx) = accum(adc).data{accum_idx} / param.proc.presums;
            end
            
          else
            % Combine wf-adc pairs in fast-time
            if accum(adc).img_comb_idx(accum_idx) == 1
              tmp2{wf_adc_idx} = zeros(param.load.wf_adc_comb.Nt_orig,1);
              tmp2{wf_adc_idx}(1:param.load.wf_adc_comb.rbins(1,out_idx)) ...
                = accum(adc).data{accum_idx}(1:param.load.wf_adc_comb.rbins(1,out_idx)) / param.proc.presums;
              %               g_data{img_idx}(1:param.load.wf_adc_comb.rbins(1,out_idx),out_idx,wf_adc_idx) ...
              %                 = accum(board+1).data{accum_idx}(1:param.load.wf_adc_comb.rbins(1,out_idx)) / param.proc.presums;
            elseif accum(adc).img_comb_idx(accum_idx) == 2
              tmp2{wf_adc_idx}(param.load.wf_adc_comb.rbins(1,out_idx)+1:end) ...
                = accum(adc).data{accum_idx}(param.load.wf_adc_comb.rbins(2,out_idx):end) / param.proc.presums;
              g_data{img_idx}(:,out_idx,wf_adc_idx) = tmp2{wf_adc_idx}(param.load.wf_adc_comb.keep_bins);
              %               g_data{img_idx}(param.load.wf_adc_comb.rbins(1,out_idx)+1:end,out_idx,wf_adc_idx) ...
              %                 = accum(board+1).data{accum_idx}(param.load.wf_adc_comb.rbins(2,out_idx):end) / param.proc.presums;
            end
          end
        end
      end
    end

    if get_last_record
      rec = rec + 1;
      if double(param.load.offset{adc_idx}(rec)) < 0 && param.load.offset{adc_idx}(rec) ~= -2^31
        % Get the last record
        fseek(fid,param.load.offset{adc_idx}(rec),'eof');
        last_record = fread(fid,-double(param.load.offset{adc_idx}(rec)),'uint8');
      else
        get_last_record = false;
        last_record = uint8([]);
      end
    end

    fclose(fid);
    
    % ===============================================================
  end
end

return;

% ===================================================================
% ===================================================================
% load_mcords_data.m Examples
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
% Builds an accumulator structure for load_mcords_data. This is required
% to run load_mcords_data.
%
% imgs = cell vector of imgs (each entry is a separate wf_adc_list)
%   Each wf_adc_list is an Nx2 array where N is the number of channels,
%   the first column is the waveform, and the second column is the adc.
%   Absolute index of wf and adc are used in this array.
% adcs = vector of valid adcs (absolute index of each adc)
%
% accum = structure vector that helps load_mcords_data load data ("accumulator")
%   Each entry in the structure vector corresponds to a specific adc.
%     accum(adc)
%   Each entry contains lists of how to accumulate and store data for the
%   adc when loading from load_mcords_data.  load_mcords_data first
%   accumulates the data (presums):
%     accum(adc).data{accumulator-instance}
%   Once the presums are finished, the data is processed and stored in
%   the output variable:
%     g_data{img_idx}(fast-time,slow-time,wf_adc_idx)
%   The three fields in accum(adc) are all length Kx1 where K is the number
%   of accumulator instances.  The number of accumulator instances is
%   determined by the length of the fields.
%  .wf = K x 1 vector indicating which waveform this instance is pulled from
%  .wf_adc_idx = an index in the final output array unless adcs are
%    combined in which case size(g_data{:}, 3) == 1.
%  .img_idx = an index in the final output array
%
% Examples: At the bottom of this file
%
% Author: John Paden
%
% See also: load_mcords_data

% Verify that all adc entries are valid
for img_idx = 1:length(imgs)
  for wf_adc_idx = 1:size(imgs{img_idx},1)
    if ~ismember(imgs{img_idx}(wf_adc_idx,2), adcs)
      error('ADC %d is invalid', imgs{img_idx}(wf_adc_idx,2));
    end
  end
end

% Build accum structure
for adc = adcs
  accum(adc).wf = [];
  accum(adc).wf_adc_idx = [];
  accum(adc).img_idx = [];
  accum(adc).img_comb_idx = [];
  for img_idx = 1:length(imgs)
    for wf_adc_idx = 1:size(imgs{img_idx},1)
      for adc_column = 2:2:size(imgs{img_idx},2)
        if imgs{img_idx}(wf_adc_idx,2) == adc
          accum(adc).wf(end+1) = imgs{img_idx}(wf_adc_idx,adc_column-1);
          accum(adc).wf_adc_idx(end+1) = wf_adc_idx;
          accum(adc).img_idx(end+1) = img_idx;
          accum(adc).img_comb_idx(end+1) = adc_column/2;
        end
      end
    end
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
for adc=1:8
  adc
  accum(adc)
end

% Example 3
imgs = {[1 1; 1 2; 1 3; 1 4; 1 5; 1 6; 1 7; 1 8], ...
  [2 1; 2 2; 2 3; 2 4; 2 5; 2 6; 2 7; 3 1; 3 2; 3 3; 3 4; 3 5; 3 6; 3 7; 3 8]};
adcs = [1 2 3 4 5 6 7 8];
accum = build_img_load_struct(imgs, adcs);
for adc=1:8
  adc
  accum(adc)
end



