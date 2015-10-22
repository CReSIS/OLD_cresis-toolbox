function [param] = load_accum2_data(param)
% [param] = load_accum2_data(param)
%
% ----------------------------------------------------------------
%
% TODO:
% - Improve readability
% - Sometimes rec_data comes out as a double; resolve this problem
%
% ----------------------------------------------------------------
%
% Function for loading accum2 data (based on mcords2 so there is some legacy
% stuff from mcords2 in the code that has nothing to do with the current accum2
% such as # of ADC "boards"). Supports loading over multiple files,
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
%  .file_idx = file indexes for each record to load (from records)
%   .offset = offset into each record to load (from records)
%   .rec_data_size = size of all the data (not header) bytes in the record
%   .filenames = cell array of strings which file_idx references
%   .imgs = cell vector of wf_adc_list's. A wf_adc_list is an Nx2 array
%     of N (wf,adc) pairs from which the image will be formed. The first
%     column is the waveform and the second column is the adc. These
%     wf/adc references are absolute (i.e. not relative indices)
%     To combine wf and adc into a single I+jQ channel, there are several
%     options:
%     1) Q is the same ADC, but the following wf, use this syntax [j*wf adc]
%     2) Q is the same WF, but the following adc, use this syntax [wf j*adc]
%     You can change the sign to -j*wf or -j*adc and it will change the
%     sign to I-jQ.  Note that IQ pairs on different ADCs must be on the
%     same board (i.e. you can combine 1&2, 2&3, 3&4, but not 4&5 since
%     4 and 5 are on boards 0 and 1 respectively).
%     IQ Note: you only specify the I-channel in .imgs (the Q-channel entry
%     is automatically added by the function build_img_load_struct).
%     Limitation: Each wf/adc pair must have equal length Nt.
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
%    This array is usually created with load_mcords2_wfs and is explained
%    in that function.
%
% Outputs:
% param: updates to this structure are made to allow loading across
%   file boundaries
%
% Global Outputs:
% g_data: cell vector of 3-D matrices of radar data (voltage scale)
%   g_data{img_idx}(bin,rline,wf_adc_idx)
%     img_idx = image index
%     bin = fast-time index (range bin)
%     rline = slow-time index (range line)
%     wf_adc_idx = waveform-adc pair index
%
% Examples at bottom of file
%
% Authors: John Paden
%
% See also get_heights.m, csarp.m, load_mcords.m, load_mcords_wfs.m

load_accum2_data_tstart = tic;

HEADER_SIZE = 0;
WF_HEADER_SIZE = 64;
WF_FOOTER_SIZE = 8;
REC_BLOCK_SIZE = 20e6;
NUM_BOARDS = 1;
NUM_ADCS_PER_BOARD = 1;

global g_data;
wfs = param.wfs;

physical_constants;

% Initialize accumulator
accum = build_img_load_struct(param.load.imgs, param.load.adcs);

% ===================================================================
% Preallocate data matrix
total_rec = length(param.load.file_idx{1});
Nx = floor(total_rec/param.proc.presums);
if ~iscell(g_data)
  g_data = cell(size(param.load.imgs));
end
for img_idx = 1:length(param.load.imgs)
  wf = abs(param.load.imgs{img_idx}(1,1));
  Nt = param.wfs(wf).Nt;
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
boards = unique(floor((param.load.adcs-1)/NUM_BOARDS));
for board_idx = 1:length(boards)
  board = boards(board_idx);
  % Since all ADCs for this board have the same offset and file_idx info
  % we just need to grab the first valid adc_idx:
  adc_idx = find(floor((param.load.adcs-1)/NUM_BOARDS) == board, 1);
  
  % ********************** %
  if param.load.offset{adc_idx}(1) < 0 && param.load.offset{adc_idx}(1) ~= -2^31
    % Get the last incomplete record from the previous file
    fn_idx = param.load.file_idx{adc_idx}(1) - 1;
    fn = param.load.filenames{adc_idx}{fn_idx};
    [fid,msg] = fopen(fn,'r','ieee-be');
    if fid<0
      error('File open failed (%s)\n%s',fn, msg);
    end
    fseek(fid,param.load.offset{adc_idx}(1),'eof');
    last_record = fread(fid,-double(param.load.offset{adc_idx}(1))/2,'int16');
    fclose(fid);
    get_last_record = true;
  else
    % Do not need the last incomplete record from the previous file
    last_record = int16([]);
    get_last_record = false;
  end
  
  num_accum = 0;
  out_idx = 0;
  % All header bytes (i.e. no data samples since param.load.rec_data_size
  % contains this information)
  all_hdr_size = HEADER_SIZE + length(wfs)*(WF_HEADER_SIZE + WF_FOOTER_SIZE);
  % Compute number of records to load per job
  NUM_RECS_LOAD = ceil(REC_BLOCK_SIZE / (all_hdr_size+NUM_ADCS_PER_BOARD*param.load.rec_data_size));
  rec = 1;
  
  while rec < total_rec;
    % Get the filename
    fn_idx = param.load.file_idx{adc_idx}(rec);
    fn = param.load.filenames{adc_idx}{fn_idx};
    
    % Get number of records to load from this file
    %   - will return empty when this is the only file left to load from
    %   - otherwise returns the number of records left to load from this
    %     file
    num_rec = find(param.load.file_idx{adc_idx}(rec:total_rec ) > fn_idx, 1) - 1;
    
    % Check to see if we crossed a file boundary
    if isempty(num_rec)
      % All the records left to load are in this file
      num_rec = total_rec-rec+1;
      get_last_record = false;
    else
      % We will be crossing a file boundary, so grab the last partial
      % record in this file (since it will be used when we start loading
      % in the next file)
      get_last_record = true;
    end
    
    if (num_rec > NUM_RECS_LOAD)
      % Number of records to load is more than the block size, so
      % truncate number to read
      num_rec = NUM_RECS_LOAD;
      get_last_record = false;
    end
    
    % ===============================================================
    % Load the records
    % ===============================================================
    
    fprintf('Loading data from %s (%.01f sec)\n', fn, toc(load_accum2_data_tstart));
    [fid,msg] = fopen(fn, 'r','ieee-be');
    if fid <= 0
      error('File open failed (%s)\n%s',fn, msg);
    end
    
    if(param.load.offset{adc_idx}(rec+num_rec-1) ~= -2^31 && param.load.offset{adc_idx}(rec) ~= -2^31)
      actual_block_size = all_hdr_size+NUM_ADCS_PER_BOARD*param.load.rec_data_size ...
        + param.load.offset{adc_idx}(rec+num_rec-1) ...
        - param.load.offset{adc_idx}(rec) - 2*length(last_record);
    else
      n = 0;
      m = 0;
      %Find correct offsets to use if the offset at the beginning/end of
      %the block is invalid
      while(param.load.offset{adc_idx}(rec+num_rec-1+n) == -2^31 || param.load.offset{adc_idx}(rec+m) == -2^31)
        if(param.load.offset{adc_idx}(rec+num_rec-1+n) == -2^31)
          n = n + 1;
        end
        if (param.load.offset{adc_idx}(rec+m) == -2^31)
          m = m + 1;
        end
      end
      
      actual_block_size = all_hdr_size+NUM_ADCS_PER_BOARD*param.load.rec_data_size ...
        + param.load.offset{adc_idx}(rec+num_rec-1+n) ...
        - param.load.offset{adc_idx}(rec+m) - 2*length(last_record);
    end
    
    %actual_block_size may be less than the ideal block size due to file
    %boundaries, but file access is minimized nonetheless
    
    offset = param.load.offset{adc_idx}(rec) + 2*length(last_record);
    fseek(fid, offset, 'bof');
    if all_hdr_size+NUM_ADCS_PER_BOARD*param.load.rec_data_size - length(last_record) > 0
      %Case where last_record is null or incomplete
      block_data = [last_record; fread(fid, (actual_block_size)/2, 'int16')];
    else
      %Case where last_record is a full record
      block_data = last_record;
    end
    
    first_rec = rec;
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
        rec_data = block_data((param.load.offset{adc_idx}(rec)/2-param.load.offset{adc_idx}(first_rec)/2+1) ...
          + (0 : (all_hdr_size+NUM_ADCS_PER_BOARD*param.load.rec_data_size)/2-1 ));
        last_record = int16([]);
      else
        % Dropped record... no data exists
        if rec == 1
          rec_data = zeros((all_hdr_size+param.load.rec_data_size)/2,1,'int16');
        else
          % Just use the data from the last record
        end
      end
      
      % ===============================================================
      % Process records
      % ===============================================================
      for accum_idx = 1:length(accum(board+1).wf)
        adc = accum(board+1).adc(accum_idx);
        rel_adc = mod(adc-1,NUM_BOARDS)+1;
        wf = accum(board+1).wf(accum_idx);
        % Convert little endian load into big endian values
        partial_hdr_size = HEADER_SIZE + wf*WF_HEADER_SIZE + (wf-1)*WF_FOOTER_SIZE;
        tmp = single(rec_data(1+mod(rel_adc-1,NUM_BOARDS) + ...
          NUM_ADCS_PER_BOARD*(0:wfs(wf).Nt_raw-1) + partial_hdr_size/2 + NUM_ADCS_PER_BOARD*wfs(wf).offset/2));
        % Convert to volts, remove DC-bias, and apply trim
        mean_tmp = mean(tmp(1+param.proc.trim_vals(1):end-param.proc.trim_vals(2)));
        tmp([1:param.proc.trim_vals(1) end-param.proc.trim_vals(2)+1:end]) = mean_tmp;
        tmp = (tmp-mean_tmp) * wfs(wf).quantization_to_V;
        % Accumulate (presum)
        if num_accum == 0
          accum(board+1).data{accum_idx} = tmp;
        else
          accum(board+1).data{accum_idx} = accum(board+1).data{accum_idx} + tmp;
        end
      end
      
      num_accum = num_accum + 1;
      if num_accum >= param.proc.presums
        num_accum = 0;
        out_idx = out_idx + 1;
        for accum_idx = 1:length(accum(board+1).wf)
          adc = accum(board+1).adc(accum_idx);
          wf = accum(board+1).wf(accum_idx);
          img_idx = accum(board+1).img_idx(accum_idx);
          wf_adc_idx = accum(board+1).wf_adc_idx(accum_idx);
          iq_mode = accum(board+1).iq_mode(accum_idx);
          
          if param.proc.pulse_rfi.en
            pdata = abs(accum(board+1).data{accum_idx}).^2;
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
                accum(board+1).data{accum_idx}(bad_idx-6 + (bad_idxs2-1)) = 0;
              end
            end
          end
          % Apply channel compensation
          chan_equal = 10.^(param.radar.wfs(wf).chan_equal_dB(param.radar.wfs(wf).rx_paths(adc))/20) ...
            .* exp(j*param.radar.wfs(wf).chan_equal_deg(param.radar.wfs(wf).rx_paths(adc))/180*pi);
          accum(board+1).data{accum_idx} = accum(board+1).data{accum_idx}/chan_equal;
          accum(board+1).data{accum_idx} = accum(board+1).data{accum_idx}/wfs(wf).adc_gains(adc);
          % Combine I&Q channels if necessary
          if iq_mode == 1
            continue;
          elseif abs(iq_mode) >= 2
            accum(board+1).data{accum_idx} = accum(board+1).data{accum_idx-1} + 1i*sign(iq_mode)*accum(board+1).data{accum_idx};
          end
          if param.proc.pulse_comp
            % ===========================================================
            % Do pulse compression
            % Apply blank (only should enable if sidelobe problems present)
            if ~isempty(wfs(wf).blank)
              accum(board+1).data{accum_idx}(wfs(wf).time_raw>wfs(wf).blank(1) & wfs(wf).time_raw<wfs(wf).blank(2)) = 0;
            end
            % Apply matched filter
            % Zero pad front: (the standard)
            accum(board+1).data{accum_idx} = fft([zeros(wfs(wf).pad_length,1); accum(board+1).data{accum_idx}]);
            % Zero pad end: (debug only)
            %accum(board+1).data{accum_idx} = fft(accum(board+1).data{accum_idx}, wfs(wf).Nt_pc);
            accum(board+1).data{accum_idx} = ifft(accum(board+1).data{accum_idx}(wfs(wf).freq_inds) ...
              .* wfs(wf).ref{adc}(wfs(wf).freq_inds));
            if wfs(wf).dc_shift ~= 0
              % Correct for small frequency offset caused by selecting bins from
              % frequency domain as the method for down conversion
              accum(board+1).data{accum_idx} = accum(board+1).data{accum_idx}.*exp(-1i*2*pi*wfs(wf).dc_shift*wfs(wf).time);
            end
          elseif param.proc.ft_dec
            accum(board+1).data{accum_idx} = fft(accum(board+1).data{accum_idx},wfs(wf).Nt_raw);
            accum(board+1).data{accum_idx} = ifft(accum(board+1).data{accum_idx}(wfs(wf).freq_inds));
            if wfs(wf).dc_shift ~= 0
              % Correct for small frequency offset caused by selecting bins from
              % frequency domain as the method for down conversion
              accum(board+1).data{accum_idx} = accum(board+1).data{accum_idx}.*exp(-1i*2*pi*wfs(wf).dc_shift*wfs(wf).time);
            end
          end
          if param.proc.combine_rx
            g_data{img_idx}(:,out_idx) = g_data{img_idx}(:,out_idx) + accum(board+1).data{accum_idx} / param.proc.presums / size(param.load.imgs{img_idx},1);
          else
            g_data{img_idx}(:,out_idx,wf_adc_idx) = accum(board+1).data{accum_idx} / param.proc.presums;
          end
        end
      end
    end
    
    rec = rec + 1;
    if get_last_record
      if double(param.load.offset{adc_idx}(rec)) < 0 && param.load.offset{adc_idx}(rec) ~= -2^31
        % Get the last record
        fseek(fid,param.load.offset{adc_idx}(rec),'eof');
        last_record = fread(fid,-double(param.load.offset{adc_idx}(rec))/2,'int16');
      else
        get_last_record = false;
        last_record = int16([]);
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
% accum = structure vector that helps load_accum2_data load data ("accumulator")
%   Each entry in the structure vector corresponds to a specific board.
%     accum(board)
%   Each entry contains lists of how to accumulate and store data for the
%   board when loading from load_accum2_data.  load_accum2_data first
%   accumulates the data (presums):
%     accum(board).data{accumulator-instance}
%   Once the presums are finished, the data is processed and stored in
%   the output variable:
%     g_data{img_idx}(fast-time,slow-time,wf_adc_idx)
%   The three fields in accum(board) are all length Kx1 where K is the number
%   of accumulator instances.  The number of accumulator instances is
%   determined by the length of the fields.
%  .adc = K x 1 vector indicating which adc this instance is pulled from
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

NUM_BOARDS = 4;

% Verify that all adc entries are valid
for img_idx = 1:length(imgs)
  for wf_adc_idx = 1:size(imgs{img_idx},1)
    if ~ismember(real(imgs{img_idx}(wf_adc_idx,2)), adcs)
      error('ADC %d is invalid', imgs{img_idx}(wf_adc_idx,2));
    end
  end
end

% Create board list
boards = unique(floor((adcs-1)/NUM_BOARDS));

% Build accum structure
for board = boards
  board_adcs = adcs(floor((adcs-1)/NUM_BOARDS) == board);
  accum(1+board).adc = [];
  accum(1+board).wf = [];
  accum(1+board).wf_adc_idx = [];
  accum(1+board).img_idx = [];
  accum(1+board).iq_mode = [];
  for adc = board_adcs
    for img_idx = 1:length(imgs)
      for wf_adc_idx = 1:size(imgs{img_idx},1)
        if abs(imgs{img_idx}(wf_adc_idx,2)) == adc
          if isreal(imgs{img_idx})
            % [wf adc] --> real only waveform
            accum(1+board).adc(end+1) = adc;
            accum(1+board).wf(end+1) = imgs{img_idx}(wf_adc_idx,1);
            accum(1+board).wf_adc_idx(end+1) = wf_adc_idx;
            accum(1+board).img_idx(end+1) = img_idx;
            accum(1+board).iq_mode(end+1) = 0;
          else
            % [j*wf adc] --> use waveform 1 and 2 to form I + j*Q
            % [-j*wf adc] --> use waveform 1 and 2 to form I - j*Q
            % [wf j*adc] --> use adc 1 and 2 to form I + j*Q
            % [wf -j*adc] --> use adc 1 and 2 to form I - j*Q
            if real(imgs{img_idx}(wf_adc_idx,1)) == 0
              % Next waveform is the Q channel
              accum(1+board).adc(end+1) = adc;
              accum(1+board).wf(end+1) = abs(imgs{img_idx}(wf_adc_idx,1));
              accum(1+board).wf_adc_idx(end+1) = wf_adc_idx;
              accum(1+board).img_idx(end+1) = img_idx;
              accum(1+board).iq_mode(end+1) = 1;
              
              accum(1+board).adc(end+1) = adc;
              accum(1+board).wf(end+1) = abs(imgs{img_idx}(wf_adc_idx,1)) + 1;
              accum(1+board).wf_adc_idx(end+1) = wf_adc_idx;
              accum(1+board).img_idx(end+1) = img_idx;
              accum(1+board).iq_mode(end+1) = 2*sign(imag(imgs{img_idx}(wf_adc_idx,1)));
            else
              % Next adc is the Q channel
              accum(1+board).adc(end+1) = adc;
              accum(1+board).wf(end+1) = imgs{img_idx}(wf_adc_idx,1);
              accum(1+board).wf_adc_idx(end+1) = wf_adc_idx;
              accum(1+board).img_idx(end+1) = img_idx;
              accum(1+board).iq_mode(end+1) = 1;
              
              accum(1+board).adc(end+1) = adc+1;
              accum(1+board).wf(end+1) = imgs{img_idx}(wf_adc_idx,1);
              accum(1+board).wf_adc_idx(end+1) = wf_adc_idx;
              accum(1+board).img_idx(end+1) = img_idx;
              accum(1+board).iq_mode(end+1) = 2*sign(imag(imgs{img_idx}(wf_adc_idx,2)));
            end
          end
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
for board = 0:length(accum)-1
  board
  accum(board+1)
end

% Example 3
imgs = {[1 1; 1 2; 1 3; 1 4; 1 5; 1 6; 1 7; 1 8], ...
  [2 1; 2 2; 2 3; 2 4; 2 5; 2 6; 2 7; 3 1; 3 2; 3 3; 3 4; 3 5; 3 6; 3 7; 3 8]};
adcs = [1 2 3 4 5 6 7 8];
accum = build_img_load_struct(imgs, adcs);
for board = 0:length(accum)-1
  board
  accum(board+1)
end



