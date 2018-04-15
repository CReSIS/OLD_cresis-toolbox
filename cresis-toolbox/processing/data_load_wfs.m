function [wfs,state] = data_load_wfs(param, records)
% [wfs,state] = data_load_wfs(param, records)
%
% Loads waveform/adc information into wfs
% Loads state information for loading raw data into state
% Both are required for data_load.m
%
% param: structure with parameter information
%  .records.file_version: raw file version
%  .load.imgs = cell vector of imgs (each entry is a separate wf_adc_list)
%    Each wf_adc_list is an Nx2 array where N is the number of channels,
%    the first column is the waveform, and the second column is the adc.
%    Absolute index of wf and adc are used in this array.
%
% wfs: 
% state: structure vector that helps data_load load data (state of loader)
%   Each entry in the structure vector corresponds to a specific board.
%     accum(board)
%   Each entry contains lists of how to accumulate and store data for the
%   board when loading from load_mcords2_data.  load_mcords2_data first
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
% See also: data_load.m

save('/tmp/data_load_wfs_test.mat')
keyboard
load('/tmp/data_load_wfs_test.mat')

%% Build raw data loading "state" structure
% =========================================================================

% Create list of wf-adc pairs to determine which boards to load
wf_list = [];
adc_list = [];
board_list = [];
for img = 1:length(param.load.imgs)
  for wf_adc_idx= 1:size(param.load.imgs{img},1)
    for adc_column = 2:2:size(imgs{img_idx},2)
      wf_list(end+1) = param.load.imgs{img}(wf_adc_idx,adc_column-1);
      adc_list(end+1) = param.load.imgs{img}(wf_adc_idx,adc_column);
      board_list(end+1) = adc_to_board(param.radar_name,adc_list(end));
    end
  end
end
boards = unique(board_list);

% Populate state structure
state = [];
for board = boards
  state(board).adc = [];
  state(board).wf = [];
  state(board).wf_adc_idx = [];
  state(board).img = [];
  state(board).wf_adc_sum = [];
  state(board).img_comb_idx = [];
  for img = 1:length(param.load.imgs) % For each image img
    for wf_adc_idx = 1:size(param.load.imgs{img},1) % For ach wf-adc pair
      for adc_column = 2:2:size(param.load.imgs{img},2) % For each combined wf-adc pair
        wf = param.load.imgs{img}(wf_adc_idx,adc_column-1); % wf stored in odd columns
        adc = param.load.imgs{img}(wf_adc_idx,adc_column); % adc stored in even columns
        if ~isfield(param.radar.wfs(wf),'wf_adc_sum') || isempty(param.radar.wfs(wf).wf_adc_sum)
          % if wf_adc_sum not specied, then this wf-adc pair is the only
          % one to load
          wf_adc_sum = [wf adc 1];
        else
          % if wf_adc_sum is specified, then extract the list out, each
          % row specifies a wf-adc pair and a complex weight
          % For zero-pi mod with two wf-adc pairs, the weight should be -0.5 or 0.5
          % For IQ mod with two wf-adc pairs, the weights should be 1, j, -j, or -1
          %   [wf adc weight] 
          wf_adc_sum = param.radar.wfs(wf).wf_adc_sum{adc};
        end
        for wf_adc_sum_idx = 1:size(wf_adc_sum,1)
          wf = wf_adc_sum(wf_adc_sum_idx,1);
          adc = wf_adc_sum(wf_adc_sum_idx,2);
          % Add wf-adc pair to state list
          state(board).adc(end+1) = adc;
          state(board).wf(end+1) = wf;
          state(board).wf_adc_idx(end+1) = wf_adc_idx;
          state(board).img(end+1) = img;
          state(board).wf_adc_sum(end+1) = wf_adc_sum(wf_adc_sum_idx,3);
          state(board).img_comb_idx(end+1) = adc_column/2;
        end
      end
    end
  end
end

%% Create wfs structure with waveform information
% =========================================================================
[output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name);
for wf = 1:length(param.radar.wfs)
  
  %% Input checks
  % =======================================================================
  if isfield(param.radar.wfs(wf),'Tpd') && ~isempty(param.radar.wfs(wf).Tpd)
    wfs(wf).Tpd = param.radar.wfs(wf).Tpd;
  elseif any(param.records.file_version == [405 406 410]) % [acords mcords]
    wfs(wf).Tpd = records.settings.wfs(1).wfs(wf).Tpd(1);
  end
  if isfield(param.radar.wfs(wf),'f0') && ~isempty(param.radar.wfs(wf).f0)
    wfs(wf).f0 = param.radar.wfs(wf).f0;
  elseif any(param.records.file_version == [405 406 410]) % [acords mcords]
    wfs(wf).f0 = records.settings.wfs(1).wfs(wf).f0(1);
  end
  if isfield(param.radar.wfs(wf),'f1') && ~isempty(param.radar.wfs(wf).f1)
    wfs(wf).f1 = param.radar.wfs(wf).f1;
  elseif any(param.records.file_version == [405 406 410]) % [acords mcords]
    wfs(wf).f1 = records.settings.wfs(1).wfs(wf).f1(1);
  end
  if isfield(param.radar.wfs(wf),'Tadc_adjust') && ~isempty(param.radar.wfs(wf).Tadc_adjust)
    wfs(wf).Tadc_adjust = param.radar.wfs(wf).Tadc_adjust;
  else
    wfs(wf).Tadc_adjust = 0;
  end
  if isfield(param.radar.wfs(wf),'Tadc') && ~isempty(param.radar.wfs(wf).Tadc)
    wfs(wf).t0    = param.radar.wfs(wf).Tadc + wfs(wf).Tadc_adjust;
  elseif any(param.records.file_version == [405 406 410]) % [acords mcords]
    wfs(wf).t0    = records.settings.wfs(1).wfs(wf).t0(1) + wfs(wf).Tadc_adjust;
  elseif isfield(records.settings.wfs(wf),'t0')
    wfs(wf).t0    = records.settings.wfs(wf).t0 + wfs(wf).Tadc_adjust;
  else
    wfs(wf).t0    = 0 + wfs(wf).Tadc_adjust;
  end
  if isfield(param.radar.wfs(wf),'blank') && ~isempty(param.radar.wfs(wf).blank)
    wfs(wf).blank   = param.radar.wfs(wf).blank;
  else
    wfs(wf).blank   = 0;
  end
  if isfield(param.radar.wfs(wf),'DDC_mode') && ~isempty(param.radar.wfs(wf).DDC_mode)
    wfs(wf).DDC_mode   = param.radar.wfs(wf).DDC_mode;
  else
    wfs(wf).DDC_mode   = 0;
  end
  if param.records.file_version == 410 % mcords
    fs = records_wfs.wfs(1).wfs(1).fs;
  else
    if wfs(wf).DDC_mode == 0
      fs = param.radar.fs;
    else
      fs = param.radar.fs / 2^(1+wfs(wf).DDC_mode);
    end
  end
  if isfield(param.radar.wfs(wf),'zero_pad') && ~isempty(param.radar.wfs(wf).zero_pad)
    wfs(wf).zero_pad   = param.radar.wfs(wf).zero_pad;
  else
    wfs(wf).zero_pad   = 0;
  end
  if isfield(param.radar.wfs(wf),'ft_dec') && ~isempty(param.radar.wfs(wf).ft_dec)
    wfs(wf).ft_dec = param.radar.wfs(wf).ft_dec;
  else
    if strcmpi(radar_type,'fmcw')
      wfs(wf).ft_dec = [1 1];
    else
      [numerator denominator] = rat((wfs(wf).f1 - wfs(wf).f0) / fs);
      wfs(wf).ft_dec = [numerator denominator];
    end
  end
  if isfield(param.radar.wfs(wf),'DDC_freq') && ~isempty(param.radar.wfs(wf).DDC_freq)
    wfs(wf).DDC_freq   = param.radar.wfs(wf).DDC_freq;
  else
    wfs(wf).DDC_freq   = 0;
  end
  if isfield(param.radar.wfs(wf),'presums') && ~isempty(param.radar.wfs(wf).presums)
    wfs(wf).presums = param.radar.wfs(wf).presums;
  elseif any(param.records.file_version == [1:8 405 406 410]) % [acords mcords]
    wfs(wf).presums = records.settings.wfs(1).wfs(wf).presums(1);
  else
    wfs(wf).presums = records.settings.wfs(wf).presums;
  end
  if isfield(param.radar.wfs(wf),'iq_mode') && ~isempty(param.radar.wfs(wf).iq_mode)
    wfs(wf).iq_mode   = param.radar.wfs(wf).iq_mode;
  else
    wfs(wf).iq_mode   = 0;
  end
  if isfield(param.radar.wfs(wf),'zero_pi_mode') && ~isempty(param.radar.wfs(wf).zero_pi_mode)
    wfs(wf).zero_pi_mode   = param.radar.wfs(wf).zero_pi_mode;
  else
    wfs(wf).zero_pi_mode   = 0;
  end
  if isfield(param.radar.wfs(wf),'BW') && ~isempty(param.radar.wfs(wf).BW)
    wfs(wf).BW_window = param.radar.wfs(wf).BW(1);
  else
    wfs(wf).BW_window = abs(wfs(wf).f1 - wfs(wf).f0);
  end
  if any(param.records.file_version == [405 406 410]) % [acords mcords]
    wfs(wf).Nt_raw = records.settings.wfs(1).wfs(wf).num_sam(1);
  elseif isfield(records.settings.wfs(wf),'num_sam')
    wfs(wf).Nt_raw = records.settings.wfs(wf).num_sam(1);
  else
    wfs(wf).Nt_raw = 0;
  end
  
  % Other fields from param.radar files
  % =======================================================================
  wfs(wf).tx_weights = param.radar.wfs(wf).tx_weights;
  wfs(wf).rx_paths = param.radar.wfs(wf).rx_paths;
  wfs(wf).adc_gains = param.radar.wfs(wf).adc_gains;

  %% Compute supporting variables
  % =======================================================================
  wfs(wf).Nt_ref  = floor(wfs(wf).Tpd * fs) + 1;
  wfs(wf).Nt_pc   = wfs(wf).Nt_raw + wfs(wf).Nt_ref + wfs(wf).zero_pad - 1;
  wfs(wf).pad_length = wfs(wf).Nt_pc - wfs(wf).Nt_raw;
  
  %% Quantization to Voltage conversion
  % =======================================================================
  if any(param.records.file_version == [405 406]) % acords
    num_bit_shifts = records.settings.wfs(1).wfs(wf).bit_shifts(1);
  elseif param.records.file_version == 410 % mcords
    num_bit_shifts = max(0,ceil(log(wfs(wf).presums)/log(2)) - 4);
  elseif isfield(records.settings.wfs(wf),'bit_shifts')
    num_bit_shifts = records.settings.wfs(wf).bit_shifts;
  else
    num_bit_shifts = 0;
  end
  wfs(wf).quantization_to_V ...
    = param.radar.Vpp_scale * 2^num_bit_shifts ...
    / (2^14*wfs(wf).presums);
  
end

% offset: bytes of data before this data channel
% rec_data_size: number of bytes of data in the record
if any(param.records.file_version == [405 406]) % [acords]
  wf_num_sam = cell2mat({settings.wfs(1).wfs.num_sam}).';
elseif any(param.records.file_version == [403 407 408]) % [mcords3 mcords5]
  wf_num_sam = cell2mat({settings.wfs.num_sam}).';
end



if any(strcmpi(radar_name,{'acords','hfrds','hfrds2','mcords','mcords2','mcords3','mcords4','mcords5','seaice'}))
  wf = 1;
  wf_offsets(wf) = 0;
  if ~wfs(wf).DDC_mode
    % Real data
    if any(strcmpi(radar_name,{'acords'}))
      num_elem = length(settings.wfs(wf).wfs(wf).adc_gains(1,:));
      rec_data_size = wf_num_sam(wf)*num_elem*sample_size;
    else
      rec_data_size = wf_num_sam(wf)*sample_size;
    end
  else
    % Complex data
    rec_data_size = 2*wf_num_sam(wf)*sample_size;
  end
  
  for wf = 2:length(param.radar.wfs)
    if isfield(param.radar.wfs(wf),'DDC_mode') && ~isempty(param.radar.wfs(wf).DDC_mode)
      wfs(wf).DDC_mode   = param.radar.wfs(wf).DDC_mode;
    else
      wfs(wf).DDC_mode   = 0;
    end
    
    if ~wfs(wf).DDC_mode
      % Real data
      if any(strcmpi(radar_name,{'acords'}))
        num_elem = length(settings.wfs(wf).wfs(wf).adc_gains(1,:));
        wf_offsets(wf) = wf_offsets(wf-1) + wf_num_sam(wf-1)*sample_size;
        rec_data_size = rec_data_size + wf_num_sam(wf)*num_elem*sample_size;
      else
        wf_offsets(wf) = wf_offsets(wf-1) + wf_num_sam(wf-1)*sample_size;
        rec_data_size = rec_data_size + wf_num_sam(wf)*sample_size;
      end
    else
      % Complex data
      wf_offsets(wf) = wf_offsets(wf-1) + 2*wf_num_sam(wf-1)*sample_size;
      rec_data_size = rec_data_size + 2*wf_num_sam(wf)*sample_size;
    end
  end
  
elseif strcmpi(radar_name,'accum2')
  wf_offsets = cumsum([0; wf_num_sam(1:end-1)]) ...
    * sample_size;
  rec_data_size = sum(wf_num_sam) * sample_size;
end

if proc_param.wf_adc_comb.en
  % Determine the waveform with the longest time record. All waveforms will
  % be increased to this length for pulse compression to guarantee alignment
  % of the frequency bins. This zero padding will be placed at the end of the
  % waveform.
  
  for wf = 1:length(param.radar.wfs)
    if isfield(param.radar.wfs(wf),'DDC_mode') && ~isempty(param.radar.wfs(wf).DDC_mode)
      wfs(wf).DDC_mode   = param.radar.wfs(wf).DDC_mode;
    else
      wfs(wf).DDC_mode   = 0;
    end
    if wfs(wf).DDC_mode == 0
      fs = param.radar.fs;
    else
      fs = param.radar.fs / 2^(1+wfs(wf).DDC_mode);
    end
    if isfield(param.radar.wfs(wf),'zero_pad') && ~isempty(param.radar.wfs(wf).zero_pad)
      wfs(wf).zero_pad   = param.radar.wfs(wf).zero_pad;
    else
      wfs(wf).zero_pad   = 0;
    end
    
    Tpd     = param.radar.wfs(wf).Tpd;
    Nt_ref  = floor(Tpd * fs) + 1;
    Nt_raw  = settings.wfs(wf).num_sam;
    Nt_pc(wf)   = Nt_raw + Nt_ref + wfs(wf).zero_pad - 1;
  end
  
  Nt_pc_max = max(Nt_pc);
end

%% Create default values for all waveforms
% =========================================================================
for wf = 1:length(param.radar.wfs)
  adc_idx = 1;
  
  if any(strcmpi(radar_name,{'acords'}))
    settings.wfs(wf).bit_shifts = settings.wfs(1).wfs(wf).bit_shifts(1);
    settings.wfs(wf).presums = settings.wfs(1).wfs(wf).presums(1);
    settings.wfs(wf).num_sam = settings.wfs(1).wfs(wf).num_sam(1);
    settings.wfs(wf).t0 = settings.wfs(1).wfs(wf).t0(1);
    settings.wfs(wf).prf = settings.wfs(1).wfs(wf).prf(1);
    settings.wfs(wf).f0 = settings.wfs(1).wfs(wf).f0(1);
    settings.wfs(wf).f1 = settings.wfs(1).wfs(wf).f1(1);
    settings.wfs(wf).wf_gen_clk = settings.wfs(1).wfs(wf).wf_gen_clk(1);
    settings.wfs(wf).daq_clk = settings.wfs(1).wfs(wf).daq_clk(1);
    settings.wfs(wf).Tpd = settings.wfs(1).wfs(wf).Tpd(1);
    settings.wfs(wf).tx_win = settings.wfs(1).wfs(wf).tx_win(1);
    settings.wfs(wf).blank = settings.wfs(1).wfs(wf).blank(1);
    settings.wfs(wf).adc_gains = settings.wfs(1).wfs(wf).adc_gains(1,:);
    if param.records.file_version == 406
      settings.wfs(wf).elem_slots = settings.wfs(1).wfs(wf).elem_slots(1,:);
    end
  end
  
  wfs(wf).Tpd     = param.radar.wfs(wf).Tpd;
  wfs(wf).f0      = param.radar.wfs(wf).f0;
  wfs(wf).f1      = param.radar.wfs(wf).f1;
  if isfield(param.radar.wfs(wf),'Tadc_adjust') && ~isempty(param.radar.wfs(wf).Tadc_adjust)
    Tadc_adjust = param.radar.wfs(wf).Tadc_adjust;
  else
    Tadc_adjust = 0;
  end
  if isfield(param.radar.wfs(wf),'Tadc') && ~isempty(param.radar.wfs(wf).Tadc)
    wfs(wf).t0    = param.radar.wfs(wf).Tadc + Tadc_adjust;
  else
    wfs(wf).t0    = settings.wfs(wf).t0 + Tadc_adjust;
  end
  if isfield(param.radar.wfs(wf),'blank') && ~isempty(param.radar.wfs(wf).blank)
    wfs(wf).blank   = param.radar.wfs(wf).blank;
  else
    wfs(wf).blank   = [];
  end
  if isfield(param.radar.wfs(wf),'DDC_mode') && ~isempty(param.radar.wfs(wf).DDC_mode)
    wfs(wf).DDC_mode   = param.radar.wfs(wf).DDC_mode;
  else
    wfs(wf).DDC_mode   = 0;
  end
  if wfs(wf).DDC_mode == 0
    fs = param.radar.fs;
  else
    fs = param.radar.fs / 2^(1+wfs(wf).DDC_mode);
  end
  if isfield(param.radar.wfs(wf),'zero_pad') && ~isempty(param.radar.wfs(wf).zero_pad)
    wfs(wf).zero_pad   = param.radar.wfs(wf).zero_pad;
  else
    wfs(wf).zero_pad   = 0;
  end
  if isfield(param.radar.wfs(wf),'ft_dec') && ~isempty(param.radar.wfs(wf).ft_dec)
    wfs(wf).ft_dec = param.radar.wfs(wf).ft_dec;
  else
    [numerator denominator] = rat((wfs(wf).f1 - wfs(wf).f0) / fs);
    wfs(wf).ft_dec = [numerator denominator];
  end
  if isfield(param.radar.wfs(wf),'DDC_freq') && ~isempty(param.radar.wfs(wf).DDC_freq)
    wfs(wf).DDC_freq   = param.radar.wfs(wf).DDC_freq;
  else
    wfs(wf).DDC_freq   = 0;
  end
  if isfield(param.radar.wfs(wf),'DC_adjust') && ~isempty(param.radar.wfs(wf).DC_adjust)
    tmp = load(fullfile(ct_filename_out(param,'noise','',1), ...
      param.radar.wfs(wf).DC_adjust),'DC_adjust');
    for adc_idx = 1:length(adcs)
      adc = adcs(adc_idx);
      wfs(wf).DC_adjust(adc_idx) = tmp.DC_adjust(adc);
    end
  else
    wfs(wf).DC_adjust   = zeros(size(adcs));
  end
  wfs(wf).presums = settings.wfs(wf).presums;
  wfs(wf).Nt_ref  = floor(wfs(wf).Tpd * fs) + 1;
  wfs(wf).Nt_raw  = settings.wfs(wf).num_sam;
  wfs(wf).Nt_pc   = wfs(wf).Nt_raw + wfs(wf).Nt_ref + wfs(wf).zero_pad - 1;
  wfs(wf).pad_length = wfs(wf).Nt_pc - wfs(wf).Nt_raw;
  wfs(wf).offset  = wf_offsets(wf);
  if proc_param.wf_adc_comb.en
    wfs(wf).Nt_pc = Nt_pc_max;
  end
  
  if isfield(param.radar.wfs(wf),'conjugate') && ~isempty(param.radar.wfs(wf).conjugate)
    wfs(wf).conjugate   = param.radar.wfs(wf).conjugate;
  else
    wfs(wf).conjugate   = 0;
  end

  if isfield(param.radar.wfs(wf),'BW') && ~isempty(param.radar.wfs(wf).BW)
    BW_window = param.radar.wfs(wf).BW(1);
  else
    BW_window = abs(wfs(wf).f1 - wfs(wf).f0);
  end
  
  % ===================================================================
  % Quantization to Voltage conversion
  wfs(wf).quantization_to_V ...
    = param.radar.Vpp_scale * 2^settings.wfs(wf).bit_shifts ...
    / (2^14*wfs(wf).presums);
  
  % Other fields from param.radar files
  wfs(wf).tx_weights = param.radar.wfs(wf).tx_weights;
  wfs(wf).rx_paths = param.radar.wfs(wf).rx_paths;
  wfs(wf).adc_gains = param.radar.wfs(wf).adc_gains;
  
  % ===================================================================
  % Create reference waveform
  Tpd = wfs(wf).Tpd;
  f0  = wfs(wf).f0;
  f1  = wfs(wf).f1;
  Nt  = wfs(wf).Nt_ref;
  t0  = wfs(wf).t0;
  
  dt = 1/fs;
  BW = f1-f0;
  time = (0:dt:(Nt-1)*dt).';
  alpha = BW / Tpd;
  fc = (f0 + f1)/2;
  
  if wfs(wf).DDC_mode == 0
    %% DDC Disabled or No DDC
    ref_function = exp(1i*2*pi*f0*time + 1i*pi*alpha*time.^2);
  else
    %% DDC Enabled
    ref_function = exp(1i*2*pi*(f0 - wfs(wf).DDC_freq)*time + 1i*pi*alpha*time.^2);
  end
  if any(strcmpi(radar_name,{'acords'}))
    Htukeywin = hamming(Nt);
  else
    Htukeywin = tukeywin(Nt+2,param.radar.wfs(wf).tukey);
    Htukeywin = Htukeywin(2:end-1);
  end
  if proc_param.ft_wind_time && ~isempty(proc_param.ft_wind)
    ref = Htukeywin.*proc_param.ft_wind(Nt).*ref_function;
  else
    ref = Htukeywin.*ref_function;
  end
  
  % Apply receiver delays to reference function
  Nt = wfs(wf).Nt_pc;
  if isempty(Nt)
    warning('Undefined waveform %d: skipping waveform.', wf);
    continue;
  end
  df = 1/(Nt*dt);
  if wfs(wf).DDC_mode == 0
    freq = fs*floor(fc/fs) + (0:df:(Nt-1)*df).';
  else
    freq = wfs(wf).DDC_freq + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
  end
  for adc = adcs
    ref_fn_name = char(param.radar.wfs(wf).ref_fn);
    ref_fn_name = regexprep(ref_fn_name,'%w',sprintf('%.0f',wf));
    ref_fn_name = regexprep(ref_fn_name,'%a',sprintf('%.0f',adc));
    ref_fn = fullfile(ct_filename_out(param,'noise','',1), [ref_fn_name '.mat']);
    
    if isempty(ref_fn_name) || ~exist(ref_fn,'file')
      wfs(wf).ref{adc} = conj(fft(ref,Nt) ...
        .* exp(-1i*2*pi*freq*param.radar.wfs(wf).Tsys(param.radar.wfs(wf).rx_paths(adc))) );
      wfs(wf).ref_windowed(adc) = false;
      
    else
      % Load reference function from collate_deconv.m (e.g. for deconvolution)
      load(ref_fn,'ref_nonnegative','ref_negative','ref_windowed','ref_window');
      ref_Nt = length(ref_nonnegative)+length(ref_negative);
      if ref_Nt > Nt
        error('Window in ref_fn %s is longer than time axis, increase zero padding to use this reference function or shorten the reference function', ref_fn);
      end
      ref_from_file = [ref_nonnegative; zeros(Nt-ref_Nt,1); ref_negative];
      wfs(wf).ref_windowed(adc) = ref_windowed;
      
      if ref_windowed && ~isequal(ref_window,proc_param.ft_wind)
        error('Window in ref_fn %s does not match ft_wind parameter', ref_fn);
      end
      
      ref_from_file = ref_from_file ./ abs(max(ref_from_file));
      wfs(wf).ref{adc} = conj(fft(ref_from_file,Nt) ...
        .* exp(-1i*2*pi*freq*param.radar.wfs(wf).Tsys(param.radar.wfs(wf).rx_paths(adc))) );
    end
    
  end

  if isfield(param.radar.wfs(wf),'gain_fn') && ~isempty(param.radar.wfs(wf).gain_fn)
    for adc = adcs
      gain_fn_name = char(param.radar.wfs(wf).gain_fn);
      gain_fn_name = regexprep(gain_fn_name,'%w',sprintf('%.0f',wf));
      gain_fn_name = regexprep(gain_fn_name,'%a',sprintf('%.0f',adc));
      gain_fn = fullfile(ct_filename_out(param,'noise','',1), [gain_fn_name '.mat']);
      
      wfs(wf).gain(adc) = load(gain_fn);
    end
  end
  
  % ===================================================================
  % Create raw data with zero padding freq axes variables
  if ~proc_param.pulse_comp
    Nt = wfs(wf).Nt_raw;
  end
  
  df = 1/(Nt*dt);
  %freq = round(fc/fs)*fs + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
  freq = fs*floor(fc/fs) + (0:df:(Nt-1)*df).';
  wfs(wf).time_raw = t0 + (0:dt:(Nt-1)*dt).';
  
  %% Create Decimation Information
  if proc_param.ft_dec
    if wfs(wf).DDC_mode ~= 0
      % DDC Enabled
      freq = wfs(wf).DDC_freq + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
    end
    wfs(wf).fc = fc;
  else
    wfs(wf).fc = fs*floor(max(f0,f1)/fs);
  end
  
  %% Create Windowing Information
  if ~proc_param.ft_wind_time && ~isempty(proc_param.ft_wind) ...
      && proc_param.pulse_comp
    ft_wind = zeros(size(freq));
    
    if wfs(wf).DDC_mode == 0
      % DDC Disabled or no DDC
      freq_inds = ifftshift(find(freq >= fc-BW_window/2 & freq <= fc+BW_window/2));
    else
      % DDC Enabled
      freq_inds = find(freq >= fc-BW_window/2 & freq <= fc+BW_window/2);
    end
    [~,sorted_freq_inds] = sort(freq(freq_inds));
    sorted_freq_inds = freq_inds(sorted_freq_inds);
    ft_wind(sorted_freq_inds) = proc_param.ft_wind(length(freq_inds));
    
    for adc = adcs
      if ~wfs(wf).ref_windowed(adc)
        wfs(wf).ref{adc} = wfs(wf).ref{adc} .* ft_wind;
      end
    end
  end
  
  %% Normalize reference function so that it is an estimator
  %  -- Accounts for pulse duration differences
  for adc = adcs
    time_domain_ref = ifft(wfs(wf).ref{adc});
    wfs(wf).ref{adc} = wfs(wf).ref{adc} ...
      ./ dot(time_domain_ref,time_domain_ref);
  end
  
  % ===================================================================
  %% Create output data time/freq axes variables
  Nt = ceil(wfs(wf).Nt_pc*wfs(wf).ft_dec(1)/wfs(wf).ft_dec(2));
  wfs(wf).Nt = Nt;
  if proc_param.ft_dec
    wfs(wf).fs = fs * wfs(wf).ft_dec(1)/wfs(wf).ft_dec(2);
    wfs(wf).dt = 1/wfs(wf).fs;
    wfs(wf).df = wfs(wf).fs/wfs(wf).Nt;
    dt = wfs(wf).dt;
    % Starts at fc goes to fc+BW/2, fc-BW/2 to fc
    wfs(wf).freq = fc + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
  else
    wfs(wf).df = df;
    wfs(wf).dt = 1/(Nt*df);
    wfs(wf).fs = Nt*df;
    dt = wfs(wf).dt;
    % Let ftnz = fast time nyquist zone
    % Starts at ftnz*fs goes to ftnz*fs+fs/2, ftnz*fs-fs/2 to ftnz*fs
    %     wfs(wf).freq = round(fc/fs)*fs ...
    %       + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
    wfs(wf).freq = fs*floor(fc/fs) + (0:df:(Nt-1)*df).';
  end
  wfs(wf).time = t0 + dt*(0:Nt-1).';
  if proc_param.pulse_comp
    % Assuming pulse compression zero pads the front of the waveform, the
    % output will start earlier by an ammount proportional to the zero
    % padding.
    wfs(wf).time = wfs(wf).time - wfs(wf).pad_length / fs;
    
    % Modify reference function so that time vector elements are multiples
    % of dt.
    wfs(wf).time_correction = dt - mod(wfs(wf).time(1),dt);
    wfs(wf).time = wfs(wf).time + wfs(wf).time_correction;
    
    for adc = adcs
      wfs(wf).ref{adc} = wfs(wf).ref{adc} .* exp(1i*2*pi*freq*wfs(wf).time_correction);
    end
  end
  
end

