function success = preprocess_task_cresis(param)
% success = preprocess_task_cresis(param)
%
% Support function for preprocess.m
%
% Example:
% Called from preprocess_task.m
%
% Author: John Paden
%
% See also: run_preprocess.m, preprocess.m, preprocess_task.m,
% preprocess_task_arena.m, preprocess_task_ni.m

% 1. Use what we already have
% 2. Use file version instead of radar name
% 3. Use last read to allow header wraps... or not.

%% Read Headers
% =========================================================================
failed_load = {};
fns_list = cell(size(param.preprocess.daq.board_map));
for board_idx = 1:numel(param.preprocess.daq.board_map)
  %% Read Headers: Filenames
  board = param.preprocess.daq.board_map{board_idx};
  board_folder_name = param.preprocess.board_folder_name;
  board_folder_name = regexprep(board_folder_name,'%b',board);
  
  get_filenames_param = struct('regexp',param.preprocess.file.regexp);
  fns = get_filenames(fullfile(param.preprocess.base_dir,board_folder_name), ...
    param.preprocess.file.prefix, param.preprocess.file.midfix, ...
    param.preprocess.file.suffix, get_filenames_param);
  fns = fns(1:100);
  if param.preprocess.online_mode == 2
    fns = fns(end);
  end
  fns_list{board_idx} = fns;
  
  if isempty(fns)
    error('No files found matching %s*%s*%s', ...
      fullfile(param.preprocess.base_dir,board_folder_name,param.preprocess.file.prefix), ...
      param.preprocess.file.midfix, param.preprocess.file.suffix);
  end
  
  % Assumption is that fns is in chronological order. Most radar systems
  % have filenames that are in chronological order with a simple
  % alphabetical sort.
  
  % ACORDS filenames are not in chronological order, resort by their
  % extension (.1, .2, .3, ..., .100, etc.)
  if any(param.preprocess.file.version == [405 406])
    basenames = {};
    file_idxs = [];
    new_fns = {};
    finfo_param.hnum = 1;
    finfo_param.preprocess.file.version = param.preprocess.file.version;
    for fidx = 1:length(fns)
      fname = fname_info_acords(fns{fidx},finfo_param);
      new_fns{fidx} = [fname.basename sprintf('.%03d',fname.file_idx)];
    end
    [new_fns,sorted_idxs] = sort(new_fns);
    fns = fns(sorted_idxs);
  end
  
  %% Read Headers: Header Info
  hdr_param = struct('file_mode','ieee-be');
  if any(param.preprocess.file.version == [1])
    hdr_param.frame_sync = uint32(hex2dec('DEADBEEF'));
    hdr_param.field_offsets = int32([4 16 20]); % epri seconds fractions
    hdr_param.field_types = {uint32(1) uint32(1) uint32(1)};
    
  elseif any(param.preprocess.file.version == [2])
    hdr_param.frame_sync = uint32(hex2dec('BADA55E5'));
    hdr_param.field_offsets = int32([4 8 12 24]); % epri seconds fractions loopback/nyquist-zone
    hdr_param.field_types = {uint32(1) uint32(1) uint32(1) uint32(1)};
    
  elseif any(param.preprocess.file.version == [4])
    hdr_param.frame_sync = uint32(hex2dec('BADA55E5'));
    hdr_param.field_offsets = int32([4 8 12 16]); % epri sec1 sec2 fractions
    hdr_param.field_types = {uint32(1) uint32(1) uint32(1) uint64(1)};
    
  elseif any(param.preprocess.file.version == [3 5])
    hdr_param.frame_sync = uint32(hex2dec('BADA55E5'));
    hdr_param.field_offsets = int32(4*[1 2 3 9 10 11]); % epri seconds fractions start/stop-index DDCfield1 DDCfield2
    hdr_param.field_types = {uint32(1) uint32(1) uint32(1) uint32(1) uint32(1) uint32(1)};
    
  elseif any(param.preprocess.file.version == [8])
    hdr_param.frame_sync = uint32(hex2dec('BADA55E5'));
    hdr_param.field_offsets = int32([4 8 12 16 33 36 38 40]);
    % epri seconds fractions counter nyquist-zone waveform-ID
    hdr_param.field_types = {uint32(1) uint32(1) uint32(1) uint64(1) uint8(1) uint16(1) uint16(1) uint64(1)};
    
  elseif any(param.preprocess.file.version == [101])
    hdr_param.frame_sync = uint32(hex2dec('DEADBEEF'));
    hdr_param.field_offsets = int32([4 8 12]); % epri seconds fractions
    hdr_param.field_types = {uint32(1) uint32(1) uint32(1)};
    
  elseif any(param.preprocess.file.version == [102])
    hdr_param.frame_sync = uint32(hex2dec('1ACFFC1D'));
    hdr_param.field_offsets = int32(4*[1 3 4 5 6]); % epri seconds fractions
    hdr_param.field_types = {uint32(1) uint32(1) uint32(1) uint32(1) uint32(1)};
    hdr_param.frame_sync = hex2dec('1ACFFC1D');
    
  elseif any(param.preprocess.file.version == [401])
    hdr_param.frame_sync = uint32(hex2dec('DEADBEEF'));
    hdr_param.field_offsets = int32([16 8 12]); % epri seconds fractions
    hdr_param.field_types = {uint32(1) uint32(1) uint32(1)};
    
  elseif any(param.preprocess.file.version == [402 403])
    hdr_param.frame_sync = uint32(hex2dec('BADA55E5'));
    hdr_param.field_offsets = int32([4 8 12 16]); % epri seconds fraction counter
    hdr_param.field_types = {uint32(1) uint32(1) uint32(1) uint64(1)};
    
  elseif any(param.preprocess.file.version == [404])
    hdr_param.frame_sync = uint32(hex2dec('1ACFFC1D'));
    hdr_param.field_offsets = int32([4 16 20]); % epri seconds fractions
    hdr_param.field_types = {uint32(1) uint32(1) uint32(1)};
    
  elseif any(param.preprocess.file.version == [405 406])
    hdr_param.file_mode = 'ieee-le';
    hdr_param.frame_sync = uint32(0);
    hdr_param.field_offsets = int32([0 4]); % epri seconds fractions
    hdr_param.field_types = {uint32(1) uint32(1)};
    
  elseif any(param.preprocess.file.version == [7 407 408])
    hdr_param.frame_sync = uint32(hex2dec('1ACFFC1D'));
    if hdr.file_version == 7
      hdr_param.field_offsets = int32([4 8 12 16]); % epri seconds fractions counter
    end
    
  else
    error('Unsupported radar %s', param.radar_name);
  end
  
  %% Read Headers: File Loop
  failed_load{board_idx} = false(size(fns));
  for fn_idx = 1:length(fns)
    
    % Create temporary filename that will store the header information for
    % this file.
    fn = fns{fn_idx};
    if any(param.preprocess.file.version == [405 406])
      [~,fn_name,ext] = fileparts(fn);
      fn_name = [fn_name,ext];
    else
      [~,fn_name] = fileparts(fn);
    end
    tmp_hdr_fn = ct_filename_ct_tmp(param,'','headers', ...
      fullfile(board_folder_name, [fn_name '.mat']));
    tmp_hdr_fn_dir = fileparts(tmp_hdr_fn);
    if ~exist(tmp_hdr_fn_dir,'dir')
      mkdir(tmp_hdr_fn_dir);
    end
    
    if param.preprocess.online_mode == 0
      % Reading in all files one time, print each out
      fprintf('%d of %d: %s (%s)\n', fn_idx, length(fns), fn, datestr(now,'HH:MM:SS'));
      if param.preprocess.reuse_tmp_files && exist(tmp_hdr_fn,'file')
        continue;
      end
    else
      % Repeatedly reading in all files as new files are added, only print
      % out the filename if it has not been processed yet
      if param.preprocess.reuse_tmp_files && exist(tmp_hdr_fn,'file')
        continue;
      else
        fprintf('%d of %d: %s (%s)\n', fn_idx, length(fns), fn, datestr(now,'HH:MM:SS'));
      end
    end
    
    try
      if any(param.preprocess.file.version == [1])
        hdr = basic_load_fmcw(fn);
        wfs = hdr.wfs;
      elseif any(param.preprocess.file.version == [2])
        hdr = basic_load_fmcw2(fn, struct('file_version',param.preprocess.file.version));
        wfs = hdr.wfs;
      elseif any(param.preprocess.file.version == [4])
        hdr = basic_load_fmcw2(fn, struct('file_version',param.preprocess.file.version,'clk',param.preprocess.daq.clk));
        wfs = hdr.wfs;
      elseif any(param.preprocess.file.version == [3 5])
        hdr = basic_load_fmcw3(fn, struct('file_version',param.preprocess.file.version));
        wfs = struct('presums',hdr.presums);
      elseif any(param.preprocess.file.version == [6])
        hdr = basic_load_fmcw4(fn, struct('file_version',param.preprocess.file.version));
        wfs = struct('presums',hdr.presums);
      elseif any(param.preprocess.file.version == [7])
        hdr = basic_load(fn);
        wfs = hdr.wfs;
      elseif any(param.preprocess.file.version == [8])
        hdr = basic_load_fmcw8(fn, struct('file_version',param.preprocess.file.version));
        wfs = struct('presums',hdr.presums);
      elseif any(param.preprocess.file.version == [101])
        hdr = basic_load_accum(fn);
        wfs = hdr.wfs;
      elseif any(param.preprocess.file.version == [102])
        hdr = basic_load_accum2(fn);
        wfs = hdr.wfs;
      elseif any(param.preprocess.file.version == [401])
        hdr = basic_load_mcords(fn);
        wfs = hdr.wfs;
      elseif any(param.preprocess.file.version == [402])
        hdr = basic_load_mcords2(fn);
        wfs = hdr.wfs;
      elseif any(param.preprocess.file.version == [403])
        hdr = basic_load_mcords3(fn);
        wfs = hdr.wfs;
      elseif any(param.preprocess.file.version == [404])
        hdr = basic_load_mcords4(fn);
        wfs = hdr.wfs;
      elseif any(param.preprocess.file.version == [405 406])
        % Load header information that never changes
        %   You need to get the record sizes
        clear hdr wfs
        hdr = basic_load_acords(fn,struct('datatype',0,'file_version',param.preprocess.file.version,'verbose',0));
        for hidx = 1:length(hdr)
          wfs{1}.num_samp(hidx) = hdr(hidx).num_samp;
          wfs{2}.num_samp(hidx) = hdr(hidx).num_samp;
          wfs{1}.presums(hidx) = hdr(hidx).presums;
          wfs{2}.presums(hidx) = hdr(hidx).presums;
          wfs{1}.shifts(hidx) = hdr(hidx).shifts;
          wfs{2}.shifts(hidx) = hdr(hidx).shifts;
          wfs{1}.prf(hidx) = hdr(hidx).prf;
          wfs{2}.prf(hidx) = hdr(hidx).prf;
          wfs{1}.f0(hidx) = hdr(hidx).f0;
          wfs{2}.f0(hidx) = hdr(hidx).f0;
          wfs{1}.f1(hidx) = hdr(hidx).f1;
          wfs{2}.f1(hidx) = hdr(hidx).f1;
          wfs{1}.wf_gen_clk(hidx) = hdr(hidx).wf_gen_clk;
          wfs{2}.wf_gen_clk(hidx) = hdr(hidx).wf_gen_clk;
          wfs{1}.daq_clk(hidx) = hdr(hidx).daq_clk;
          wfs{2}.daq_clk(hidx) = hdr(hidx).daq_clk;
          wfs{1}.tpd(hidx) = hdr(hidx).tpd;
          wfs{2}.tpd(hidx) = hdr(hidx).tpd;
          wfs{1}.tx_win(hidx) = hdr(hidx).tx_win;
          wfs{2}.tx_win(hidx) = hdr(hidx).tx_win;
          wfs{1}.t0(hidx) = hdr(hidx).samp_win_delay;
          wfs{2}.t0(hidx) = hdr(hidx).samp_win_delay;
          if param.preprocess.file.version == 406
            wfs{1}.elem_slots(hidx,:) = [hdr(hidx).elem_1 hdr(hidx).elem_2 hdr(hidx).elem_3 hdr(hidx).elem_4];
            wfs{2}.elem_slots(hidx,:) = [hdr(hidx).elem_1 hdr(hidx).elem_2 hdr(hidx).elem_3 hdr(hidx).elem_4];
            wfs{1}.blank(hidx) = hdr(hidx).rx_blank_end;
            wfs{1}.adc_gains(hidx,:) = 10.^((44-hdr(hidx).low_gain_atten).*ones(1,hdr(hidx).num_elem+1)./20);
            wfs{2}.blank(hidx) = hdr(hidx).hg_blank_end;
            wfs{2}.adc_gains(hidx,:) = 10.^((80-hdr(hidx).high_gain_atten).*ones(1,hdr(hidx).num_elem+1)./20);
          elseif param.preprocess.file.version == 405
            wfs{1}.blank(hidx) = hdr(hidx).rx_blank_end;
            wfs{1}.adc_gains(hidx,:) = 10.^(44-hdr(hidx).low_gain_atten./20);
            wfs{2}.blank(hidx) = hdr(hidx).hg_blank_end;
            wfs{2}.adc_gains(hidx,:) = 10.^(80-hdr(hidx).high_gain_atten./20);
          end
        end
      elseif any(param.preprocess.file.version == [407 408])
        try
          hdr = basic_load_mcords5(fn,struct('presum_bug_fixed',presum_bug_fixed));
          hdr_param.frame_sync = uint32(hex2dec('1ACFFC1D'));
          if hdr.file_version == 407
            hdr_param.field_offsets = int32([4 16 20 24]); % epri seconds fractions counter
          elseif hdr.file_version == 408
            hdr_param.field_offsets = int32([4 32 36 48]); % epri seconds fractions counter
          end
          hdr_param.field_types = {uint32(1) uint32(1) uint32(1) uint64(1)};
        catch ME
          if 0
            error(ME);
          else
            fprintf('Warning HACK enabled for mcords5 without frame sync field\n');
            fn_hack = '/mnt/HDD10/1805101801/UWB/chan6/mcords5_06_20180510_112936_00_0000.bin';
            hdr = basic_load_mcords5(fn_hack,struct('presum_bug_fixed',presum_bug_fixed));
            hdr_param.frame_sync = uint32(hex2dec('01600558')); % Used for 20180510 Greenland Polar6 recovery
            hdr_param.field_offsets = int32([4 16 20 24]-36); % epri seconds fractions counter % Used for 20180511 Greenland Polar6 recovery
            hdr_param.field_types = {uint32(1) uint32(1) uint32(1) uint64(1)};
          end
        end
        wfs = hdr.wfs;
        for wf=1:length(wfs); wfs(wf).file_version = hdr.file_version; end;
      end
    catch ME
      ME
      warning('  Failed to load... skipping.\n');
      failed_load{board_idx}(fn_idx) = true;
      continue;
    end
    
    if any(param.preprocess.file.version == [1])
      [file_size offset epri seconds fraction] = basic_load_hdr_mex(fn,hdr_param.frame_sync,hdr_param.field_offsets,hdr_param.field_types,hdr_param.file_mode);
      
      % Find bad records by checking their size (i.e. the distance between
      % frame syncs which should be constant).
      expected_rec_size = median(diff(offset));
      meas_rec_size = diff(offset);
      bad_mask = meas_rec_size ~= expected_rec_size;
      bad_mask(end+1) = file_size < offset(end) + expected_rec_size;
      
      % Remove bad records (i.e. ones with sizes that are not expected
      offset = double(offset(~bad_mask));
      epri = double(epri(~bad_mask));
      seconds = double(seconds(~bad_mask));
      fraction = double(fraction(~bad_mask));
      
      save(tmp_hdr_fn,'offset','epri','seconds','fraction','wfs');
      
    elseif any(param.preprocess.file.version == [2])
      [file_size offset epri seconds fraction tmp] = basic_load_hdr_mex(fn,hdr_param.frame_sync,hdr_param.field_offsets,hdr_param.field_types,hdr_param.file_mode);
      loopback = mod(floor(tmp/2^16),2^2) - 1;
      nyquist_zone = mod(tmp,2^3) - 1;
      
      % Find bad records by checking their size (i.e. the distance between
      % frame syncs which should be constant).
      expected_rec_size = median(diff(offset));
      meas_rec_size = diff(offset);
      bad_mask = meas_rec_size ~= expected_rec_size;
      bad_mask(end+1) = file_size < offset(end) + expected_rec_size;
      
      % Remove bad records (i.e. ones with sizes that are not expected
      offset = double(offset(~bad_mask));
      epri = double(epri(~bad_mask));
      seconds = double(seconds(~bad_mask));
      fraction = double(fraction(~bad_mask));
      loopback = double(loopback(~bad_mask));
      nyquist_zone = double(nyquist_zone(~bad_mask));
      
      save(tmp_hdr_fn,'offset','epri','seconds','fraction','loopback','nyquist_zone','wfs');
      
    elseif any(param.preprocess.file.version == [4])
      [file_size offset epri sec1 sec2 fraction] = basic_load_hdr_mex(fn,hdr_param.frame_sync,hdr_param.field_offsets,hdr_param.field_types,hdr_param.file_mode);
      
      % Convert seconds from NMEA ASCII string
      %   64 bits: 00 HH MM SS
      %   ASCII zero is "48"
      seconds = ((floor(sec1/2^8)-48)*10 + mod(sec1,2^8)-48) * 3600 ...
        + ((floor(sec2/2^24)-48)*10 + mod(floor(sec2/2^16),2^8)-48) * 60 ...
        + ((floor(mod(sec2,2^16)/2^8)-48)*10 + mod(mod(sec2,2^16),2^8)-48);
      
      % Find bad records by checking their size (i.e. the distance between
      % frame syncs which should be constant).
      expected_rec_size = median(diff(offset));
      meas_rec_size = diff(offset);
      bad_mask = meas_rec_size ~= expected_rec_size;
      bad_mask(end+1) = file_size < offset(end) + expected_rec_size;
      
      % Remove bad records (i.e. ones with sizes that are not expected
      offset = double(offset(~bad_mask));
      epri = double(epri(~bad_mask));
      seconds = double(seconds(~bad_mask));
      fraction = double(fraction(~bad_mask));
      
      save(tmp_hdr_fn,'offset','epri','seconds','fraction','wfs');
      
    elseif any(param.preprocess.file.version == [3 5 6])
      [file_size offset epri seconds fraction hdr9 hdr10 hdr11] = basic_load_hdr_mex(fn,hdr_param.frame_sync,hdr_param.field_offsets,hdr_param.field_types,hdr_param.file_mode);
      
      start_idx = floor(hdr9/2^16);
      stop_idx = mod(hdr9,2^16);
      NCO_freq_step = mod(hdr10,2^16);
      if param.preprocess.file.version == 6
        nadir_or_sidelooking_select = mod(floor(hdr11/2^24),2^8);
      else
        nyquist_zone = mod(floor(hdr11/2^24),2^8);
      end
      if param.preprocess.file.version == 3
        DDC_filter_select = mod(floor(hdr11/2^16),2^8) + 1;
      else
        DDC_filter_select = mod(floor(hdr11/2^16),2^8);
      end
      if param.preprocess.file.version == 3 || param.preprocess.file.version == 6
        DDC_or_raw_select = mod(hdr11,2^8);
      else
        DDC_or_raw_select = mod(hdr11,2^8);
        DDC_filter_select(DDC_or_raw_select == 1) = DDC_filter_select(DDC_or_raw_select == 1) - 1;
        DDC_or_raw_select(DDC_or_raw_select == 1) = 0;
      end
      if DDC_or_raw_select
        % Raw data
        num_sam = stop_idx - start_idx;
      else
        % DDC data
        num_sam = floor(((stop_idx - start_idx) ./ 2.^(1+DDC_filter_select)));
      end
      HEADER_SIZE = 48;
      SAMPLE_SIZE = 2;
      expected_rec_size = HEADER_SIZE + SAMPLE_SIZE*double(num_sam).*(1+(DDC_or_raw_select==0));
      meas_rec_size = diff(offset);
      bad_mask = meas_rec_size ~= expected_rec_size(1:end-1);
      bad_mask(end+1) = file_size < offset(end) + expected_rec_size(end);
      if any(bad_mask)
        warning('Found %d of %d record size errors', sum(bad_mask), length(bad_mask));
      end
      offset = offset(~bad_mask);
      epri = epri(~bad_mask);
      seconds = seconds(~bad_mask);
      fraction = fraction(~bad_mask);
      
      seconds = BCD_to_seconds(seconds);
      save(tmp_hdr_fn,'offset','epri','seconds','fraction','wfs', ...
        'start_idx','stop_idx','DDC_filter_select','DDC_or_raw_select', ...
        'num_sam','nyquist_zone','NCO_freq_step');
      
    elseif any(param.preprocess.file.version == [8])
      [file_size offset epri seconds fraction counter nyquist_zone start_idx stop_idx waveform_ID] = basic_load_hdr_mex(fn,hdr_param.frame_sync,hdr_param.field_offsets,hdr_param.field_types,hdr_param.file_mode);
      
      HEADER_SIZE = 48;
      SAMPLE_SIZE = 2;
      num_sam = 2*(stop_idx - start_idx);
      expected_rec_size = HEADER_SIZE + SAMPLE_SIZE*double(num_sam);
      meas_rec_size = diff(offset);
      bad_mask = meas_rec_size ~= expected_rec_size(1:end-1);
      % Last record in file is probably bad if its size does not match any
      % of the good records.
      bad_mask(end+1) = all(expected_rec_size(end) ~= expected_rec_size(~bad_mask));
      if sum(bad_mask) > 0
        warning('Found %d of %d record size errors', sum(bad_mask), length(bad_mask));
      end
      offset = offset(~bad_mask);
      epri = epri(~bad_mask);
      seconds = seconds(~bad_mask);
      fraction = fraction(~bad_mask);
      counter = counter(~bad_mask);
      nyquist_zone = nyquist_zone(~bad_mask);
      start_idx = start_idx(~bad_mask);
      stop_idx = stop_idx(~bad_mask);
      
      seconds = BCD_to_seconds(seconds);
      save(tmp_hdr_fn,'offset','epri','seconds','fraction','wfs', ...
        'counter','nyquist_zone','start_idx','stop_idx');
      
    elseif any(param.preprocess.file.version == [101])
      [file_size offset unknown seconds fraction] = basic_load_hdr_mex(fn,hdr_param.frame_sync,hdr_param.field_offsets,hdr_param.field_types,hdr_param.file_mode);
      
      % Find bad records by checking their size (i.e. the distance between
      % frame syncs which should be constant).
      expected_rec_size = median(diff(offset));
      meas_rec_size = diff(offset);
      bad_mask = meas_rec_size ~= expected_rec_size;
      bad_mask(end+1) = file_size < offset(end) + expected_rec_size;
      
      % Remove bad records (i.e. ones with sizes that are not expected
      offset = double(offset(~bad_mask));
      unknown = double(unknown(~bad_mask));
      seconds = double(seconds(~bad_mask));
      fraction = double(fraction(~bad_mask));
      
      save(tmp_hdr_fn,'offset','unknown','seconds','fraction');
      
    elseif any(param.preprocess.file.version == [102])
      [file_size offset radar_time_ms radar_time_ls radar_time_1pps_ms radar_time_1pps_ls] ...
        = basic_load_hdr_mex(fn,hdr_param.frame_sync,hdr_param.field_offsets,hdr_param.field_types,hdr_param.file_mode);
      radar_time = (hdr_data(3,:)*2^32 + hdr_data(4,:)) / (param.preprocess.daq.clk/100);
      radar_time_1pps = (hdr_data(5,:)*2^32 + hdr_data(6,:)) / (param.preprocess.daq.clk/100);
      
      save(tmp_hdr_fn,'offset','radar_time','radar_time_1pps','wfs');
      
    elseif any(param.preprocess.file.version == [401])
      error('Not supported');
      
    elseif any(param.preprocess.file.version == [402])
      error('Not supported');
      
    elseif any(param.preprocess.file.version == [403])
      hdr_param.field_offsets = int32([4 8 12 16]); % epri seconds fraction counter
      hdr_param.field_types = {uint32(1) uint32(1) uint32(1) uint64(1)};
      [file_size offset epri seconds fraction counter] = basic_load_hdr_mex(fn,hdr_param.frame_sync,hdr_param.field_offsets,hdr_param.field_types,hdr_param.file_mode);
      seconds = BCD_to_seconds(seconds);
      
      % Find bad records by checking their size
      % The distance between frame syncs should be constant
      expected_rec_size = median(diff(offset));
      meas_rec_size = diff(offset);
      bad_mask = all(bsxfun(@(x,y) x ~= y, meas_rec_size, expected_rec_size(:)),1);
      % Note that we always assume that the last record in the file is
      % good (since it is a partial record and we would have to look at
      % the next file to see if the complete record is there)
      bad_mask(end+1) = false;
      
      % Remove bad records (i.e. ones with sizes that are not expected
      offset = double(offset(~bad_mask));
      epri = double(epri(~bad_mask));
      seconds = double(seconds(~bad_mask));
      fraction = double(fraction(~bad_mask));
      counter = double(counter(~bad_mask));
      
      save(tmp_hdr_fn,'offset','epri','seconds','fraction','counter','wfs');
      
    elseif any(param.preprocess.file.version == [404])
      
    elseif any(param.preprocess.file.version == [405 406])
      % Load header information that can change on every record AND
      % is required for records generation (radar time)
      %     radar_time =
      file_size = 0;
      offset = [];
      seconds = [];
      % Get header timestamps and offsets
      [hdr htime hoffset] = basic_load_acords(fn,struct('datatype',0,'file_version',param.preprocess.file.version,'verbose',0));
      raw_file_time = htime(1);
      % Get data records timestamps and offsets
      [data seconds offset] = basic_load_acords(fn,struct('datatype',2,'file_version',param.preprocess.file.version,'verbose',0));
      fractions = zeros(size(seconds));
      save(tmp_hdr_fn,'offset','seconds','hdr','hoffset','htime','wfs','raw_file_time');
      
    elseif any(param.preprocess.file.version == [7 407 408])
      hdr_param.field_types = {uint32(1) uint32(1) uint32(1) uint64(1)};
      [file_size offset epri seconds fraction counter] = basic_load_hdr_mex(fn,hdr_param.frame_sync,hdr_param.field_offsets,hdr_param.field_types,hdr_param.file_mode);
      seconds = BCD_to_seconds(seconds);
      
      % Find bad records by checking their size
      if any(param.preprocess.file.version == [407 408])
        % The distance between frame syncs should be constant
        expected_rec_size = median(diff(offset));
        meas_rec_size = diff(offset);
        bad_mask = all(bsxfun(@(x,y) x ~= y, meas_rec_size, expected_rec_size(:)),1);
        % Note that we always assume that the last record in the file is
        % good (since it is a partial record and we would have to look at
        % the next file to see if the complete record is there)
        bad_mask(end+1) = false;
        
      elseif any(param.preprocess.file.version == [7])
        % User must supply the valid record sizes
        if ~exist('expected_rec_sizes','var')
          fprintf('Record sizes found in this file:\n');
          fprintf('  %d', unique(diff(offset)));
          fprintf('\n');
          error('For file version 7 expected record sizes must be supplied.');
        end
        expected_rec_size = expected_rec_sizes;
        meas_rec_size = diff(offset);
        bad_mask = all(bsxfun(@(x,y) x ~= y, meas_rec_size, expected_rec_size(:)),1);
        bad_mask(end+1) = 1;
      end
      
      % Remove bad records (i.e. ones with sizes that are not expected
      offset = double(offset(~bad_mask));
      epri = double(epri(~bad_mask));
      seconds = double(seconds(~bad_mask));
      fraction = double(fraction(~bad_mask));
      counter = double(counter(~bad_mask));
      
      save(tmp_hdr_fn,'offset','epri','seconds','fraction','counter','wfs');
      
    end
  end
end

%% List bad files
for board_idx = 1:numel(param.preprocess.daq.board_map)
  if any(failed_load{board_idx})
    warning('Some files failed to load, consider deleting these to avoid problems.');
    for fn_idx = find(failed_load{board_idx})
      fprintf('  %s\n', fns_list{fn_idx});
    end
  end
end

%% Load temp files
% (only do this for the first board)

epri = [];
seconds = [];
fraction = [];
counter = [];
radar_time = [];
radar_time_1pps = [];
file_idxs = [];
unknown = [];
hdr_log = [];
htime = [];
hdr_raw = [];
hoffset = 0;
offset = 0;
board_idx = 1;
board = param.preprocess.daq.board_map{board_idx};
board_folder_name = param.preprocess.board_folder_name;
board_folder_name = regexprep(board_folder_name,'%b',board);
for fn_idx = 1:length(fns_list{board_idx})
  % Skip files that failed to load
  if failed_load{1}(fn_idx)
    continue;
  end
  
  % Create temporary header filename
  fn = fns_list{1}{fn_idx};
  if any(param.preprocess.file.version == [405 406])
    [~,fn_name,ext] = fileparts(fn);
    fn_name = [fn_name,ext];
  else
    [~,fn_name] = fileparts(fn);
  end
  
  tmp_hdr_fn = ct_filename_ct_tmp(param,'','headers', ...
    fullfile(board_folder_name, [fn_name '.mat']));
  tmp_hdr_fn_dir = fileparts(tmp_hdr_fn);
  
  if any(param.preprocess.file.version == [101])
    hdr = load(tmp_hdr_fn);
    unknown = cat(2,unknown,reshape(hdr.unknown,[1 length(hdr.unknown)]));
    seconds = cat(2,seconds,reshape(hdr.seconds,[1 length(hdr.seconds)]));
    fraction = cat(2,fraction,reshape(hdr.fraction,[1 length(hdr.fraction)]));
  elseif any(param.preprocess.file.version == [102])
    hdr = load(tmp_hdr_fn);
    radar_time = cat(2,epri,hdr.radar_time);
    radar_time_1pps = cat(2,epri,hdr.radar_time_1pps);
  elseif any(param.preprocess.file.version == [405 406])
    hdr = load(tmp_hdr_fn);
    hdr_log = [hdr_log,hdr.hdr];
    hdr_raw = [hdr_raw fn_idx*ones(1,length(hdr.hdr))];
    htime = [htime hdr.htime];
    seconds = cat(2,seconds,reshape(hdr.seconds,[1 length(hdr.seconds)]));
    fraction = cat(2,0,reshape(fraction,[1 length(fraction)]));
  else
    hdr = load(tmp_hdr_fn);
    epri = cat(2,epri,reshape(hdr.epri,[1 length(hdr.epri)]));
    seconds = cat(2,seconds,reshape(hdr.seconds,[1 length(hdr.seconds)]));
    fraction = cat(2,fraction,reshape(hdr.fraction,[1 length(hdr.fraction)]));
    if any(param.preprocess.file.version == [403 407 408])
      counter = cat(2,counter,reshape(hdr.counter,[1 length(hdr.counter)]));
    end
  end
  file_idxs = cat(2,file_idxs,fn_idx*ones([1 length(hdr.offset)]));
end

if param.preprocess.online_mode
  epri_jumps = diff(double(epri));
  fprintf('List of up to 10 last EPRI jumps of >100 records:\n');
  bad_jumps = epri_jumps(abs(epri_jumps) > 100);
  fprintf('  %.0f jumps. Up to last 10: ', length(bad_jumps));
  fprintf(' %.0f', bad_jumps(max(1,end-9):end));
  fprintf(' record jumps\n');
  
  utc_time_sod = double(seconds) + double(fraction) / param.preprocess.daq.clk;
  fprintf('UTC time SOD jumps of >0.5 sec:\n');
  utc_time_sod_jumps = diff(utc_time_sod);
  bad_jumps = utc_time_sod_jumps(abs(utc_time_sod_jumps) > 0.5);
  fprintf('  %.0f jumps. Up to last 10: ', length(bad_jumps));
  fprintf(' %.1f', bad_jumps(max(1,end-9):end));
  fprintf(' second jumps\n');
  
  return;
end

if 0
  % Remove large jumps in time
  % - 2016_Greenland_P3 jump to large value and then back down again
  big_sec_jump_idxs = find(abs(diff(double(seconds)))>1e5);
  fraction_wrap_idxs = find(diff(double(fraction))<0);
  mid_jumps = find(~ismember(big_sec_jump_idxs,fraction_wrap_idxs));
  if ~isempty(mid_jumps)
    [tmp,tmp_idxs] = min(abs(fraction_wrap_idxs - big_sec_jump_idxs(mid_jumps)));
    big_sec_jump_idxs(mid_jumps) = fraction_wrap_idxs(tmp_idxs);
  end
  big_sec_jump_idxs = big_sec_jump_idxs(ismember(big_sec_jump_idxs,fraction_wrap_idxs));
  if ~isempty(big_sec_jump_idxs)
    warning('Header seconds jump more than 1e5 sec, correcting jumps');
    if 0
      figure(101);plot(fraction);
      max_fraction = max(fraction);
      for idx = 1:length(big_sec_jump_idxs)
        figure(101);hold on;plot([big_sec_jump_idxs(idx),big_sec_jump_idxs(idx)]+1,[0,1.2*max_fraction],'r--');
      end
    end
    for idx = 1:2:length(big_sec_jump_idxs)
      jump_start = big_sec_jump_idxs(idx)+1;
      wrap_idxs = [find(fraction_wrap_idxs == big_sec_jump_idxs(idx)):find(fraction_wrap_idxs == big_sec_jump_idxs(idx+1))];
      wrap_idxs(1) = [];
      for wrap_idx = wrap_idxs
        seconds(jump_start:fraction_wrap_idxs(wrap_idx)) = seconds(jump_start-1)+1;
        jump_start = fraction_wrap_idxs(wrap_idx) + 1;
      end
    end
  end
end

%% Correct time variable
if any(param.preprocess.file.version == [1 2 3 4 5 6 7 8 101 403 407 408])
  utc_time_sod = double(seconds) + double(fraction) / param.preprocess.daq.clk;
  
  if 0
    % Test sequences
    utc_time_sod = [0 1 2 3 10000 5 6 7 8 9 10 11 12 13 24 25 26 27 28 19 20 21 22 19 20 21 22]
    utc_time_sod = utc_time_sod + 0.0001*randn(size(utc_time_sod))
    epri = 100 + [1:23, 20:23]
    epri(15) = 5000;
  end
  
  % Estimate the pulse repetition interval, PRI
  PRI = median(diff(utc_time_sod));
  
  % Create an EPRI sequence from the time record
  time_epri = utc_time_sod / PRI;
  [~,good_time_idx] = min(abs(utc_time_sod - median(utc_time_sod)));
  time_epri = time_epri - time_epri(good_time_idx);
  
  % Find the difference of the time-generated epri and the recorded epri
  dtime_epri = diff(time_epri);
  depri = diff(epri);
  
  % Find good/bad differences. Mask values are:
  %  0: both differences are bad
  %  1: EPRI good
  %  2: Time-generated EPRI good
  %  3: EPRI and time-generated EPRI good
  dtime_epri_threshold = 0.1; % Allow for 10% PRI error
  mask = (depri == 1) + (2*(abs(dtime_epri-1) < dtime_epri_threshold));
  % If the EPRI's both indicate the same number of skipped records,
  % consider it a good difference.
  mask(mask ~= 3 & depri == round(dtime_epri)) = 3;
  
  % Fix differenced time-generated EPRIs using differenced EPRIs
  dtime_epri(mask==1) = depri(mask==1);
  % Fix differenced EPRIs using differenced time-generated EPRIs
  depri(mask==2) = round(dtime_epri(mask==2));
  
  % Find sequences of good records (where mask > 0) and deal with each
  % segment separately.
  good_out_mask = false(size(utc_time_sod));
  start_idx = find(mask ~= 0,1);
  while ~isempty(start_idx)
    stop_idx = start_idx-1 + find(mask(start_idx+1:end) == 0,1);
    if isempty(stop_idx)
      stop_idx = numel(mask);
    end
    
    % Find a median point in each segment and assume this value is good
    [~,good_time_idx] = min(abs(utc_time_sod(start_idx:stop_idx+1) - median(utc_time_sod(start_idx:stop_idx+1))));
    [~,good_epri_idx] = min(abs(epri(start_idx:stop_idx+1) - median(epri(start_idx:stop_idx+1))));
    
    % Reconstruct epri
    tmp = [0 cumsum(depri(start_idx:stop_idx))];
    tmp = tmp - tmp(good_epri_idx) + epri(start_idx-1+good_epri_idx);
    epri_new(start_idx:stop_idx+1) = tmp;
    
    % Reconstruct time from time-generated EPRIs
    tmp = [0 cumsum(dtime_epri(start_idx:stop_idx))*PRI];
    tmp = tmp - tmp(good_time_idx) + utc_time_sod(start_idx-1+good_time_idx);
    utc_time_sod_new(start_idx:stop_idx+1) = tmp;
    
    % Mark these records as good outputs
    good_out_mask(start_idx:stop_idx+1) = true;
    
    % Find the next sequence
    start_idx = stop_idx + find(mask(stop_idx+1:end) ~= 0,1);
  end
  
  h_fig = figure(1); clf(h_fig); h_axes = axes('parent',h_fig);
  plot(h_axes,utc_time_sod);
  hold(h_axes,'on');
  plot(h_axes,utc_time_sod_new,'r');
  hold(h_axes,'off');
  xlabel(h_axes,'Record number');
  ylabel(h_axes,'UTC Time SOD (sec)');
  legend(h_axes,'Original','Corrected','location','best');
  title(h_axes,sprintf('UTC time original and corrected should\nmatch except at outliers'),'fontsize',10);
  
  UTC_MAX_ERROR = 0.1;
  h_fig = figure(2); clf(h_fig); h_axes = axes('parent',h_fig);
  plot(h_axes,utc_time_sod - utc_time_sod_new);
  xlabel(h_axes,'Record number');
  ylabel(h_axes,'Time correction (sec)');
  ylim(h_axes,[-UTC_MAX_ERROR UTC_MAX_ERROR]);
  title(h_axes,sprintf('Time correction should be within limits\nexcept for a few outliers.'),'fontsize',10);
  
  h_fig = figure(3); clf(h_fig);
  h_axes = subplot(2,1,1);
  plot(h_axes,diff(epri),'.');
  ylabel(h_axes,'Diff EPRI');
  h_axes = subplot(2,1,2);
  plot(h_axes,diff(epri),'.');
  ylim(h_axes,[-3 5]);
  xlabel(h_axes,'Record number');
  ylabel(h_axes,'Diff EPRI');
  title(h_axes,'Should be 1 except occasional record drops.','fontsize',10);
  
  fprintf('\n\n');
  warning(sprintf('==> Ensure that corrected time in figure 1 is good since this is used to create segment boundaries. Set utc_time_sod_new to correct value if not. Run "dbcont" when done reviewing (and applying corrections if they were needed).\n\nFigure 2 shows the applied time correction and if the time correction is outside the y-limits except for a few outliers it may indicate that there is a problem.\nFigure 3 shows the EPRI counter difference (how much it changes between each record). Many large jumps may indicate a problem in recording or in the headers themselves.\n'));
  fn_fig = ct_filename_ct_tmp(param,'','headers', fullfile(param.preprocess.date_str,'create_segments_utc_time_sod.fig'));
  fprintf('Saving %s\n', fn_fig);
  saveas(1,fn_fig);
  fn_fig = ct_filename_ct_tmp(param,'','headers', fullfile(param.preprocess.date_str,'create_segments_utc_time_sod_error.fig'));
  fprintf('Saving %s\n', fn_fig);
  saveas(2,fn_fig);
  fn_fig = ct_filename_ct_tmp(param,'','headers', fullfile(param.preprocess.date_str,'create_segments_epri.fig'));
  fprintf('Saving %s\n', fn_fig);
  saveas(3,fn_fig);
  
  utc_time_sod = utc_time_sod_new;
  
  % Check for day wraps in the UTC time seconds of day
  day_wrap_idxs = find(diff(utc_time_sod) < -50000);
  day_wrap_offset = zeros(size(utc_time_sod));
  for day_wrap_idx = day_wrap_idxs
    day_wrap_offset(day_wrap_idx+1:end) = day_wrap_offset(day_wrap_idx+1:end) + 86400;
  end
  utc_time_sod = utc_time_sod + day_wrap_offset;
  
  % Look for time gaps (this is used later for segmentation)
  time_gaps = [1 find(abs(diff(utc_time_sod)) > param.preprocess.max_time_gap)];
  
  % Plot results
  figure(1); clf;
  plot(utc_time_sod);
  ylabel('UTC time seconds of day');
  xlabel('Record');
  grid on;
  hold on;
  plot([1 time_gaps], utc_time_sod([1 time_gaps]),'ro');
  hold off;
  legend('UTC Time','Start of each gap','location','best');
  
elseif any(param.preprocess.file.version == [405 406])
  utc_time_sod = seconds; % this is actually comp_time but doesn't need to
  % be converted to actual utc_time_sod since it's only looking at gaps in
  % the data
  
  day_wrap_idxs = find(diff(utc_time_sod) < -50000);
  day_wrap_offset = zeros(size(utc_time_sod));
  for day_wrap_idx = day_wrap_idxs
    day_wrap_offset(day_wrap_idx+1:end) = day_wrap_offset(day_wrap_idx+1:end) + 86400;
  end
  utc_time_sod = utc_time_sod + day_wrap_offset;
  
  time_gaps = find(abs(diff(utc_time_sod)) > MAX_TIME_GAP);
  
  reason = {};
  
  figure(1); clf;
  plot(utc_time_sod);
  ylabel('UTC time seconds of day');
  xlabel('Record');
  grid on;
  hold on;
  plot(time_gaps, utc_time_sod(time_gaps),'ro');
  hold off;
  
  names = fieldnames(hdr_log(1));
  if param.preprocess.file.version == 406
    change_fields = [2 5 6 7 8 9 10 11 12 17 18 19 20 21];
  elseif param.preprocess.file.version == 405
    change_fields = [2 5 6 7 8 9 10 11 12];
  end
  
  % Detect changes in header file information that would necessitate
  % creating a new segment.
  hdr_gaps = [];
  change_log = {};
  for n_idx = change_fields
    %     fprintf('Checking %s...\n',names{n_idx})
    hdr_changes = eval(sprintf('find(diff(squeeze(horzcat(hdr_log(:).%s))) ~= 0)',names{n_idx}));
    hdr_gaps = [hdr_gaps, hdr_changes];
    if ~isempty(hdr_changes)
      reason = [reason; repmat({names{n_idx}},eval(sprintf('length(find(diff(squeeze(horzcat(hdr_log(:).%s))) ~= 0))',names{n_idx})),1)];
      for cl_idx=1:length(hdr_changes)
        change_log{end+1} = eval(sprintf('[hdr_log(hdr_changes(cl_idx)).%s hdr_log(hdr_changes(cl_idx)+1).%s]',names{n_idx},names{n_idx}));
      end
    end
  end
  
  [hdr_gaps I J] = unique(hdr_gaps);
  comb_reason = cell(1,length(I));
  comb_change = cell(1,length(I));
  for uidx=1:length(J)
    comb_reason{J(uidx)} = [comb_reason{J(uidx)} reason{uidx}];
    comb_change{J(uidx)} = [comb_change{J(uidx)} change_log{uidx}];
  end
  
  % Convert from file index to file offset
  hdr_idxs = [];
  time_idxs = [];
  for gaps_idx = 1:length(hdr_gaps)
    hdr_idxs(gaps_idx) = find(utc_time_sod == htime(hdr_gaps(gaps_idx)+1),1,'last')-1;
  end
  
  if ~isempty(hdr_gaps)
    figure(1);
    hold on;
    plot(hdr_idxs, utc_time_sod(hdr_idxs),'go');
    hold off;
  end
  
  [time_gaps I J] = unique([time_gaps hdr_idxs]);
  if ~isempty(reason)
    extend_unique = [repmat({'time'},length(time_gaps),1); comb_reason.'];
    reason_unique = extend_unique(I);
    extend_change_log = [repmat({'N/A'},length(time_gaps),1); comb_change.'];
    change_log_unique = extend_change_log(I);
  else
    reason_unique = repmat({''},length(I),1);
    change_log_unique = repmat({''},length(I),1);
  end
else
  error('Not supported');
end

% if any(strcmpi(radar_name,{'accum'}))
%   %% Accum filter of bad records
%   % Repeated data often occurs in accum during a FIFO overflow. Remove
%   % these blocks.
%   bad_mask = zeros(size(utc_time-sod));
%   start_segment_time = utc_time_sod(1);
%   cur_time = utc_time_sod(1);
%   for idx = 2:utc_time_sod
%     diff_time = utc_time_sod(idx) - utc_time_sod(idx-1);
%     if diff_time > MAX_TIME_GAP
%       start_segment_time = utc_time_sod(idx);
%       cur_time = utc_time_sod(idx);
%     else
%       if utc_time_sod(idx) < start_segment_time
%       end
%       if utc_time_sod(idx) < cur_time
%         bad_mask(idx) = 1;
%       end
%     end
%   end
% end

%% Create Segments
% (only do this for the first channel)

if any(param.preprocess.file.version == [403 404 407 408])
  % NI XML settings files available, break segments based on settings files
  % and header information
  
  xml_version = param.preprocess.daq.xml_version;
  cresis_xml_mapping;
  
  settings_fn_dir = fullfile(param.preprocess.base_dir,param.preprocess.config_folder_name);
  fprintf('%s\n', settings_fn_dir);
  
  % Read XML files in this directory
  [settings,settings_enc] = read_ni_xml_directory(settings_fn_dir,xml_file_prefix,false);
  
  % Get the date information out of the filename
  board_idx = 1;
  fn_datenums = [];
  for data_fn_idx = 1:length(fns_list{board_idx})
    fname = fname_info_mcords2(fns_list{board_idx}{data_fn_idx});
    fn_datenums(end+1) = fname.datenum;
  end
  
  %% Print out settings from each XML file (and plot if enabled)
  oparams = {};
  for set_idx = 1:length(settings)
    % Print out settings
    [~,settings_fn_name] = fileparts(settings(set_idx).fn);
    fprintf('===================== Setting %d =================\n', set_idx);
    fprintf('%s: %d waveforms\n', settings_fn_name, length(settings(set_idx).(config_var).Waveforms));
    if isfield(settings(set_idx),'XML_File_Path')
      fprintf('  %s\n', settings(set_idx).XML_File_Path{1}.values{1});
    end
    fprintf('   PRF:'); fprintf(' %g', settings(set_idx).(config_var).(prf_var)); fprintf('\n');
    fprintf('   Amp:'); fprintf(' %g', settings(set_idx).(config_var).(ram_amp_var)); fprintf('\n');
    fprintf('   Tukey:'); fprintf(' %g', settings(set_idx).(config_var).RAM_Taper); fprintf('\n');
    Tpd = double(settings(set_idx).(config_var).Waveforms(1).Len_Mult)*settings(set_idx).(config_var).Base_Len;
    fprintf('   f0-f1:'); fprintf(' %g-%g MHz %g us', settings(set_idx).(config_var).Waveforms(1).Start_Freq(1)/1e6, ...
      settings(set_idx).(config_var).Waveforms(1).Stop_Freq(1)/1e6, Tpd*1e6); fprintf('\n');
    fprintf('   Tx Mask:'); fprintf(' %g', settings(set_idx).(config_var).Waveforms(1).TX_Mask); fprintf('\n');
    for wf = 1:length(settings(set_idx).(config_var).Waveforms)
      fprintf('    WF %d Atten:', wf); fprintf(' %g', settings(set_idx).(config_var).Waveforms(wf).Attenuator_2); fprintf('\n');
      fprintf('    WF %d Len:', wf); fprintf(' %.1f us', 1e6*settings(set_idx).(config_var).Base_Len*settings(set_idx).(config_var).Waveforms(wf).Len_Mult); fprintf('\n');
    end
    
    if set_idx < length(settings)
      settings(set_idx).file_matches = find(fn_datenums >= settings(set_idx).datenum & fn_datenums < settings(set_idx+1).datenum);
    else
      settings(set_idx).file_matches = find(fn_datenums >= settings(set_idx).datenum);
    end
    
    % Use file header results to remove bad files
    settings(set_idx).day_wrap_offset = 0;
    
    % Associate default parameters with each settings
    default = default_radar_params_settings_match(param.preprocess.defaults,settings(set_idx));
    oparams{end+1} = default;
    oparams{end} = rmfield(oparams{end},'config_regexp');
    oparams{end} = rmfield(oparams{end},'name');
    
    % Parameter spreadsheet
    % =======================================================================
    oparams{end}.day_seg = sprintf('%s_%02d',param.preprocess.date_str,length(oparams));
    oparams{end}.cmd.notes = default.name;
    
    oparams{end}.records.file.base_dir = param.preprocess.base_dir;
    oparams{end}.records.file.board_folder_name = param.preprocess.board_folder_name;
    oparams{end}.records.gps.time_offset = oparams{end}.records.gps.time_offset+settings(set_idx).day_wrap_offset;
    if ~isempty(oparams{end}.records.file.board_folder_name) ...
        && oparams{end}.records.file.board_folder_name(1) ~= filesep
      % Ensures that board_folder_name is not a text number which Excel
      % will misinterpret as a numeric type
      oparams{end}.records.file.board_folder_name = ['/' oparams{end}.records.file.board_folder_name];
    end
    oparams{end}.records.file.clk = param.preprocess.daq.clk;
    oparams{end}.radar.prf = settings(set_idx).(config_var).(prf_var);

    % Usually the default.radar.wfs structure only has one waveform
    % entry which is to be copied to all the waveforms so we keep "wf"
    % and "wf" separate.
    if numel(oparams{end}.radar.wfs) == 1
      oparams{end}.radar.wfs = repmat(oparams{end}.radar.wfs,[1 numel(settings(set_idx).(config_var).Waveforms)]);
    end
      
    for wf = 1:numel(settings(set_idx).(config_var).Waveforms)
      oparams{end}.radar.wfs(wf).Tpd = double(settings(set_idx).(config_var).Waveforms(wf).Len_Mult)*settings(set_idx).(config_var).Base_Len;
      oparams{end}.radar.wfs(wf).f0 = settings(set_idx).(config_var).Waveforms(wf).Start_Freq(1);
      oparams{end}.radar.wfs(wf).f1 = settings(set_idx).(config_var).Waveforms(wf).Stop_Freq(1);
      oparams{end}.radar.wfs(wf).tukey = settings(set_idx).(config_var).RAM_Taper;
      % Transmit weights
      if any(param.preprocess.file.version == [403 407 408])
        tx_mask_inv = fliplr(~(dec2bin(double(settings(set_idx).(config_var).Waveforms(wf).TX_Mask),8) - '0'));
        tx_weights = double(settings(set_idx).(config_var).(ram_var)) .* tx_mask_inv / param.preprocess.daq.max_wg_counts*param.preprocess.daq.max_wg_voltage;
      else
        tx_mask_inv = ~(dec2bin(double(settings(set_idx).(config_var).Waveforms(wf).TX_Mask),8) - '0');
        tx_weights = double(settings(set_idx).(config_var).(ram_var)) .* tx_mask_inv / param.preprocess.daq.max_wg_counts*param.preprocess.daq.max_wg_voltage;
      end
      tx_weights = tx_weights(logical(param.preprocess.daq.tx_mask));
      oparams{end}.radar.wfs(wf).tx_weights = tx_weights;
      
      % ADC Gains
      atten = double(settings(set_idx).(config_var).Waveforms(wf).Attenuator_1(1)) ...
        + double(settings(set_idx).(config_var).Waveforms(wf).Attenuator_2(1));
      oparams{end}.radar.wfs(wf).adc_gains = 10.^((param.preprocess.daq.rx_gain - atten(1)*ones(1,length(oparams{end}.radar.wfs(wf).rx_paths)))/20);
      
      % DDC mode and frequency
      oparams{end}.radar.wfs(wf).DDC_dec = 2^(2+settings(set_idx).DDC_Ctrl.DDC_sel.Val);
      oparams{end}.radar.wfs(wf).DDC_freq = settings(set_idx).DDC_Ctrl.(NCO_freq)*1e6;
    end
    
    % Adjust start/stop files for this segment if there are time gaps in
    % the start/stop files.
    start_idx = settings(set_idx).file_matches(1);
    while start_idx <= settings(set_idx).file_matches(end)
      mask = file_idxs==start_idx;
      if all(diff(utc_time_sod(mask)) <= param.preprocess.max_time_gap)
        % Found a good start file
        break;
      end
      start_idx = start_idx + 1;
    end
    
    stop_idx = settings(set_idx).file_matches(end);
    while stop_idx >= 1
      mask = file_idxs==stop_idx;
      if all(diff(utc_time_sod(mask)) <= param.preprocess.max_time_gap)
        % Found a good end file
        break;
      end
      stop_idx = stop_idx - 1;
    end
    
    settings(set_idx).file_matches = start_idx:stop_idx;
    
    mask = file_idxs >= start_idx & file_idxs <= stop_idx;
    
    time_gaps = diff(utc_time_sod(mask)) > param.preprocess.max_time_gap;
    if any(time_gaps)
      % Create multiple segments because there are time jumps in this
      % segment (probably due to recording errors) which exceed the allowed
      % amount.
      
      % NOT SUPPORT YET...
      keyboard
    end
    
    oparams{end}.records.file.start_idx = start_idx;
    oparams{end}.records.file.stop_idx = stop_idx;
  end
  
elseif any(param.preprocess.file.version == [410])
  % MCRDS headers available, break segments based on filenames
  
else
  % No heading information, break segments based on time or radar counter
  % information.
  % Using time
  bad_mask = logical(zeros(size(fns_list{1})));
  segments = [];
  segment_start = file_idxs(1);
  start_time = utc_time_sod(1);
  start_day_wrap_offset = day_wrap_offset(1);
  seg_idx = 0;
  for gap_idx = 1:length(time_gaps)
    time_gap = time_gaps(gap_idx);
    if file_idxs(time_gap) == file_idxs(time_gap + 1)
      bad_mask(file_idxs(time_gap)) = 1;
    end
    
    if bad_mask(file_idxs(time_gap))
      segment_stop = file_idxs(time_gap)-1;
    else
      segment_stop = file_idxs(time_gap);
    end
    
    if segment_stop - segment_start + 1 >= param.preprocess.min_seg_size
      seg_idx = seg_idx + 1;
      segments(seg_idx).start_time = start_time;
      segments(seg_idx).start_idx = segment_start;
      segments(seg_idx).stop_idx = segment_stop;
      segments(seg_idx).day_wrap_offset = start_day_wrap_offset;
      [~,fn_start_name,fn_start_name_ext] = fileparts(fns_list{1}{segment_start});
      [~,fn_stop_name,fn_stop_name_ext] = fileparts(fns_list{1}{segment_stop});
      fprintf('%2d: %s %4d-%4d %s - %s\n', seg_idx, ...
        datestr(epoch_to_datenum(start_time)), segment_start, segment_stop,...
        [fn_start_name fn_start_name_ext], [fn_stop_name fn_stop_name_ext]);
    end
    
    segment_start = file_idxs(time_gap)+1;
    start_time = utc_time_sod(time_gap+1);
    start_day_wrap_offset = day_wrap_offset(time_gap+1);
  end
  segment_stop = length(fns_list{1});
  if segment_stop - segment_start + 1 >= param.preprocess.min_seg_size
    seg_idx = seg_idx + 1;
    segments(seg_idx).start_time = start_time;
    segments(seg_idx).start_idx = segment_start;
    segments(seg_idx).stop_idx = segment_stop;
    segments(seg_idx).day_wrap_offset = start_day_wrap_offset;
    [~,fn_start_name,fn_start_name_ext] = fileparts(fns_list{1}{segment_start});
    [~,fn_stop_name,fn_stop_name_ext] = fileparts(fns_list{1}{segment_stop});
    fprintf('%2d: %s %4d-%4d %s - %s\n', seg_idx, ...
      datestr(epoch_to_datenum(start_time)), segment_start, segment_stop,...
      [fn_start_name fn_start_name_ext], [fn_stop_name fn_stop_name_ext]);
  end
  
  [~,sort_idxs] = sort(cell2mat({segments.start_time}));
  segments = segments(sort_idxs);
  
  oparams = {};
  for segment_idx = 1:length(segments)
    segment = segments(segment_idx);
    
    % Determine which default parameters to use
    % =======================================================================
    match_idx = 1;
    oparams{end+1} = param.preprocess.defaults{match_idx};
    oparams{end} = rmfield(oparams{end},'config_regexp');
    oparams{end} = rmfield(oparams{end},'name');
    
    % Parameter spreadsheet
    % =======================================================================
    oparams{end}.day_seg = sprintf('%s_%02d',param.preprocess.date_str,segment_idx);
    oparams{end}.cmd.notes = param.preprocess.defaults{match_idx}.name;
    
    for board_idx = 1:numel(param.preprocess.daq.board_map)
      oparams{end}.records.file.start_idx(board_idx) = segment.start_idx;
      oparams{end}.records.file.stop_idx(board_idx) = segment.stop_idx;
    end
    oparams{end}.records.file.base_dir = param.preprocess.base_dir;
    oparams{end}.records.file.board_folder_name = param.preprocess.board_folder_name;
    if ~isempty(oparams{end}.records.file.board_folder_name) ...
        && oparams{end}.records.file.board_folder_name(1) ~= filesep
      % Ensures that board_folder_name is not a text number which Excel
      % will misinterpret as a numeric type
      oparams{end}.records.file.board_folder_name = ['/' oparams{end}.records.file.board_folder_name];
    end
    oparams{end}.records.gps.time_offset = oparams{end}.records.gps.time_offset+segment.day_wrap_offset;
    if ~isnan(str2double(oparams{end}.records.file.board_folder_name))
      oparams{end}.records.file.board_folder_name = ['/' oparams{end}.records.file.board_folder_name];
    end
    oparams{end}.records.file.clk = param.preprocess.daq.clk;
    
    for wf_idx = 1:length(hdr.wfs)
      wf = hdr.wfs(wf_idx);
    end
  end
  
end

%% Print out segments
% =========================================================================
if ~isempty(param.preprocess.param_fn)
  % Print parameter spreadsheet values
  % =========================================================================
  fprintf('<strong>%s\n','='*ones(1,80)); fprintf('  cmd\n'); fprintf('%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.preprocess.param_fn,'cmd',oparams);
  fprintf('<strong>%s\n','='*ones(1,80)); fprintf('  records\n'); fprintf('%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.preprocess.param_fn,'records',oparams);
  fprintf('<strong>%s\n','='*ones(1,80)); fprintf('  qlook\n'); fprintf('%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.preprocess.param_fn,'qlook',oparams);
  fprintf('<strong>%s\n','='*ones(1,80)); fprintf('  sar\n'); fprintf('%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.preprocess.param_fn,'sar',oparams);
  fprintf('<strong>%s\n','='*ones(1,80)); fprintf('  array\n'); fprintf('%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.preprocess.param_fn,'array',oparams);
  fprintf('<strong>%s\n','='*ones(1,80)); fprintf('  radar\n'); fprintf('%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.preprocess.param_fn,'radar',oparams);
  fprintf('<strong>%s\n','='*ones(1,80)); fprintf('  post\n'); fprintf('%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.preprocess.param_fn,'post',oparams);
end

%% Exit task
% =========================================================================
fprintf('%s done %s\n', mfilename, datestr(now));

success = true;
