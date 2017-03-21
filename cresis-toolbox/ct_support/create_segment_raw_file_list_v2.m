% script create_segment_raw_file_list_v2.m
%
% Extracts headers from radar files and stores these in temporary
% files. These temporary files are used by create_records_*.m.
%
% Also helps find the file_prefix, start_idx, and stop_idx fields for
% creating segments in the parameter spreadsheet. This function will
% be typically run at the end of each flight day when the segments are being
% added to the spreadsheet for the first time.
%
% See run_create_segment_raw_file_list_v2.m for how to run.
%
% Author: John Paden

% =========================================================================
%% Automated Section
% =========================================================================
dbstack_info = dbstack;
fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', dbstack_info(1).name, day_string, datestr(now,'HH:MM:SS'));
fprintf('=====================================================================\n');

[output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name);

if isempty(file_prefix_override)
  file_prefix = radar_name;
else
  file_prefix = file_prefix_override;
end
if ~exist('file_regexp','var')
  file_regexp = '';
end
get_fns_param.regexp = file_regexp;

if ~exist('union_time_epri_gaps','var')
  union_time_epri_gaps = false;
end

if ~exist('online_mode','var')
  online_mode = false;
end

if ~isfield('param','gps_time_offset') || isempty(param.gps_time_offset)
  param.gps_time_offset = 1;
end

if ~exist('adcs','var')
  warning('adcs not set, assuming that adcs should be "1"');
  adcs = 1;
end

%% Load in the headers one ADC board at a time
failed_load = {};
for adc_idx = 1:length(adcs)
  %% Get the files for this ADC
  adc = adcs(adc_idx);
  board = adc_to_board(param.radar_name,adc);
  adc_folder_name = param.adc_folder_name;
  adc_folder_name = regexprep(adc_folder_name,'%02d',sprintf('%02.0f',adc));
  adc_folder_name = regexprep(adc_folder_name,'%d',sprintf('%.0f',adc));
  adc_folder_name = regexprep(adc_folder_name,'%b',sprintf('%.0f',board));
  
  fns = get_filenames(fullfile(base_dir,adc_folder_name), file_prefix, file_midfix, raw_file_suffix, get_fns_param);
  fns_list{adc_idx} = fns;

  if isempty(fns)
    error('No files found matching %s*%s*%s', ...
      fullfile(base_dir,adc_folder_name,file_prefix), file_midfix, ...
      raw_file_suffix);
  end
  
  % Sort ACORDS filenames because the extenions are not a standard length
  if any(strcmpi(radar_name,{'acords'}))
    basenames = {};
    file_idxs = [];
    new_fns = {};
    finfo_param.hnum = 1;
    finfo_param.file_version = param.file_version;
    for fidx = 1:length(fns)
      fname = fname_info_acords(fns{fidx},finfo_param);
      new_fns{fidx} = [fname.basename sprintf('.%03d',fname.file_idx)];
    end
    [new_fns,sorted_idxs] = sort(new_fns);
    fns = fns(sorted_idxs);
  end
  
  hdr_param = struct('file_mode','ieee-be');
  %% Setup the header information for this radar
  if any(strcmpi(radar_name,{'accum'}))
    hdr_param.frame_sync = uint32(hex2dec('DEADBEEF'));
    hdr_param.field_offsets = uint32([4 8 12]); % epri seconds fractions
    hdr_param.field_types = {uint32(1) uint32(1) uint32(1)};
    param.clk = 1e9/16;
  elseif any(strcmpi(radar_name,{'accum2'}))
    hdr_param.frame_sync = uint32(hex2dec('1ACFFC1D'));
    hdr_param.field_offsets = uint32(4*[1 3 4 5 6]); % epri seconds fractions
    hdr_param.field_types = {uint32(1) uint32(1) uint32(1) uint32(1) uint32(1)};
    hdr_param.field_offsets = [1 3 4 5 6];
    hdr_param.frame_sync = hex2dec('1ACFFC1D');
    param.clk = 1e9;
  elseif any(strcmpi(radar_name,{'acords'}))
    hdr_param.file_mode = 'ieee-le';
    hdr_param.frame_sync = uint32(0);
    hdr_param.field_offsets = uint32([0 4]); % epri seconds fractions
    hdr_param.field_types = {uint32(1) uint32(1)};
    param.clk = 55e6;
  elseif any(strcmpi(radar_name,{'mcords'}))
    hdr_param.frame_sync = uint32(hex2dec('DEADBEEF'));
    hdr_param.field_offsets = uint32([16 8 12]); % epri seconds fractions
    hdr_param.field_types = {uint32(1) uint32(1) uint32(1)};
    param.clk = 1e9/9;
  elseif any(strcmpi(radar_name,{'mcords2','mcords3'}))
    hdr_param.frame_sync = uint32(hex2dec('BADA55E5'));
    hdr_param.field_offsets = uint32([4 8 12]); % epri seconds fractions
    hdr_param.field_types = {uint32(1) uint32(1) uint32(1)};
    param.clk = 1e9/9;
  elseif any(strcmpi(radar_name,{'mcords4'}))
    hdr_param.frame_sync = uint32(hex2dec('1ACFFC1D'));
    hdr_param.field_offsets = uint32([4 16 20]); % epri seconds fractions
    hdr_param.field_types = {uint32(1) uint32(1) uint32(1)};
    param.clk = 1e9/2;
  elseif any(strcmpi(radar_name,{'mcords5','snow5'}))
    hdr_param.frame_sync = uint32(hex2dec('1ACFFC1D'));
  elseif any(strcmpi(radar_name,{'snow','kuband'}))
    hdr_param.frame_sync = uint32(hex2dec('DEADBEEF'));
    hdr_param.field_offsets = uint32([4 16 20]); % epri seconds fractions
    hdr_param.field_types = {uint32(1) uint32(1) uint32(1)};
    param.clk = 1e9/16;
  elseif any(strcmpi(radar_name,{'snow2','kuband2'}))
    if param.file_version == 2
      hdr_param.frame_sync = uint32(hex2dec('BADA55E5'));
      hdr_param.field_offsets = uint32([4 8 12 24]); % epri seconds fractions loopback/nyquist-zone
      hdr_param.field_types = {uint32(1) uint32(1) uint32(1) uint32(1)};
    elseif param.file_version == 4
      hdr_param.frame_sync = uint32(hex2dec('BADA55E5'));
      hdr_param.field_offsets = uint32([4 8 12 16]); % epri sec1 sec2 fractions
      hdr_param.field_types = {uint32(1) uint32(1) uint32(1) uint64(1)};
    else
      error('File version %d not supported for this radar %s.', param.file_version, param.radar_name);
    end
    param.clk = 1e9/8;
    
  elseif any(strcmpi(radar_name,{'snow3','kuband3','kaband3'}))
    hdr_param.frame_sync = uint32(hex2dec('BADA55E5'));
    hdr_param.field_offsets = uint32(4*[1 2 3 9 10 11]); % epri seconds fractions start/stop-index DDCfield1 DDCfield2
    hdr_param.field_types = {uint32(1) uint32(1) uint32(1) uint32(1) uint32(1) uint32(1)};
    param.clk = 1e9/8;
    
  elseif any(strcmpi(radar_name,{'snow8'}))
    hdr_param.frame_sync = uint32(hex2dec('BADA55E5'));
    hdr_param.field_offsets = uint32([4 8 12 16 33 36 38 40]);
    % epri seconds fractions counter nyquist-zone waveform-ID
    hdr_param.field_types = {uint32(1) uint32(1) uint32(1) uint64(1) uint8(1) uint16(1) uint16(1) uint64(1)};
    param.clk = 1e9/8;
    
  else
    error('Unsupported radar %s', param.radar_name);
  end
  
  %% Search for all the file prefixes
  failed_load{adc_idx} = zeros(size(fns));
  for fn_idx = 1:length(fns)
    fn = fns{fn_idx};
    if strcmp(radar_name,'acords')
      [~,fn_name,ext] = fileparts(fn);
      fn_name = [fn_name,ext];
    else
      [~,fn_name] = fileparts(fn);
    end
    
    if tmp_fn_uses_adc_folder_name
      tmp_hdr_fn = ct_filename_ct_tmp(param,'','headers', ...
        fullfile(adc_folder_name, [fn_name '.mat']));
    else
      tmp_hdr_fn = ct_filename_ct_tmp(param,'','headers', [fn_name '.mat']);
    end
    tmp_hdr_fn_dir = fileparts(tmp_hdr_fn);
    if ~exist(tmp_hdr_fn_dir,'dir')
      mkdir(tmp_hdr_fn_dir);
    end
    
    if ~online_mode
      fprintf('%d of %d: %s (%s)\n', fn_idx, length(fns), fn, datestr(now,'HH:MM:SS'));
      if reuse_tmp_files && exist(tmp_hdr_fn,'file')
        continue;
      end
    else
      if reuse_tmp_files && exist(tmp_hdr_fn,'file')
        continue;
      else
        fprintf('%d of %d: %s (%s)\n', fn_idx, length(fns), fn, datestr(now,'HH:MM:SS'));
      end
    end
    
    try
      if strcmp(radar_name,'accum')
        hdr = basic_load_accum(fn);
        wfs = hdr.wfs;
      elseif strcmp(radar_name,'accum2')
        hdr = basic_load_accum2(fn);
        wfs = hdr.wfs;
      elseif strcmp(radar_name,'acords')
        % Load header information that never changes
        %   You need to get the record sizes
        clear hdr wfs
        hdr = basic_load_acords(fn,struct('datatype',0,'file_version',param.file_version,'verbose',0));
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
          if param.file_version == 406
            wfs{1}.elem_slots(hidx,:) = [hdr(hidx).elem_1 hdr(hidx).elem_2 hdr(hidx).elem_3 hdr(hidx).elem_4];
            wfs{2}.elem_slots(hidx,:) = [hdr(hidx).elem_1 hdr(hidx).elem_2 hdr(hidx).elem_3 hdr(hidx).elem_4];
            wfs{1}.blank(hidx) = hdr(hidx).rx_blank_end;
            wfs{1}.adc_gains(hidx,:) = 10.^((44-hdr(hidx).low_gain_atten).*ones(1,hdr(hidx).num_elem+1)./20);
            wfs{2}.blank(hidx) = hdr(hidx).hg_blank_end;
            wfs{2}.adc_gains(hidx,:) = 10.^((80-hdr(hidx).high_gain_atten).*ones(1,hdr(hidx).num_elem+1)./20);
          elseif param.file_version == 405
            wfs{1}.blank(hidx) = hdr(hidx).rx_blank_end;
            wfs{1}.adc_gains(hidx,:) = 10.^(44-hdr(hidx).low_gain_atten./20);
            wfs{2}.blank(hidx) = hdr(hidx).hg_blank_end;
            wfs{2}.adc_gains(hidx,:) = 10.^(80-hdr(hidx).high_gain_atten./20);
          end
        end
      elseif strcmp(radar_name,'mcords')
        hdr = basic_load_mcords(fn);
        wfs = hdr.wfs;
      elseif strcmp(radar_name,'mcords2')
        hdr = basic_load_mcords2(fn);
        wfs = hdr.wfs;
      elseif strcmp(radar_name,'mcords3')
        hdr = basic_load_mcords3(fn);
        wfs = hdr.wfs;
      elseif strcmp(radar_name,'mcords4')
        hdr = basic_load_mcords4(fn);
        wfs = hdr.wfs;
      elseif strcmp(radar_name,'mcords5')
        hdr = basic_load_mcords5(fn,struct('presum_bug_fixed',presum_bug_fixed));
        wfs = hdr.wfs;
        for wf=1:length(wfs); wfs(wf).file_version = hdr.file_version; end;
      elseif any(strcmp(radar_name,{'snow','kuband'}))
        hdr = basic_load_fmcw(fn);
        wfs = hdr.wfs;
      elseif any(strcmp(radar_name,{'snow2','kuband2'}))
        hdr = basic_load_fmcw2(fn, struct('file_version',param.file_version));
        wfs = hdr.wfs;
      elseif any(strcmp(radar_name,{'snow3','kuband3','kaband3'}))
        if param.file_version == 6
          hdr = basic_load_fmcw4(fn, struct('file_version',param.file_version));
        else
          hdr = basic_load_fmcw3(fn, struct('file_version',param.file_version));
        end
        wfs = struct('presums',hdr.presums);
      elseif any(strcmp(radar_name,{'snow5'}))
        hdr = basic_load(fn);
        wfs = hdr.wfs;
      elseif any(strcmp(radar_name,{'snow8'}))
        hdr = basic_load_fmcw8(fn, struct('file_version',param.file_version));
        wfs = struct('presums',hdr.presums);
      end
    catch ME
      ME
      warning('  Failed to load... skipping.\n');
      failed_load{adc_idx}(fn_idx) = 1;
      continue;
    end
    
    if any(strcmpi(radar_name,{'accum'}))
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
      
    elseif any(strcmpi(radar_name,{'accum2'}))
      [file_size offset radar_time_ms radar_time_ls radar_time_1pps_ms radar_time_1pps_ls] ...
        = basic_load_hdr_mex(fn,hdr_param.frame_sync,hdr_param.field_offsets,hdr_param.field_types,hdr_param.file_mode);
      radar_time = (hdr_data(3,:)*2^32 + hdr_data(4,:)) / (param.clk/100);
      radar_time_1pps = (hdr_data(5,:)*2^32 + hdr_data(6,:)) / (param.clk/100);
      
      save(tmp_hdr_fn,'offset','radar_time','radar_time_1pps','wfs');
      
    elseif strcmp(radar_name,'acords')
      % Load header information that can change on every record AND
      % is required for records generation (radar time)
      %     radar_time =
      file_size = 0;
      offset = [];
      seconds = [];
      % Get header timestamps and offsets
      [hdr htime hoffset] = basic_load_acords(fn,struct('datatype',0,'file_version',param.file_version,'verbose',0));
      raw_file_time = htime(1);
      % Get data records timestamps and offsets
      [data seconds offset] = basic_load_acords(fn,struct('datatype',2,'file_version',param.file_version,'verbose',0));
      fractions = zeros(size(seconds));
      save(tmp_hdr_fn,'offset','seconds','hdr','hoffset','htime','wfs','raw_file_time');
      
    elseif strcmp(radar_name,'mcords')
      
    elseif strcmp(radar_name,'mcords2')
      
    elseif strcmp(radar_name,'mcords3')
      hdr_param.field_offsets = uint32([4 8 12 16]); % epri seconds fraction counter
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
      
    elseif strcmp(radar_name,'mcords4')
      
    elseif any(strcmp(radar_name,{'mcords5','snow5'}))
      if hdr.file_version == 407
        hdr_param.field_offsets = uint32([4 16 20 24]); % epri seconds fractions counter
      elseif hdr.file_version == 408
        hdr_param.field_offsets = uint32([4 32 36 48]); % epri seconds fractions counter
      elseif hdr.file_version == 7
        hdr_param.field_offsets = uint32([4 8 12 16]); % epri seconds fractions counter
      end
      hdr_param.field_types = {uint32(1) uint32(1) uint32(1) uint64(1)};
      [file_size offset epri seconds fraction counter] = basic_load_hdr_mex(fn,hdr_param.frame_sync,hdr_param.field_offsets,hdr_param.field_types,hdr_param.file_mode);
      seconds = BCD_to_seconds(seconds);
      
      % Find bad records by checking their size
      if hdr.file_version == 407 || hdr.file_version == 408
        % The distance between frame syncs should be constant
        expected_rec_size = median(diff(offset));
        meas_rec_size = diff(offset);
        bad_mask = all(bsxfun(@(x,y) x ~= y, meas_rec_size, expected_rec_size(:)),1);
        % Note that we always assume that the last record in the file is
        % good (since it is a partial record and we would have to look at
        % the next file to see if the complete record is there)
        bad_mask(end+1) = false;
      elseif hdr.file_version == 7
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
      
    elseif any(strcmp(radar_name,{'snow','kuband'}))
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
      
    elseif any(strcmp(radar_name,{'snow2','kuband2'}))
      if param.file_version == 2
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
        
      else
        [file_size offset epri sec1 sec2 fraction] = basic_load_hdr_mex(fn,hdr_param.frame_sync,hdr_param.field_offsets,hdr_param.field_types,hdr_param.file_mode);
        
        % Convert seconds from NMEA ASCII string
        %   64 bits: 00 HH MM SS
        %   ASCII zero is "48"
        seconds = ((floor(sec1/2^8)-48)*10 + mod(sec1,2^8)-48) * 3600 ...
          + ((floor(sec2/2^24)-48)*10 + mod(floor(sec2/2^16),2^8)-48) * 60 ...
          + ((floor(sec2/2^8)-48)*10 + mod(sec2,2^8)-48);
        
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
      end
      
    elseif any(strcmp(radar_name,{'snow3','kuband3','kaband3'}))
      [file_size offset epri seconds fraction hdr9 hdr10 hdr11] = basic_load_hdr_mex(fn,hdr_param.frame_sync,hdr_param.field_offsets,hdr_param.field_types,hdr_param.file_mode);
      
      start_idx = floor(hdr9/2^16);
      stop_idx = mod(hdr9,2^16);
      NCO_freq_step = mod(hdr10,2^16);
      if param.file_version == 6
        nadir_or_sidelooking_select = mod(floor(hdr11/2^24),2^8);
      else
        nyquist_zone = mod(floor(hdr11/2^24),2^8);
      end
      if param.file_version == 3
        DDC_filter_select = mod(floor(hdr11/2^16),2^8) + 1;
      else
        DDC_filter_select = mod(floor(hdr11/2^16),2^8);
      end
      if param.file_version == 3 | param.file_version == 6
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

      
    elseif any(strcmp(radar_name,{'snow8'}))
      [file_size offset epri seconds fraction counter nyquist_zone start_idx stop_idx waveform_ID] = basic_load_hdr_mex(fn,hdr_param.frame_sync,hdr_param.field_offsets,hdr_param.field_types,hdr_param.file_mode);
      
      HEADER_SIZE = 48;
      SAMPLE_SIZE = 2;
      num_sam = 2*(stop_idx - start_idx);
      expected_rec_size = HEADER_SIZE + SAMPLE_SIZE*double(num_sam);
      meas_rec_size = diff(offset);
      bad_mask = meas_rec_size ~= expected_rec_size(1:end-1);
      bad_mask(end+1) = file_size < offset(end) + expected_rec_size(end);
      if sum(bad_mask) > 1
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
    end
  end
end

%% Warn about any bad files
for adc_idx = 1:length(adcs)
  if any(failed_load{adc_idx})
    warning('Some files failed to load, consider deleting these to avoid problems.');
    for fn_idx = find(failed_load{adc_idx})
      fprintf('  %s\n', fns{fn_idx});
    end
  end
end

%% Load the parsed header data from temporary files
% (only do this for the first channel)
adc = adcs(1);
board = adc_to_board(param.radar_name,adc);
adc_folder_name = param.adc_folder_name;
adc_folder_name = regexprep(adc_folder_name,'%d',sprintf('%.0f',adc));
adc_folder_name = regexprep(adc_folder_name,'%b',sprintf('%.0f',board));

fns = get_filenames(fullfile(base_dir,adc_folder_name), file_prefix, file_midfix, raw_file_suffix, get_fns_param);

% Sort ACORDS filenames because the extenions are not a standard length
if any(strcmpi(radar_name,{'acords'}))
  basenames = {};
  file_idxs = [];
  new_fns = {};
  finfo_param.hnum = 1;
  finfo_param.file_version = param.file_version;
  for fidx = 1:length(fns)
    fname = fname_info_acords(fns{fidx},finfo_param);
    new_fns{fidx} = [fname.basename sprintf('.%03d',fname.file_idx)];
  end
  [new_fns,sorted_idxs] = sort(new_fns);
  fns = fns(sorted_idxs);
end

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
for fn_idx = 1:length(fns_list{1})
  if failed_load{1}(fn_idx)
    continue;
  end
  fn = fns_list{1}{fn_idx};
  if strcmp(radar_name,'acords')
    [~,fn_name,ext] = fileparts(fn);
    fn_name = [fn_name,ext];
  else
    [~,fn_name] = fileparts(fn);
  end
  
  if tmp_fn_uses_adc_folder_name
    tmp_hdr_fn = ct_filename_ct_tmp(param,'','headers', ...
      fullfile(adc_folder_name, [fn_name '.mat']));
  else
    tmp_hdr_fn = ct_filename_ct_tmp(param,'','headers', [fn_name '.mat']);
  end
  tmp_hdr_fn_dir = fileparts(tmp_hdr_fn);
  if ~exist(tmp_hdr_fn_dir,'dir')
    mkdir(tmp_hdr_fn_dir);
  end
  
  if any(strcmpi(radar_name,{'accum'}))
    hdr = load(tmp_hdr_fn);
    unknown = cat(2,unknown,reshape(hdr.unknown,[1 length(hdr.unknown)]));
    seconds = cat(2,seconds,reshape(hdr.seconds,[1 length(hdr.seconds)]));
    fraction = cat(2,fraction,reshape(hdr.fraction,[1 length(hdr.fraction)]));
  elseif any(strcmpi(radar_name,{'accum2'}))
    hdr = load(tmp_hdr_fn);
    radar_time = cat(2,epri,hdr.radar_time);
    radar_time_1pps = cat(2,epri,hdr.radar_time_1pps);
  elseif any(strcmpi(radar_name,{'acords'}))
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
    if any(strcmpi(radar_name,{'mcords3','mcords5'}))
      counter = cat(2,counter,reshape(hdr.counter,[1 length(hdr.counter)]));
    end
  end
  file_idxs = cat(2,file_idxs,fn_idx*ones([1 length(hdr.offset)]));
end

if online_mode
  epri_jumps = diff(double(epri));
  fprintf('List of EPRI jumps:\n');
  epri_jumps(abs(epri_jumps) > 100)
  
  utc_time_sod = double(seconds) + double(fraction) / param.clk;
  fprintf('List of UTC time SOD jumps:\n');
  utc_time_sod_jumps = diff(utc_time_sod);
  utc_time_sod_jumps(abs(utc_time_sod_jumps) > 0.5)
  
  return;
end

%% Correct and process time variable
%% Check to see if there are big soconds jumps (those encountered in 2016_Greenland_P3,
%% jump to a large number and then drop back to correct values)?
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
if any(strcmpi(radar_name,{'accum','snow','kuband','snow2','kuband2','snow3','kuband3','kaband3','snow5','mcords3','mcords5','snow8'}))
  utc_time_sod = double(seconds) + double(fraction) / param.clk;
  utc_time_sod = medfilt1(double(utc_time_sod));
  
  if counter_correction_en
    warning('You have enabled counter correction. Normally, this should not be necessary. Set correction parameters and then dbcont to continue');
    % This is an index into hdr.utc_time_sod that is correct
    counter_bad_threshold = 0.01;
    counter_min_freq = 100;
    counter_bin = 0.01;
    anchor_idx = 1;
    if any(strcmpi(radar_name,{'mcords5'}))
      counter_clk = param.clk;
      keyboard
    elseif any(strcmpi(radar_name,{'snow3','snow5','kuband3','accum'}))
      % counter_clk should be the EPRF (effective PRF after hardware presumming)
      % set anchor_idx to a record that you believe has the correct time
      counter_clk = 3906.250/2/8;
      counter = epri;
    elseif any(strcmpi(radar_name,{'mcords3'}))
      counter_clk = param.clk;
      keyboard
    elseif any(strcmpi(radar_name,{'kuband3','snow3'}))
      % counter_clk should be the EPRF (effective PRF after hardware presumming)
      % set anchor_idx to a record that you believe has the correct time
      counter_clk = 3906.250/8;
      counter = epri;
      keyboard
    elseif any(strcmpi(radar_name,{'snow8'}))
      % counter_clk should be the EPRF (effective PRF after hardware presumming)
      % set anchor_idx to a record that you believe has the correct time
      counter_clk = 4000/8;
      counter = epri;
      keyboard
    else
      keyboard
    end
    
    % Test example
    % counter_clk = 1;
    % counter_bad_threshold = 0.01;
    % counter_min_freq = 3;
    % counter_bin = 0.01;
    % anchor_idx = 5;
    % utc_time_sod = [4 5 6 10 11 12 15 16 17];
    % counter = [0 2:9];
    % anchor_idx = 5;
    
    % Find the differences of UTC and counter
    diff_utc = diff(double(utc_time_sod));
    diff_counter = diff(double(counter));
    
    % Find the valid counter differences based on frequency of occurance
    % Assumption is that invalid differences happen infrequently
    % (counter_min_freq)
     modes = unique(diff_counter);
     correctable_records = zeros(size(diff_counter));
     for idx=1:length(modes)
       if sum(diff_counter == modes(idx)) > counter_min_freq
         % Allow for a little slop in the differences (counter_bin)
         correctable_records = correctable_records ...
           | abs(diff_counter - modes(idx))/modes(idx) < counter_bin; 
       end
     end
     
     % Detect bad records when differences of counter and utc do not align
     % with the expected relationship (counter is running at counter_clk
     % frequency). Allow for some slop (counter_bad_threshold).
     bad_records = abs(diff_counter ./ diff_utc - counter_clk) > counter_clk*counter_bad_threshold;
     % Only correct records that had valid counter differences
     bad_records = bad_records & correctable_records;
     % Correct diff UTC using counter
     diff_utc(bad_records) = diff_counter(bad_records) / counter_clk;
     
     % Recreate UTC by integrating diff UTC around the anchor_idx
     % Integrate from the anchor_idx to 1
     first = cumsum(-diff_utc(anchor_idx-1:-1:1));
     first = first(end:-1:1);
     % Integrate from the anchor_idx to the end
     last = cumsum(diff_utc(anchor_idx:end));
    
     % Create the corrected utc_time_sod
     utc_time_sod_new = utc_time_sod(anchor_idx) + [first 0 last];
     
     figure(1); clf;
     plot(utc_time_sod);
     hold on;
     plot(utc_time_sod_new,'r');
     hold off;
     xlabel('Record');
     ylabel('UTC Time SOD (sec)');
     legend('Original','Corrected','location','best');
     figure(2); clf;
     plot(utc_time_sod - utc_time_sod_new);
     xlabel('Record');
     ylabel('Time correction (sec)');
     title('Ideally abs() is less than a few milliseconds');
     figure(3); clf;
     subplot(2,1,1);
     plot(diff(epri),'.');
     subplot(2,1,2);
     plot(diff(epri),'.');
     ylim([-3 5]);
     xlabel('Record');
     ylabel('Diff EPRI');
     title('Should be 1 except at segment boundaries');
     warning('Please check the corrected utc_time_sod (red) in figure 1 and the correction in figure 2. If correct, run "dbcont" to continue.');
     keyboard
     
     utc_time_sod = utc_time_sod_new;
  end

  % Check for day wraps in the UTC time seconds of day
  day_wrap_idxs = find(diff(utc_time_sod) < -50000);
  day_wrap_offset = zeros(size(utc_time_sod));
  for day_wrap_idx = day_wrap_idxs
    day_wrap_offset(day_wrap_idx+1:end) = day_wrap_offset(day_wrap_idx+1:end) + 86400;
  end
  utc_time_sod = utc_time_sod + day_wrap_offset;

  % Look for time gaps (this is used later for segmentation)
  time_gaps = find(abs(diff(utc_time_sod)) > MAX_TIME_GAP);

  % Look for EPRI gaps (this may be used later for segmentation)
  if union_time_epri_gaps
    MAX_EPRI_GAP = 2000;
    epri_gaps = find(abs(diff(medfilt1(double(epri),11))) > MAX_EPRI_GAP);
    time_gaps = sort(union(time_gaps,epri_gaps));
  end

  % Plot results
  figure(1); clf;
  plot(utc_time_sod);
  ylabel('UTC time seconds of day');
  xlabel('Record');
  grid on;
  hold on;
  plot(time_gaps, utc_time_sod(time_gaps),'ro');
  hold off;
elseif any(strcmpi(radar_name,{'acords'}))
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
  if param.file_version == 406
    change_fields = [2 5 6 7 8 9 10 11 12 17 18 19 20 21];
  elseif param.file_version == 405
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

%% Break into segments
if 1
  % Using time and optionally EPRI
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
    
    if segment_stop - segment_start + 1 >= MIN_SEG_SIZE
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
  if segment_stop - segment_start + 1 >= MIN_SEG_SIZE
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

elseif 0
  % Break into segments with EPRI
  % This code required for 2013 Antarctica P3 20131127 because time record is bad
  EPRI_JUMP_MIN = 0;
  EPRI_JUMP_MAX = 2e3;
  time_gaps = find(diff(epri) < EPRI_JUMP_MIN | diff(epri) > EPRI_JUMP_MAX);
  
  bad_mask = logical(zeros(size(fns_list{1})));
  segments = [];
  segment_start = file_idxs(1);
  finfo = fname_info_fmcw(fns_list{1});
  start_time = datenum_to_epoch(finfo.datenum);
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
    
    if segment_stop - segment_start + 1 >= MIN_SEG_SIZE
      seg_idx = seg_idx + 1;
      finfo = fname_info_fmcw(fns_list{1}{segment_start});
      segments(seg_idx).start_time = datenum_to_epoch(finfo.datenum);
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
    finfo = fname_info_fmcw(fns_list{1}{segment_start});
    start_time = datenum_to_epoch(finfo.datenum);
    start_day_wrap_offset = day_wrap_offset(time_gap+1);
  end
  segment_stop = length(fns_list{1});
  if segment_stop - segment_start + 1 >= MIN_SEG_SIZE
    seg_idx = seg_idx + 1;
    finfo = fname_info_fmcw(fns_list{1}{segment_start});
    segments(seg_idx).start_time = datenum_to_epoch(finfo.datenum);
    segments(seg_idx).start_idx = segment_start;
    segments(seg_idx).stop_idx = segment_stop;
    segments(seg_idx).day_wrap_offset = start_day_wrap_offset;
    [~,fn_start_name,fn_start_name_ext] = fileparts(fns_list{1}{segment_start});
    [~,fn_stop_name,fn_stop_name_ext] = fileparts(fns_list{1}{segment_stop});
    fprintf('%2d: %s %4d-%4d %s - %s\n', seg_idx, ...
      datestr(epoch_to_datenum(start_time)), segment_start, segment_stop,...
      [fn_start_name fn_start_name_ext], [fn_stop_name fn_stop_name_ext]);
  end

end

[~,sort_idxs] = sort(cell2mat({segments.start_time}));
segments = segments(sort_idxs);

fprintf('Done %s\n', datestr(now));

%% Vector worksheet of param spreadsheet print out
fprintf('Copy and paste the following into the vector worksheet of the parameter spreadsheet:\n');
fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'Date', ...
  'Segment', 'file.start_idx', 'file.stop_idx', 'file.basedir', 'file.adc_folder_name', 'file.prefix', 'file.midfix','file.regexp','gps.time_offset');
for seg_idx = 1:length(segments)
  fprintf('%s\t%02d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%g\n', day_string, ...
    seg_idx, segments(seg_idx).start_idx, segments(seg_idx).stop_idx, base_dir, param.adc_folder_name, file_prefix, file_midfix, file_regexp, param.gps_time_offset+segments(seg_idx).day_wrap_offset);
end

if any(strcmpi(radar_name,{'acords'}))
  % Print out some results that can be copied and pasted easily
  fprintf('\n')
  for seg_idx = 1:length(segments)
%     [hdr htime hoffset] = basic_load_acords(sprintf('%s/%s/%s.%d',base_dir,adc_folder_name,file_prefix_override,segments(seg_idx).start_idx-1),struct('datatype',0,'file_version',param.file_version,'verbose',0));
    [hdr htime hoffset] = basic_load_acords(fns_list{1}{segments(seg_idx).start_idx},struct('datatype',0,'file_version',param.file_version,'verbose',0));
    if param.file_version == 406
      if hdr(1).num_elem == 0
        fprintf('%s\t%02d\t%e\t%d\t12\t1\t2\t%4.4e\t\t\t%3.2e\t%3.2e\t\t0\t[1]\t[%d %d %d %d]\t10.^((44-%d*ones(1,1))/20)\t[0]\t[0]\t[0]/1e9\t%4.4e\t\t\t%3.2e\t%3.2e\t\t0\t[1]\t[%d %d %d %d]\t10.^((80-%d*ones(1,1))/20)\t[0]\t[0]\t[0]/1e9\n',...
        day_string,seg_idx,hdr(1).daq_clk,hdr(1).prf,hdr(1).tpd,hdr(1).f0,hdr(1).f1,hdr(1).elem_1+1,hdr(1).elem_2+1,hdr(1).elem_3+1,hdr(1).elem_4+1,hdr(1).low_gain_atten,...
        hdr(1).tpd,hdr(1).f0,hdr(1).f1,hdr(1).elem_1+1,hdr(1).elem_2+1,hdr(1).elem_3+1,hdr(1).elem_4+1,hdr(1).high_gain_atten);
      else
        fprintf('%s\t%02d\t%e\t%d\t12\t1\t2\t%4.4e\t\t\t%3.2e\t%3.2e\t\t0\t[1 1 1 1]\t[%d %d %d %d]\t10.^((44-%d*ones(1,4))/20)\t[0 0 0 0]\t[0 0 0 0]\t[0 0 0 0]/1e9\t%4.4e\t\t\t%3.2e\t%3.2e\t\t0\t[1 1 1 1]\t[%d %d %d %d]\t10.^((80-%d*ones(1,4))/20)\t[0 0 0 0]\t[0 0 0 0]\t[0 0 0 0]/1e9\n',...
        day_string,seg_idx,hdr(1).daq_clk,hdr(1).prf,hdr(1).tpd,hdr(1).f0,hdr(1).f1,hdr(1).elem_1+1,hdr(1).elem_2+1,hdr(1).elem_3+1,hdr(1).elem_4+1,hdr(1).low_gain_atten,...
        hdr(1).tpd,hdr(1).f0,hdr(1).f1,hdr(1).elem_1+1,hdr(1).elem_2+1,hdr(1).elem_3+1,hdr(1).elem_4+1,hdr(1).high_gain_atten);
      end
    elseif param.file_version == 405
      fprintf('%s\t%02d\t%e\t%d\t12\t1\t2\t%4.4e\t\t\t%3.2e\t%3.2e\t\t0\t[1]\t[1]\t10.^((44-%d*ones(1,1))/20)\t[0]\t[0]\t[0]/1e9\t%4.4e\t\t\t%3.2e\t%3.2e\t\t0\t[1]\t[1]\t10.^((80-%d*ones(1,1))/20)\t[0]\t[0]\t[0]/1e9\n',...
        day_string,seg_idx,hdr(1).daq_clk,hdr(1).prf,hdr(1).tpd,hdr(1).f0,hdr(1).f1,hdr(1).low_gain_atten,...
        hdr(1).tpd,hdr(1).f0,hdr(1).f1,hdr(1).high_gain_atten);
    end
  end
end

%% Save segment results
adc_idx = 1;
adc = adcs(adc_idx);
board = adc_to_board(param.radar_name,adc);
adc_folder_name = param.adc_folder_name;
adc_folder_name = regexprep(adc_folder_name,'%02d',sprintf('%02.0f',adc));
adc_folder_name = regexprep(adc_folder_name,'%d',sprintf('%.0f',adc));
adc_folder_name = regexprep(adc_folder_name,'%b',sprintf('%.0f',board));
tmp_hdr_fn = ct_filename_ct_tmp(param,'','headers', ...
  fullfile(adc_folder_name, 'create_segment_raw_file_list_v2.mat'));
fprintf('Saving %s\n', tmp_hdr_fn);
save(tmp_hdr_fn,'day_string','base_dir','param','file_prefix','file_midfix','file_regexp','segments');

return;


