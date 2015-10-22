% script basic_rx_chan_equalization
%
% RUN THIS FUNCTION FROM "run_basic_rx_chan_equalization"
%
% This script is for helping with setting the receiver coefficients
% from raw data. It requires loading one waveform and N receive channels
% and then analyzing these data.
%
% 1. Collect data with so that receiver gain settings are the same as
%    the ordinary data collection, but the air/ice surface response
%    is unsaturated (ideally from high altitude over the ocean)
%    If the TR switch is slow, you will need to collect the data from
%    a high enough altitude where the fast-time gain changes of the TR
%    switch are stable (this is 10 us after the switch change for the
%    original MCoRDS system).
% 2. Pulse duration should not matter, but probably affects results
%    slightly
% 3. Time gate should be large enough to include noise-only data which
%    will be used to determine the SNR of the surface return (SNR
%    threshold is used to exclude low SNR data points from the measurement)
%
% Author: John Paden, Logan Smith
%
% See Also: run_basic_rx_chan_equalization

physical_constants;
basic_rx_chan_equalization_tstart = tic;

clear td_out phase_out amp_out;

for file_idx = 1:length(file_nums)
  file_num = file_nums(file_idx);
  
  % =======================================================================
  % Load data
  % =======================================================================
  fprintf('========================================================\n');
  fprintf('Loading data %d (%.1f sec)\n', file_num, toc(basic_rx_chan_equalization_tstart));
  % Load the data (disable if you have already loaded)
  clear data;
  clear num_rec;
  epri_intersect = [];
  if strcmpi(radar_name,'mcords')
    adcs = unique(param.img(:,2));
    for adc_idx = 1:length(adcs)
      adc1 = adcs(adc_idx);
      
      % May need to adjust base_path for non-standard directory structures
      fn_dir = fullfile(base_path, sprintf(chan,adc1), seg);
      fn = get_filename(fn_dir,'mcords','','*.dat');
      fname = fname_info_mcords(fn);
      if fname.rev == 2
        file_suffix = sprintf('.%04d.r%d-%d.*.dat',file_num,radar_id,adc1);
        fprintf('  Path: %s\n', fn_dir);
        fprintf('  Match: mcords*%s*%s\n', '', file_suffix);
        fn = get_filename(fn_dir,fn_start,'',file_suffix);
      elseif fname.rev == 3
        file_midfix = sprintf('r%d-%d.',radar_id,adc1);
        file_suffix = sprintf('.%04d.dat',file_num);
        fprintf('  Path: %s\n', fn_dir);
        fprintf('  Match: mcords*%s*%s\n', file_midfix, file_suffix);
        fn = get_filename(fn_dir,fn_start,file_midfix,file_suffix);
      else
        error('Unsupported file version');
      end
      file_name_list{file_idx} = fn;
      if isempty(fn)
        fprintf('  Could not find any files which match\n');
        return;
      end
      [hdr,data_tmp] = basic_load_mcords(fn, struct('clk',1e9/9,'first_byte',2^26));
      
      for wf_adc_idx = 1:size(param.img,1)
        % wf,adc: pair of values for this entry in param.img
        wf = param.img(wf_adc_idx,1);
        adc = param.img(wf_adc_idx,2);
        % Check to see if this pair was loaded
        if adc == adc1
          if isempty(epri_intersect)
            epri_intersect = hdr.epri;
            data(:,:,wf_adc_idx) = data_tmp{wf};
            hdr_utc_time_sod = hdr.utc_time_sod;
          else
            [epri_intersect data_idx data_tmp_idx] = intersect(epri_intersect,hdr.epri);
            data = data(:,data_idx,:);
            hdr_utc_time_sod = hdr.utc_time_sod(data_tmp_idx);
            data(:,:,wf_adc_idx) = data_tmp{wf}(:,data_tmp_idx);
          end
        end
      end
    end
    hdr.utc_time_sod = hdr_utc_time_sod;
    % Trim away first and last samples (e.g. last sample is known to always be bad)
    data = data(trim(1):end-trim(2),:,:);
    % Remove coherent noise on center elements
    %     for wf_adc_idx = 1:size(data,3)
    %       adc = param.img(wf_adc_idx,2);
    %       if adc >= 2 && adc <= 4
    %         data(:,:,wf_adc_idx) = data(:,:,wf_adc_idx) - repmat(mean(data(:,:,wf_adc_idx),2),[1 size(data,2)]);
    %       end
    %     end
    
  elseif strcmpi(radar_name,'mcords2')
    % adcs: a list of the adcs that we are loading
    adcs = unique(param.img(:,2));
    
    % Consider each ADC board
    for board = 0:3
      % Check to see if any adcs are needed from this board
      if any(board == floor((adcs-1)/4))
        % get_adcs: The list of adcs associated with this board
        get_adcs = board*4 + (1:4);
        % Build file path and load
        file_prefix = sprintf('mcords2_%d_',board);
        if isempty(acquisition_num)
          file_suffix = sprintf('%04d.bin',file_num);
        else
          file_suffix = sprintf('%02d_%04d.bin',acquisition_num,file_num);
        end
        fn_dir = fullfile(base_path, sprintf('board%d',board), seg);
        fprintf('  Path: %s\n', fn_dir);
        fprintf('  Match: %s*%s\n', file_prefix, file_suffix);
        fn = get_filenames(fn_dir, file_prefix, '', file_suffix);
        if isempty(fn)
          error('  Could not find any files which match %s/%s*%s\n', fn_dir, file_prefix, file_suffix);
        elseif length(fn) == 1
          fn = fn{1};
        else
          fprintf('Select a specific file (fn = fn{1}), then dbcont\n');
          fn
          %keyboard
          fn = fn{3}
        end
        file_name_list{file_idx} = fn;
        fprintf('  Loading file %s\n', fn);
        % Load the data file
        [hdr,data_tmp] = basic_load_mcords2(fn,struct('clk',fs));
        % Go through each adc that was loaded and store it in the output
        % matrix as required by param.adcs
        for get_adc_idx = 1:length(get_adcs)
          for wf_adc_idx = 1:size(param.img,1)
            % wf,adc: pair of values for this entry in param.img
            wf = param.img(wf_adc_idx,1);
            adc = param.img(wf_adc_idx,2);
            % Check to see if this pair was loaded
            if adc == get_adcs(get_adc_idx)
              % This pair was loaded, insert into output array... handle
              % mismatched EPRIs using intersect function.
              if isempty(epri_intersect)
                epri_intersect = hdr.epri;
                data(:,:,wf_adc_idx) = data_tmp{wf}(:,:,get_adc_idx);
                hdr_utc_time_sod = hdr.utc_time_sod;
              else
                [epri_intersect data_idx data_tmp_idx] = intersect(epri_intersect,hdr.epri);
                data = data(:,data_idx,:);
                hdr_utc_time_sod = hdr.utc_time_sod(data_tmp_idx);
                data(:,:,wf_adc_idx) = data_tmp{wf}(:,data_tmp_idx,get_adc_idx);
              end
            end
          end
        end
      end
    end
    hdr.utc_time_sod = hdr_utc_time_sod;
    
  elseif strcmpi(radar_name,'mcords3')
    % adcs: a list of the adcs that we are loading
    adcs = unique(param.img(:,2));
    % Consider each ADC board
    for board = 0:3
      % Check to see if any adcs are needed from this board
      if any(board == floor((adcs-1)/4))
        % get_adcs: The list of adcs associated with this board
        get_adcs = board*4 + (1:4);
        % Build file path and load
        file_prefix = sprintf('mcords3_%d_',board);
        if isempty(acquisition_num)
          file_suffix = sprintf('%04d.bin',file_num);
        else
          file_suffix = sprintf('%02d_%04d.bin',acquisition_num,file_num);
        end
        fn_dir = fullfile(base_path, sprintf('board%d',board), seg);
        fprintf('  Path: %s\n', fn_dir);
        fprintf('  Match: %s*%s\n', file_prefix, file_suffix);
        fn = get_filenames(fn_dir, file_prefix, '', file_suffix);
        if isempty(fn)
          error('  Could not find any files which match\n');
        elseif length(fn) == 1
          fn = fn{1};
        else
          fprintf('Select a specific file (fn = fn{1}), then dbcont\n');
          fn
          keyboard
        end
        file_name_list{file_idx} = fn;
        fprintf('  Loading file %s\n', fn);
        % Load the data file
        [hdr,data_tmp] = basic_load_mcords3(fn,struct('clk',fs));
        % Go through each adc that was loaded and store it in the output
        % matrix as required by param.adcs
        for get_adc_idx = 1:length(get_adcs)
          for wf_adc_idx = 1:size(param.img,1)
            % wf,adc: pair of values for this entry in param.img
            wf = param.img(wf_adc_idx,1);
            adc = param.img(wf_adc_idx,2);
            % Check to see if this pair was loaded
            if adc == get_adcs(get_adc_idx)
              % This pair was loaded, insert into output array... handle
              % mismatched EPRIs using intersect function.
              if isempty(epri_intersect)
                epri_intersect = hdr.epri;
                data(:,:,wf_adc_idx) = data_tmp{wf}(:,:,get_adc_idx);
                hdr_utc_time_sod = hdr.utc_time_sod;
              else
                [epri_intersect data_idx data_tmp_idx] = intersect(epri_intersect,hdr.epri);
                data = data(:,data_idx,:);
                hdr_utc_time_sod = hdr.utc_time_sod(data_tmp_idx);
                data(:,:,wf_adc_idx) = data_tmp{wf}(:,data_tmp_idx,get_adc_idx);
              end
            end
          end
        end
      end
    end
    hdr.utc_time_sod = hdr_utc_time_sod;
    
  elseif strcmpi(radar_name,'mcords4')
    % adcs: a list of the adcs that we are loading
    adcs = unique(param.img(:,2));
    
    for adc = reshape(adcs,[1 length(adcs)])
      file_prefix = sprintf('mcords4_%02d_',adc);
      if isempty(acquisition_num)
        file_suffix = sprintf('%04d.bin',file_num);
      else
        file_suffix = sprintf('%s_%04d.bin',acquisition_num,file_num);
      end
      fn_dir = fullfile(base_path, sprintf('chan%d',adc), seg);
      fprintf('  Path: %s\n', fn_dir);
      fprintf('  Match: %s*%s\n', file_prefix, file_suffix);
      fn = get_filename(fn_dir, file_prefix, '', file_suffix);
      file_name_list{file_idx} = fn;
      fprintf('  Loading file %s\n', fn);
      % Load the data file
      [hdr,data_tmp] = basic_load_mcords4(fn,struct('clk',fs/4));
      % Map each of the read waveforms needed into the correct output
      for wf_adc_idx = 1:size(param.img,1)
        % wf,adc: pair of values for this entry in param.img
        wf = param.img(wf_adc_idx,1);
        if adc == abs(param.img(wf_adc_idx,2));
          % This pair needs to be loaded, insert into output array... handle
          % mismatched EPRIs using intersect function.
          if isempty(epri_intersect)
            epri_intersect = hdr.epri;
            if imag(wf) == 0
              data(:,:,wf_adc_idx) = data_tmp{wf}(:,:);
            else
              data(:,:,wf_adc_idx) = data_tmp{abs(wf)}(:,:) ...
                + sign(wf) * data_tmp{abs(wf)+1}(:,:);
            end
            hdr_utc_time_sod = hdr.utc_time_sod;
          else
            [epri_intersect data_idx data_tmp_idx] = intersect(epri_intersect,hdr.epri);
            data = data(:,data_idx,:);
            hdr_utc_time_sod = hdr.utc_time_sod(data_tmp_idx);
            if imag(wf) == 0
              data(:,:,wf_adc_idx) = data_tmp{wf}(:,data_tmp_idx);
            else
              data(:,:,wf_adc_idx) = data_tmp{abs(wf)}(:,data_tmp_idx) ...
                + sign(wf) * data_tmp{abs(wf)+1}(:,data_tmp_idx);
            end
          end
        end
      end
    end
    hdr.utc_time_sod = hdr_utc_time_sod;
    
  elseif strcmpi(radar_name,'mcrds')
    % May need to adjust base_path for non-standard directory structures
    fn_dir = fullfile(base_path,file_adc_folder_name);
    file_suffix = sprintf('%04d.raw',file_num);
    
        fn = get_filenames(fn_dir, fn_prefix, '', file_suffix);
        if isempty(fn)
          error('  Could not find any files which match %s/%s*%s\n', fn_dir, fn_prefix, file_suffix);
        elseif length(fn) == 1
          fn = fn{1};
        else
          fprintf('Select a specific file (fn = fn{1}), then dbcont\n');
          fn
          %keyboard
          fn = fn{3}
        end
        file_name_list{file_idx} = fn;    
    
    [hdr,data_tmp] = basic_load_mcrds(fn);
    for adc_idx = 1:length(unique(param.img(:,2)))
      for wf_adc_idx = 1:size(param.img,1)
        % wf,adc: pair of values for this entry in param.img
        wf = param.img(wf_adc_idx,1);
        adc = param.img(wf_adc_idx,2);
        data(:,:,wf_adc_idx) = data_tmp.wf{wf}(:,:,adc);
      end
    end
    
    % Subtract average from data (use median as a robust way to find the
    % mean)
    data = data - median(data(:,1,1));
  end
  clear data_tmp;
  fprintf('  %s\n', fn);
  
  if strcmpi(radar_name,'mcords')

    if 0
      % =======================================================================
      % Lots of bad records (typical of transmit calibration data)
      % =======================================================================
      
      % Subtract average from data (use median as a robust way to find the
      % mean)
      data = data - median(data(:,1,1));
      
      end_max_threshold = 500;
      peak_max_threshold = 3e4;

      end_max_collect = [];
      for wf_adc_idx = 1:size(data,3)
        end_max_good_idxs ...
          = intersect(intersect(find(max(abs(data(end-100:end,:,wf_adc_idx))) < end_max_threshold).', ...
          find(max(abs(data(1:100,:,wf_adc_idx))) < end_max_threshold).'), ...
          find(max(abs(data(:,:,wf_adc_idx))) < peak_max_threshold).');
        if isempty(end_max_collect)
          end_max_collect = end_max_good_idxs;
        else
          end_max_collect = intersect(end_max_good_idxs,end_max_collect);
        end
      end
      data2(:,:,:) = data(:,end_max_collect,:);
      hdr_utc_time_sod = hdr.utc_time_sod(end_max_collect);
      hdr.utc_time_sod = hdr_utc_time_sod;
      data = data2;
      clear data2 hdr_utc_time_sod
    else
      % =======================================================================
      % Digital error removal only
      % =======================================================================
      digital_noise = [];
      for wf_adc_idx = 1:size(data,3)
        for rline = 1:size(data,2)
          noise_idx = strfind(data(:,rline,wf_adc_idx).',[0 13415 0]);
          if ~isempty(noise_idx)
            digital_noise = cat(1,digital_noise,[noise_idx rline, wf_adc_idx]);
            %fprintf('%d %d %d %d\n', adc, rline, noise_idx, mod(data(noise_idx-1,rline,wf_adc_idx),256));
            %figure;
            %plot(data(noise_idx-1:noise_idx+2,rline,adc));
          end
        end
      end
      
      data = data - median(data(:,1,1));
      
      for row = 1:size(digital_noise)
        noise_idx = digital_noise(row,1);
        rline = digital_noise(row,2);
        wf_adc_idx = digital_noise(row,3);
        data(noise_idx-1:noise_idx+2,rline,wf_adc_idx) = 0;
      end
    end
    
    % =======================================================================
    % Presumming/coherent averaging
    % =======================================================================
    if presums > 1
      fprintf('Coherent averaging (%.1f sec)\n', toc(basic_rx_chan_equalization_tstart));
      clear data_out;
      for adc_idx = 1:size(data,3)
        data_out(:,:,adc_idx) = fir_dec(data(:,:,adc_idx),presums);
      end
      data = data_out;
      clear data_out;
    end
    
    % =========================================================================
    % Pulse Compression
    % =========================================================================
    fprintf('Pulse compression (%.1f sec)\n', toc(basic_rx_chan_equalization_tstart));
    pc_param.decimate = true;
    pc_param.zero_pad = true;
    orig_data_size = size(data);
    clear data_pc;
    for wf_adc_idx = 1:size(data,3)
      % wf,adc: The wf/adc pair of the current wf_adc_idx
      wf = param.img(wf_adc_idx,1);
      adc = param.img(wf_adc_idx,2);
      % pc_param.time: fast time axis for this set of raw data
      pc_param.time = hdr.wfs(wf).t0 + (0:size(data,1)-1)/fs;
      % Pulse compress
      [data_pc(:,:,wf_adc_idx),hdr.wfs(wf).time,hdr.wfs(wf).freq] = pulse_compress(data(:,:,wf_adc_idx),pc_param);
    end
    clear data;
    
    % =========================================================================
    % Select ROI using ginput
    % =========================================================================
    figure(101); clf;
    imagesc(lp(data_pc(:,:,param.ref_wf_adc_idx)));
    if use_ginput
      fprintf('\nSelect Signal Region (2 clicks)\n');
      title('Select Signal Region (2 clicks)');
      [param.rlines,param.rbins] = ginput(2);
      fprintf('\nSelect Noise Region (2 clicks)\n');
      title('Select Noise Region (2 clicks)');
      [param.noise_rlines,param.noise_rbins] = ginput(2);
    end
    
  elseif any(strcmpi(radar_name,{'mcords2','mcords3','mcords4'}))
    % =======================================================================
    % Presumming/coherent averaging
    % =======================================================================
    if presums > 1
      fprintf('Coherent averaging (%.1f sec)\n', toc(basic_rx_chan_equalization_tstart));
      clear data_out;
      for adc_idx = 1:size(data,3)
        data_out(:,:,adc_idx) = fir_dec(data(:,:,adc_idx),presums);
      end
      data = data_out;
      clear data_out;
    end
    
    % =========================================================================
    % Pulse Compression
    % =========================================================================
    fprintf('Pulse compression (%.1f sec)\n', toc(basic_rx_chan_equalization_tstart));
    pc_param.decimate = true;
    pc_param.zero_pad = true;
    orig_data_size = size(data);
    clear data_pc;
    for wf_adc_idx = 1:size(data,3)
      % wf,adc: The wf/adc pair of the current wf_adc_idx
      wf = abs(param.img(wf_adc_idx,1));
      adc = param.img(wf_adc_idx,2);
      % pc_param.time: fast time axis for this set of raw data
      pc_param.time = hdr.wfs(wf).t0 + (0:size(data,1)-1)/fs;
      % Pulse compress
      [data_pc(:,:,wf_adc_idx),hdr.wfs(wf).time,hdr.wfs(wf).freq] = pulse_compress(data(:,:,wf_adc_idx),pc_param);
    end
    clear data;
    
    % =========================================================================
    % Select ROI using ginput
    % =========================================================================
    figure(101); clf;
    imagesc(lp(data_pc(:,:,param.ref_wf_adc_idx)));
    if use_ginput
      fprintf('\nSelect Signal Region (2 clicks)\n');
      title('Select Signal Region (2 clicks)');
      [param.rlines,param.rbins] = ginput(2);
      fprintf('\nSelect Noise Region (2 clicks)\n');
      title('Select Noise Region (2 clicks)');
      [param.noise_rlines,param.noise_rbins] = ginput(2);
    end
    
  elseif strcmpi(radar_name,'mcrds')
    % =======================================================================
    % Presumming/coherent averaging
    % =======================================================================
    if presums > 1
      fprintf('Coherent averaging (%.1f sec)\n', toc(basic_rx_chan_equalization_tstart));
      clear data_out;
      for adc_idx = 1:size(data,3)
        data_out(:,:,adc_idx) = fir_dec(data(:,:,adc_idx),presums);
      end
      data = data_out;
      clear data_out;
    end
    
    % =========================================================================
    % Pulse Compression
    % =========================================================================
    fprintf('Pulse compression (%.1f sec)\n', toc(basic_rx_chan_equalization_tstart));
    pc_param.decimate = true;
    pc_param.zero_pad = true;
    orig_data_size = size(data);
    clear data_pc;
    for wf_adc_idx = 1:size(data,3)
      % wf,adc: The wf/adc pair of the current wf_adc_idx
      wf = param.img(wf_adc_idx,1);
      adc = param.img(wf_adc_idx,2);
      % pc_param.time: fast time axis for this set of raw data
      pc_param.time = (hdr.wf(wf).sampleDelayCount(adc)-10)/hdr.sampFreq + (0:size(data,1)-1)/fs;
      % Pulse compress
      [data_pc(:,:,wf_adc_idx),hdr.wfs(wf).time,hdr.wfs(wf).freq] = pulse_compress(data(:,:,wf_adc_idx),pc_param);
    end
    clear data;
    
    % =========================================================================
    % Select ROI using ginput
    % =========================================================================
    figure(101); clf;
    imagesc(lp(data_pc(:,:,param.ref_wf_adc_idx)));
    if use_ginput
      fprintf('\nSelect Signal Region (2 clicks)\n');
      title('Select Signal Region (2 clicks)');
      [param.rlines,param.rbins] = ginput(2);
      fprintf('\nSelect Noise Region (2 clicks)\n');
      title('Select Noise Region (2 clicks)');
      [param.noise_rlines,param.noise_rbins] = ginput(2);
    end
    
  end
  
  % =======================================================================
  % Load GPS
  % =======================================================================
  if ~isempty(param.gps_fn)
    fprintf('Loading GPS (%.1f sec)\n', toc(basic_rx_chan_equalization_tstart));
    gps = load(param.gps_fn);
    
    if strcmpi(radar_name,'mcrds')
      % Isolate the section of radar time from gps.radar_time that will be
      % used to interpolate with.
      guard_time = 5;
      good_idxs = find(gps.comp_time >= hdr.comp_time(1)-guard_time ...
        & gps.comp_time <= hdr.comp_time(end)+guard_time);
      good_radar_time = gps.radar_time(good_idxs);
      good_sync_gps_time = gps.sync_gps_time(good_idxs);
      % From these good indexes, remove any repeat radar times (usually caused
      % by there being more than one NMEA string every 1 PPS
      good_idxs = 1+find(diff(good_radar_time) ~= 0);
      good_radar_time = good_radar_time(good_idxs);
      good_sync_gps_time = good_sync_gps_time(good_idxs);
      
      
      if exist('utc_time_correction','var')
        if ischar(utc_time_correction)
          % records file
          records = load(utc_time_correction);
          % Apply correction
          hdr.radar_time_corr = hdr.radar_time ...
            + interp1(records.raw.radar_time-records.raw.radar_time_correction, ...
            records.raw.radar_time_correction, hdr.radar_time);
        else
          % double scalar offset
          hdr.radar_time_corr = hdr.radar_time + utc_time_correction;
        end
      else
        hdr.radar_time_corr = hdr.radar_time;
      end
      
      % Interpolate gps.sync_gps_time to records.gps_time using gps.radar_time
      % and records.radar_time
      hdr.gps_time = interp1(good_radar_time, good_sync_gps_time, ...
        hdr.radar_time_corr,'linear','extrap');
      
      hdr.lat = double(interp1(gps.gps_time,gps.lat,hdr.gps_time));
      hdr.lon = double(mod(interp1(gps.gps_time,unwrap(gps.lon/180*pi),hdr.gps_time)*180/pi+180, 360)-180);
      hdr.elev = double(interp1(gps.gps_time,gps.elev,hdr.gps_time));
      hdr.roll = double(interp1(gps.gps_time,gps.roll,hdr.gps_time));
      hdr.pitch = double(interp1(gps.gps_time,gps.pitch,hdr.gps_time));
      hdr.heading = double(mod(interp1(gps.gps_time,unwrap(gps.heading),hdr.gps_time)+pi,2*pi)-pi);
      
    else
      finfo = fname_info_mcords2(fn);
      [year,month,day] = datevec(finfo.datenum);
      hdr.gps_time = datenum_to_epoch(datenum(year,month,day,0,0,hdr.utc_time_sod)); % Still UTC time
      hdr.gps_time = hdr.gps_time + utc_leap_seconds(hdr.gps_time(1)); % Convert from UTC to GPS
      
      if exist('utc_time_correction','var')
        hdr.gps_time_corr = hdr.gps_time + utc_time_correction;
      else
        hdr.gps_time_corr = hdr.gps_time;
      end

      hdr.lat = interp1(gps.gps_time, gps.lat, hdr.gps_time_corr);
      hdr.lon = interp1(gps.gps_time, gps.lon, hdr.gps_time_corr);
      hdr.elev = interp1(gps.gps_time, gps.elev, hdr.gps_time_corr);
      hdr.roll = interp1(gps.gps_time, gps.roll, hdr.gps_time_corr);
      hdr.pitch = interp1(gps.gps_time, gps.pitch, hdr.gps_time_corr);
      hdr.heading = interp1(gps.gps_time, gps.heading, hdr.gps_time_corr);
    end
    hdr.gps_time = fir_dec(hdr.gps_time,presums);
    hdr.lat = fir_dec(hdr.lat,presums);
    hdr.lon = fir_dec(hdr.lon,presums);
    hdr.elev = fir_dec(hdr.elev,presums);
    hdr.roll = fir_dec(hdr.roll,presums);
    hdr.pitch = fir_dec(hdr.pitch,presums);
    hdr.heading = fir_dec(hdr.heading,presums);
    hdr.gps_source = gps.gps_source;

    figure(100); clf;
    plot(gps.gps_time,gps.roll*180/pi,'r');
    hold on;
    plot(hdr.gps_time,hdr.roll*180/pi);
    hold off;
    xlim([hdr.gps_time(1)-50 hdr.gps_time(end)+50])
    grid on;
    ylabel('Roll (deg)');
  else
    hdr.gps_time = zeros(1,size(data_pc,2));
    hdr.lat = zeros(1,size(data_pc,2));
    hdr.lon = zeros(1,size(data_pc,2));
    hdr.elev = zeros(1,size(data_pc,2));
    hdr.roll = zeros(1,size(data_pc,2));
    hdr.pitch = zeros(1,size(data_pc,2));
    hdr.heading = zeros(1,size(data_pc,2));
    hdr.gps_source = '';
  end
  
  % Perform receiver channel equalization
  [td_out(:,file_idx),amp_out(:,file_idx),phase_out(:,file_idx), full_out] = rx_chan_equal(data_pc,param,hdr);
  
  if file_idx == 1
    peak_ref = full_out.peak_ref;
    roll = hdr.roll;
    gps_time = hdr.gps_time;
    peak_offset = full_out.peak_offset;
    roll_est_theta = full_out.roll_est_theta;
    roll_est_val = full_out.roll_est_val;
    roll_est_gps_time = full_out.gps_time;
  else
    peak_ref = cat(2,peak_ref,full_out.peak_ref);
    roll = cat(2,roll,hdr.roll);
    gps_time = cat(2,gps_time,hdr.gps_time);
    peak_offset = cat(2,peak_offset,full_out.peak_offset);
    roll_est_theta = cat(2,roll_est_theta,full_out.roll_est_theta);
    roll_est_val = cat(2,roll_est_val,full_out.roll_est_val);
    roll_est_gps_time = cat(2,roll_est_gps_time,full_out.gps_time);
  end
  
end

fprintf('========================================================\n');
fprintf('Recommended equalization coefficients (averaged results)\n');

fprintf('  Date of processing: %s, mocomp %d, wf/adc %d/%d bins %d-%d\n', ...
  datestr(now), param.mocomp_type, param.img(param.ref_wf_adc_idx,1), ...
  param.img(param.ref_wf_adc_idx,2), param.rbins(1), param.rbins(end));
fprintf('td settings\n');
for file_idx = 1:length(file_nums)
  [~,fn] = fileparts(file_name_list{file_idx});
  fprintf('%s', fn);
  fprintf('\t%.2f', td_out(:,file_idx)*1e9);
  fprintf('\n');
end
fprintf('amp settings\n');
for file_idx = 1:length(file_nums)
  [~,fn] = fileparts(file_name_list{file_idx});
  fprintf('%s', fn);
  fprintf('\t%.1f', amp_out(:,file_idx));
  fprintf('\n');
end
fprintf('phase settings\n');
for file_idx = 1:length(file_nums)
  [~,fn] = fileparts(file_name_list{file_idx});
  fprintf('%s', fn);
  fprintf('\t%.1f', phase_out(:,file_idx));
  fprintf('\n');
end

td_ave = mean(td_out,2);
amp_ave = mean(amp_out,2);
phase_ave = angle(mean(exp(j*phase_out/180*pi),2))*180/pi;

fprintf('Rx Path\n');
for wf_adc_idx = 1:size(param.img,1)
  wf = abs(param.img(wf_adc_idx,1));
  adc = param.img(wf_adc_idx,2);
  if wf_adc_idx < size(param.img,1)
    fprintf('%d\t', param.rx_paths{wf}(adc));
  else
    fprintf('%d', param.rx_paths{wf}(adc));
  end
end
fprintf('\n');

fprintf('Original/Recommended/Difference td settings (ns):\n');
fprintf('%.2f\t', param.td(1:end-1)*1e9);
fprintf('%.2f', param.td(end)*1e9);
fprintf('\n');
fprintf('%.2f\t', td_ave(1:end-1)*1e9);
fprintf('%.2f', td_ave(end)*1e9);
fprintf('\n');
fprintf('%.2f\t', (td_ave(1:end-1) - param.td(1:end-1))*1e9);
fprintf('%.2f', (td_ave(end) - param.td(end))*1e9);
fprintf('\n');

fprintf('Original/Recommended/Difference amp settings (dB):\n');
fprintf('%.1f\t', param.amp(1:end-1));
fprintf('%.1f', param.amp(end));
fprintf('\n');
fprintf('%.1f\t', amp_ave(1:end-1));
fprintf('%.1f', amp_ave(end));
fprintf('\n');
fprintf('%.1f\t', amp_ave(1:end-1) - param.amp(1:end-1));
fprintf('%.1f', amp_ave(end) - param.amp(end));
fprintf('\n');

% Rewrap each phase so that the output does not print +355 deg and -5 deg
fprintf('Original/Recommended/Difference phase settings (deg):\n');
phase_rewrapped = angle(exp(j*param.phase/180*pi)) * 180/pi;
fprintf('%.1f\t', phase_rewrapped(1:end-1));
fprintf('%.1f', phase_rewrapped(end));
fprintf('\n');
phase_rewrapped = angle(exp(j*phase_ave/180*pi)) * 180/pi;
fprintf('%.1f\t', phase_rewrapped(1:end-1));
fprintf('%.1f', phase_rewrapped(end));
fprintf('\n');
phase_rewrapped = angle(exp(j*(phase_ave - param.phase)/180*pi)) * 180/pi;
fprintf('%.1f\t', phase_rewrapped(1:end-1));
fprintf('%.1f', phase_rewrapped(end));
fprintf('\n');

return;



