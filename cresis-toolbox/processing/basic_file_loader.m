function [data,fn,settings,default,gps,hdr,pc_param] = basic_file_loader(param,defaults)
% [data,fn,settings,default,gps,hdr,pc_param] = basic_file_loader(param,defaults)
%
% Support script for loading raw data files for the in flight scripts. Used by
% basic_noise_analysis.m, basic_rx_equalization.m, and
% basic_radiometric_impulse_response.m
%
% Author: John Paden

global g_basic_noise_analysis_fn;
g_basic_noise_analysis_fn = char(g_basic_noise_analysis_fn); % Force type to be string

% Assume the first default parameters until we know which is the correct
default = defaults{1};

global g_file_search_mode;
if isempty(param.file_search_mode) ...
    || all(~strcmpi(param.file_search_mode,{'default','default+1','specific','last_file'}))
  fprintf('Select file search mode:\n');
  fprintf(' 1) default (loads whatever was last loaded)');
  if strcmpi(g_file_search_mode,'default')
    fprintf(' *');
  end
  fprintf('\n');
  fprintf(' 2) specific (lets you specify a specific file)');
  if strcmpi(g_file_search_mode,'specific')
    fprintf(' *');
  end
  fprintf('\n');
  fprintf(' 3) last_file (lists the last 10 recorded files to select from)');
  if strcmpi(g_file_search_mode,'last_file')
    fprintf(' *');
  end
  fprintf('\n');
  done = false;
  while ~done
    file_search_mode = input('Selection (1, 2, or 3): ');
    try
      if isempty(file_search_mode) && ~isempty(g_file_search_mode)
        param.file_search_mode = g_file_search_mode;
        done = true;
      elseif file_search_mode==1
        param.file_search_mode = 'default';
        done = true;
      elseif file_search_mode==2
        param.file_search_mode = 'specific';
        done = true;
      elseif file_search_mode==3
        param.file_search_mode = 'last_file';
        done = true;
      end
    end
  end
end
g_file_search_mode = param.file_search_mode;

%% Determine which file to load
if ~strncmpi(param.file_search_mode,'default',length('default'))
  good_mask = logical(zeros(size(param.base_dir_search)));
  for base_dir_idx = 1:length(param.base_dir_search)
    base_dir = param.base_dir_search{base_dir_idx};
    if exist(base_dir,'dir')
      good_mask(base_dir_idx) = true;
    end
  end
  param.base_dir_search = param.base_dir_search(good_mask);
  
  global g_basic_noise_analysis_base_dir;
  base_dir = [];
  if length(param.base_dir_search) >= 1
    default_base_dir_idx = [];
    if isempty(g_basic_noise_analysis_base_dir)
      g_basic_noise_analysis_base_dir = param.base_dir_search{1};
    end
    for base_dir_idx = 1:length(param.base_dir_search)
      fprintf('(%d): %s', base_dir_idx, param.base_dir_search{base_dir_idx});
      if strcmp(param.base_dir_search{base_dir_idx},g_basic_noise_analysis_base_dir)
        fprintf(' *');
        default_base_dir_idx = base_dir_idx;
      end
      fprintf('\n');
    end
    fprintf('(%d): Custom', base_dir_idx+1);
    if isempty(default_base_dir_idx)
      fprintf(' *');
    end
    fprintf('\n');
    base_dir_idx = input('More than one base directory exists, choose one: ');
    if isempty(base_dir_idx) && ~isempty(default_base_dir_idx)
      base_dir_idx = default_base_dir_idx;
    end
    if base_dir_idx <= length(param.base_dir_search)
      base_dir = param.base_dir_search{base_dir_idx};
    end
  end
  
  if isempty(base_dir)
    while ~exist(base_dir,'dir')
      base_dir = input('Enter custom directory path: ','s');
      if ~exist(base_dir,'dir')
        warning('Does not exist: %s', base_dir);
      end
    end
  end
  g_basic_noise_analysis_base_dir = base_dir;
  
  global g_file_select;
  if strcmpi(param.file_search_mode,'last_file')
    fns = get_filenames(base_dir,'','','.bin',struct('recursive',true));
    if isempty(fns)
      error('No data files: %s\n', base_dir);
    end
    fns_idxs = max(1,length(fns)-9) : length(fns);
    if isempty(g_file_select)
      g_file_select = max(1,length(fns_idxs)-1);
    end
    for fn_idx = 1:length(fns_idxs)
      fprintf('(%d): %s', fn_idx, fns{fns_idxs(fn_idx)});
      if g_file_select == fn_idx
        fprintf(' *');
      end
      fprintf('\n');
    end
    done = false;
    while ~done
      try
        user_fn_idx = input('Choose one: ');
        if isempty(user_fn_idx)
          user_fn_idx = g_file_select;
        end
        fn = fns{fns_idxs(user_fn_idx)};
        done = true;
      catch
        user_fn_idx = g_file_select;
        fn = fns{fns_idxs(user_fn_idx)};
        done = true;
      end
    end
    g_file_select = user_fn_idx;
  else
    fn = '';
    while ~exist(fn,'file')
      fn = input(sprintf('Filename [%s]: ', g_basic_noise_analysis_fn),'s');
      if isempty(fn)
        [~,fn] = fileparts(g_basic_noise_analysis_fn);
      end
      fn = get_filename(base_dir,'',fn,'',struct('recursive',true));
    end
  end
  g_basic_noise_analysis_fn = fn;
else
  if strcmp(param.file_search_mode,'default+1')
    % Get the next file index after the current file
    if strcmpi(param.radar_name,'mcords5')
      [fn_dir,fn_name,fn_ext] = fileparts(g_basic_noise_analysis_fn);
      cur_file_idx = str2double(fn_name(end-3:end));
      fn_name(end-3:end) = sprintf('%04d',cur_file_idx + 1);
      fn_name = [fn_name(1:end-14) '*' fn_name(end-7:end)]
      g_basic_noise_analysis_fn = get_filename(fn_dir,'',fn_name,fn_ext);
    end
  end
end

%% Load the chosen file(s)
tstart = tic;
fprintf('Loading data (%.1f sec)\n', toc(tstart));
% Load the data (disable if you have already loaded)
clear data;
clear num_rec;
if strcmpi(param.radar_name,'mcords')
  for adc_idx = 1:length(param.adcs)
    adc = param.adcs(adc_idx);

    % May need to adjust base_path for non-standard directory structures
    base_path = fullfile(param.base_path, sprintf('chan%d',adc), ...
      param.seg);
    file_midfix = sprintf('r%d-%d.',param.radar_num,adc);
    file_suffix = sprintf('.%04d.dat',param.data_file_num);
    fprintf('  Path: %s\n', base_path);
    fprintf('  Match: mcords*%s*%s\n', file_midfix, file_suffix);
    fn = get_filename(base_path,'mcords',file_midfix,file_suffix);
    if isempty(fn)
      fprintf('  Could not find any files which match\n');
      return;
    end
    fprintf('  Loading file %s\n', fn);
    [hdr,data_tmp] = basic_load_mcords(fn, struct('clk',1e9/9,'first_byte',2^26));
    data(:,:,adc_idx) = data_tmp{param.wf}(1:end-1,1:min(size(data_tmp{param.wf},2),param.rlines(2)));
  end
  data = data - median(data(:,1));
%   basic_remove_mcords_digital_errors;
elseif any(strcmpi(param.radar_name,{'mcords2','mcords3'}))
  % test1_1.dat0
  %   testA_N.datB
  %   A = acquisition number
  %   N = file number
  %   B = board number
  % Map ADCs to board numbers
  for board = 0:3
    if any(board == floor((param.img(:,2)-1)/4))
      get_adcs = board*4 + (1:4);
      file_prefix = sprintf('%s_%d_',param.radar_name,board);
      if isempty(acquisition_num)
        file_suffix = sprintf('%04d.bin',file_num);
      else
        file_suffix = sprintf('%s_%04d.bin',acquisition_num,file_num);
      end
      base_path = fullfile(param.base_path, sprintf('board%d',board), ...
        seg);
      fprintf('  Path: %s\n', base_path);
      fprintf('  Match: %s*%s\n', file_prefix, file_suffix);
      fn = get_filenames(base_path, file_prefix, '', file_suffix);
      if isempty(fn)
        fprintf('  Could not find any files which match\n');
        return;
      end
      fn = fn{end};
      fprintf('  Loading file %s\n', fn);
      % Fix get_filenames     'The filename, directory name, or volume label syntax is incorrect.'
      if strcmpi(param.radar_name,'mcords2')
        [hdr,data_tmp] = basic_load_mcords2(fn,struct('clk',fs));
      else
        [hdr,data_tmp] = basic_load_mcords3(fn,struct('clk',fs));
      end
      for get_adc_idx = 1:length(get_adcs)
        adc = get_adcs(get_adc_idx);
        for wf_adc_idx = find(param.img(:,2) == adc)
          wf = param.img(wf_adc_idx,1);
          fprintf('  Loading wf %d, adc %d\n', wf, adc);
          if ~exist('num_rec','var')
            % Since each file may have slightly different numbers of
            % records we do this
            num_rec = size(data_tmp{wf},2) - 1;
          end
          data(:,:,wf_adc_idx) = data_tmp{wf}(:,1:num_rec,get_adc_idx);
        end
      end
    end
  end
  
elseif any(strcmpi(param.radar_name,{'mcords4','mcords5'}))
  file_idx = 1;
  epri_intersect = [];
  % adcs: a list of the adcs that we are loading
  adcs = unique(param.img(:,2));
  
  for adc = reshape(adcs,[1 length(adcs)])
    
    [fn_dir,fn_name] = fileparts(g_basic_noise_analysis_fn);
    fn_dir = fileparts(fn_dir);
    fn_name(9:10) = sprintf('%02d',adc);
    fn = fullfile(fn_dir,sprintf('chan%d',adc),[fn_name,'.bin']);
    
    if ~exist(fn,'file')
      warning('File %s not found. Filename may have a slightly different time stamp, so trying a wild character search for the file.',fn);
      [fn_dir,fn_name] = fileparts(fn);
      fn = get_filename(fn_dir,fn_name(1:7),fn_name(end-6:end),'');
    end
    fprintf('  Loading file %s\n', fn);
    % Load the data file
    if strcmp(param.radar_name,'mcords4')
      [hdr,data_tmp] = basic_load_mcords4(fn,struct('clk',default.radar.fs/4,'recs',param.recs));
    else
      [hdr,data_tmp] = basic_load_mcords5(fn,struct('clk',default.radar.fs,'recs',param.recs,'presum_bug_fixed',1));
    end
    % Remove extra records to help reduce total memory usage
%     if isfield(param,'rlines') && ~isempty(param.rlines)
%       for wf = 1:length(data_tmp)
%         data_tmp{wf} = data_tmp{wf}(:,param.rlines,:);
%       end
%       hdr.utc_time_sod = hdr.utc_time_sod(param.rlines);
%       hdr.epri = hdr.epri(param.rlines);
%     end
    % Map each of the read waveforms needed into the correct output
    for wf_adc_idx = 1:size(param.img,1)
      % wf,adc: pair of values for this entry in param.img
      wf = param.img(wf_adc_idx,1);
      if wf > length(data_tmp)
        error('Requested waveform (%d) is larger than the number of waveforms in the file (%d in the file).', wf, length(data_tmp));
      end
      if adc == abs(param.img(wf_adc_idx,2));
        % This pair needs to be loaded, insert into output array... handle
        % mismatched EPRIs using intersect function. Throw away any records
        % that do not have matching EPRI in all channels.
        if isempty(epri_intersect)
          epri_intersect = hdr.epri;
          if imag(wf) == 0
%             data(:,1:size(data_tmp{wf},2),wf_adc_idx) = data_tmp{wf}(:,1:end);
            data(:,:,wf_adc_idx) = data_tmp{wf}(:,:);
          else
            data(:,:,wf_adc_idx) = data_tmp{abs(wf)}(:,:) ...
              + sign(wf) * data_tmp{abs(wf)+1}(:,:);
          end
          hdr_utc_time_sod = hdr.utc_time_sod;
        else
          [new_epri_intersect data_idx data_tmp_idx] = intersect(epri_intersect,hdr.epri);
          if isempty(new_epri_intersect)
            warning('No matching EPRI is this load. This script does not handle file recording offsets between channels very well. Just assuming all data is good.');
            new_epri_intersect = epri_intersect;
            data_idx = 1:length(epri_intersect);
            data_tmp_idx = 1:length(epri_intersect);
          end
          epri_intersect = new_epri_intersect;
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
  
  xml_version = 2.0;
  cresis_xml_mapping;
  
  %% Read XML files in this directory
  settings = read_ni_xml_directory(fn_dir,'',false);
  finfo = fname_info_mcords2(fn);
  
  settings_idx = find(cell2mat({settings.datenum}) < finfo.datenum,1,'last');
  if isempty(settings_idx)
    settings_idx = 1;
  end
  settings = settings(settings_idx);

  default = default_radar_params_settings_match(defaults,settings);
  
  %% Format settings into pc_param
  DDC_freq = double(settings.DDC_Ctrl.NCO_freq)*1e6;
  DDC_mode = double(settings.DDC_Ctrl.DDC_sel.Val);
  if DDC_mode == 0
    hdr.fs = default.radar.fs;
  else
    hdr.fs = default.radar.fs / 2^(1+DDC_mode);
  end
  f0 = settings.DDS_Setup.Waveforms(wf).Start_Freq(1);
  f1 = settings.DDS_Setup.Waveforms(wf).Stop_Freq(1);
  fc = (f0+f1)/2;
  Tpd = settings.DDS_Setup.Base_Len * double(settings.DDS_Setup.Waveforms(wf).Len_Mult);
  hdr.BW = f1-f0;
  if DDC_mode == 0
    hdr.BW_noise = 450e6;
  elseif DDC_mode == 1
    hdr.BW_noise = 350e6;
  elseif DDC_mode == 2
    hdr.BW_noise = 175e6;
  end
  atten = double(settings.DDS_Setup.Waveforms(wf).Attenuator_1(1) + settings.DDS_Setup.Waveforms(wf).Attenuator_2(1));
  hdr.rx_gain = default.radar.rx_gain .* 10.^(-atten/20);
  t0 = hdr.wfs(wf).t0 + default.radar.Tadc_adjust;
  tukey = settings.DDS_Setup.RAM_Taper;
  
  dt = 1/hdr.fs;
  Nt = size(data,1);
  clear pc_param;
  pc_param.DDC_mode = DDC_mode;
  pc_param.DDC_freq = DDC_freq;
  pc_param.f0 = f0;
  pc_param.f1 = f1;
  pc_param.Tpd = Tpd;
  pc_param.zero_pad = 1;
  pc_param.decimate = true;
  pc_param.window_func = @hanning;
  pc_param.time = t0 + (0:dt:(Nt-1)*dt).';
  pc_param.tukey = tukey;
  
  finfo = fname_info_mcords2(fn);
  [year,month,day] = datevec(finfo.datenum);
  hdr.radar_time = utc_to_gps(datenum_to_epoch(datenum(year,month,day,0,0,hdr.utc_time_sod))) + default.vectors.gps.time_offset;
  
  %% Read GPS files in this directory
  param.day_seg = sprintf('%04d%02d%02d_01',year,month,day);
  gps_fn = ct_filename_support(param,'','gps',1);
  if exist(gps_fn,'file')
    gps = load(gps_fn);
  else
    try
      gps_fns = get_filenames(fn_dir,'GPS_','','.txt');
      
      for gps_fn_idx = 1:length(gps_fns)
        gps_fn = gps_fns{gps_fn_idx};
        fprintf('  GPS file: %s\n', gps_fn);
        [~,gps_fn_name] = fileparts(gps_fn);
        gps_params.year = str2double(gps_fn_name(5:8));
        gps_params.month = str2double(gps_fn_name(9:10));
        gps_params.day = str2double(gps_fn_name(11:12));
        gps_params.time_reference = 'utc';
        if gps_fn_idx == 1
          gps = read_gps_nmea(gps_fns{gps_fn_idx},gps_params);
        else
          gps_tmp = read_gps_nmea(gps_fns{gps_fn_idx},gps_params);
          gps.gps_time = [gps.gps_time, gps_tmp.gps_time];
          gps.lat = [gps.lat, gps_tmp.lat];
          gps.lon = [gps.lon, gps_tmp.lon];
          gps.elev = [gps.elev, gps_tmp.elev];
          gps.roll = [gps.roll, gps_tmp.roll];
          gps.pitch = [gps.pitch, gps_tmp.pitch];
          gps.heading = [gps.heading, gps_tmp.heading];
        end
      end
    end
  end
  hdr.gps_time = hdr.radar_time;
  try
    hdr.lat = interp1(gps.gps_time,gps.lat,hdr.radar_time);
    hdr.lon = interp1(gps.gps_time,gps.lon,hdr.radar_time);
    hdr.elev = interp1(gps.gps_time,gps.elev,hdr.radar_time);
    hdr.roll = interp1(gps.gps_time,gps.roll,hdr.radar_time);
    hdr.pitch = interp1(gps.gps_time,gps.pitch,hdr.radar_time);
    hdr.heading = interp1(gps.gps_time,gps.heading,hdr.radar_time);
    hdr.gps_source = gps.gps_source;
  catch ME
    ME.getReport
    hdr.lat = zeros(size(hdr.radar_time));
    hdr.lon = zeros(size(hdr.radar_time));
    hdr.elev = zeros(size(hdr.radar_time));
    hdr.roll = zeros(size(hdr.radar_time));
    hdr.pitch = zeros(size(hdr.radar_time));
    hdr.heading = zeros(size(hdr.radar_time));
    hdr.gps_source = '';
  end
  
end
