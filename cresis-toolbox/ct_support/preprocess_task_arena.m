function success = preprocess_task_arena(param)
% success = preprocess_task_arena(param)
%
% Strips data out of Arena packets and creates:
% 1. Files of just radar data
% 2. Files with just header data
%   a. binary flat file .hdr
%   b. Matlab .mat format (if supported)
% NOTE: This script only works on files with fixed header lengths. More
% specifically, every header must have the payload length field in the same spot.
%
% Example:
% Called from run_arena_packet_strip.m
%
% Author: John Paden
%
% See also: basic_load_arena.m, run_arena_packet_strip.m,
% arena_packet_strip.m, arena_packet_strip_task.m

% Pull in the inputs from param struct
base_dir = fullfile(param.config.base_dir);
config_folder_name = fullfile(param.config.config_folder_name);
reuse_tmp_files = param.config.reuse_tmp_files;
mat_or_bin_hdr_output = param.config.mat_or_bin_hdr_output;

%% Read each config/system XML file pair into a configs structure
% =========================================================================
% config_fns = get_filenames(fullfile(base_dir,config_folder_name),'','','config.xml',struct('recursive',true));
config_fns_dir = fullfile(base_dir,config_folder_name);
config_fns = get_filenames(config_fns_dir,'','','config.xml');

if isempty(config_fns)
  error('No configuration files found in %s.', config_fns_dir);
end

clear configs;
bad_mask = false(size(config_fns)); % Keep track of bad config files
for config_idx = 1:length(config_fns)
  config_fn = config_fns{config_idx};

  try
    configs(config_idx) = read_arena_xml(config_fn,'',param.config.board_map,param.config.tx_map);
  catch ME
    % If the error is one that can be ignored, just mark this config file
    % as bad and continue. Otherwise rethrow the error.
    if strcmpi(ME.identifier,'READ_ARENA_XML:MISSING_SYSTEM_XML')
      warning('%s\n', ME.getReport);
      bad_mask(config_idx) = true;
    else
      rethrow(ME)
    end
  end
end
% Remove bad config files
config_fns = config_fns(~bad_mask);
configs = configs(~bad_mask);

%% Process each segment of data
% =========================================================================
% With each config XML file:
% 1. Read in the corresponding system XML file and combine with the config
%   structure (already done above).
% 2. Associate data files with the system structure based on the timestamps
% 3. Packet strip using the system/config XML information
% 4. Create segment information
for config_idx = 1:length(configs)
  %% Initialize variables
  arena_radar_header_type; % Load radar header types
  
  last_bytes_m = [];
  last_bytes = zeros(64,1,'uint8');
  last_bytes_len = int32(0);
  num_expected = int32(-1);
  pkt_counter = int32(-1);
  if strcmpi(configs(config_idx).radar_name,'KUSnow') ...
      || strcmpi(configs(config_idx).radar_name,'SnowRadar2')
    radar_header_type = snow_radar_header_type; % arena_radar_header_type
    min_num_expected = int32(0);
    max_num_expected = int32(configs(config_idx).max_num_bins);
    default_num_expected = int32(512);
    num_header_fields = int32(9);
    length_field_offset = int32(68);
  elseif strcmpi(configs(config_idx).radar_name,'TOHFSounder')
    radar_header_type = hf_sounder_radar_header_type; % arena_radar_header_type
    min_num_expected = int32(0);
    max_num_expected = int32(configs(config_idx).max_num_bins);
    default_num_expected = int32(512);
    num_header_fields = int32(9);
    length_field_offset = int32(68);
  elseif strcmpi(configs(config_idx).radar_name,'DopplerScat')
    min_num_expected = int32(0);
    max_num_expected = int32(configs(config_idx).max_num_bins);
    default_num_expected = int32(512);
    num_header_fields = int32(33);
    length_field_offset = int32(260);
  elseif strcmpi(configs(config_idx).radar_name,'ku0001') || strcmpi(configs(config_idx).radar_name,'ku0002')
    radar_header_type = ku0001_radar_header_type; % arena_radar_header_type
    min_num_expected = int32(0);
    max_num_expected = int32(configs(config_idx).max_num_bins);
    default_num_expected = int32(512);
    num_header_fields = int32(9);
    length_field_offset = int32(68);
    
  else
    warning('Not a supported arena config XML radar name %s. Could indicate a file error. If a new system, then an entry may need to be added here. Run "dbcont" to continue and use default values.', configs(config_idx).radar_name);
    keyboard
    min_num_expected = int32(0);
    max_num_expected = int32(configs(config_idx).max_num_bins);
    default_num_expected = int32(512);
    num_header_fields = int32(9);
    length_field_offset = int32(68);
  end
  
  %% Get Files for each board
  for board_idx = 1:length(param.config.board_map)
    board = param.config.board_map(board_idx);
    
    board_folder_name = fullfile(param.config.board_folder_name);
    
    % Replace all "%b" in board_folder_name with the board number
    board_folder_name = regexprep(board_folder_name,'%b',board);
    
    % =========================================================================
    %% Get Data File list for this board
    fns = get_filenames(fullfile(base_dir,board_folder_name),'','','.dat',struct('recursive',true,'regexp','[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]'));
    fns_datenum = zeros(size(fns));
    for fn_idx = 1:length(fns)
      config_fname_info = fname_info_arena(fns{fn_idx});
      fns_datenum(fn_idx) = config_fname_info.datenum;
    end
    
    fns_mask = configs(config_idx).config_fname_info.datenum == fns_datenum;
    configs(config_idx).fns{board_idx} = fns(fns_mask);
    
    
    %% Iterate packet_strip through file list
    old_fn_dir = [];
    for fn_idx = 1:length(configs(config_idx).fns{board_idx})
      fn = fullfile(configs(config_idx).fns{board_idx}{fn_idx});
      
      [fn_dir,fn_name] = fileparts(fn);
      if ~strcmpi(fn_dir,old_fn_dir)
        % New data directory: assume that this is from a different Arena 313
        % board and state vectors should be reset.
        last_bytes_m = [];
        last_bytes = zeros(64,1,'uint8');
        last_bytes_len = int32(0);
        num_expected = int32(-1);
        pkt_counter = int32(-1);
        old_fn_dir = fn_dir;
      end
      
      % Create output filenames
      out_fn = ct_filename_ct_tmp(param,'','headers', ...
        fullfile(board_folder_name, fn_name));
      [out_fn_dir,out_fn_name] = fileparts(out_fn);
      out_fn = fullfile(out_fn_dir,[out_fn_name,'.dat']);
      if strcmpi(mat_or_bin_hdr_output,'.mat')
        out_hdr_fn = fullfile(out_fn_dir,[out_fn_name,'.mat']);
      else
        out_hdr_fn = fullfile(out_fn_dir,[out_fn_name,'.hdr']);
      end
      
      % Print status
      fprintf('arena_pkt_strip %d/%d %d/%d %s (%s)\n    %s\n', config_idx, ...
        length(configs), fn_idx, length(configs(config_idx).fns{board_idx}), fn, datestr(now), out_fn);
      
      % Check to make sure output directory exists
      if ~exist(out_fn_dir,'dir')
        mkdir(out_fn_dir);
      end
      
      % Copy XML files
      if fn_idx == 1 && board_idx == 1
        [~,out_config_fn_name,out_config_fn_ext] = fileparts(configs(config_idx).config_fn);
        out_config_fn = ct_filename_ct_tmp(param,'','headers', ...
          fullfile(config_folder_name, [out_config_fn_name '.xml']));
        fprintf('Copy %s\n  %s\n', configs(config_idx).config_fn,out_config_fn);
        try
          copyfile(configs(config_idx).config_fn,out_config_fn);
          if ispc
            fileattrib(out_config_fn,'+w');
          else
            fileattrib(out_config_fn,'+w -x','ug');
          end
        catch ME
          warning(ME.getReport)
        end
        [~,out_config_fn_name,out_config_fn_ext] = fileparts(configs(config_idx).system_fn);
        out_config_fn = ct_filename_ct_tmp(param,'','headers', ...
          fullfile(config_folder_name, [out_config_fn_name '.xml']));
        fprintf('Copy %s\n  %s\n', configs(config_idx).system_fn,out_config_fn);
        try
          copyfile(configs(config_idx).system_fn,out_config_fn);
          if ispc
            fileattrib(out_config_fn,'+w');
          else
            fileattrib(out_config_fn,'+w -x','ug');
          end
        catch ME
          warning(ME.getReport)
        end
      end
      
      % Copy Log Files
      if fn_idx == 1 && board_idx == 1
        [out_config_fn_dir] = fileparts(configs(config_idx).config_fn);
        log_files = fullfile(out_config_fn_dir,'logs/*');
        out_log_dir = fullfile(param.data_support_path, param.season_name, param.config.config_folder_name);
        try
          fprintf('Copy %s\n  %s\n', log_files, out_log_dir);
          if ~exist(out_log_dir,'dir')
            mkdir(out_log_dir)
          end
          copyfile(log_files, out_log_dir);
        catch ME
          warning('Error while copying log files:\n%s\n', ME.getReport);
        end
      end
      
      % Check to see if outputs already exist
      if reuse_tmp_files && exist(out_hdr_fn,'file')
        % For older systems that used UDP (datastream_type == 'udp'), the
        % UDP packet headers are in the raw data files and preprocess also
        % creates a copy of the data files without the packet header. Look
        % for this file in that case. Newer systems can override this check
        % by setting param.config.arena.daq.udp_packet_headers = false.
        if ~strcmpi(configs(config_idx).datastream_type,'udp') ...
            || exist(out_fn,'file') ...
            || isfield(param.config.arena,'daq') ...
               && isfield(param.config.arena.daq,'udp_packet_headers') ...
               && ~param.config.arena.daq.udp_packet_headers
          % Load the file to ensure it is not corrupted:
          % * If not corrupted, then execution continues onto the next
          %   file.
          % * If the file is corrupted, matlab will produce an error here
          %   and stop. The corrupt file should be deleted and the
          %   preprocess run again.
          load(out_hdr_fn,'last_bytes','last_bytes_len','num_expected','pkt_counter');
          continue;
        end
      end
      
      %% Read in headers from data file and create network packet stripped data file
      if strcmpi(configs(config_idx).datastream_type,'udp') ...
          || ~isfield(param.config.arena,'daq') ...
          || ~isfield(param.config.arena.daq,'udp_packet_headers') ...
          || param.config.arena.daq.udp_packet_headers
        % In old arena systems, choosing the UDP datastream created raw
        % data files with UDP packet headers in them and a copy of the raw
        % data file without the packet headers will be made in "out_fn". If
        % using a new system, then the UDP packet headers are not in the
        % files, but param.config.arena.daq.udp_packet_headers must be set
        % to true for preprocess to process them correctly.
        [hdr,last_bytes_len,num_expected,pkt_counter] = arena_packet_strip_mex(fn,out_fn,last_bytes,last_bytes_len, ...
          num_expected,pkt_counter,min_num_expected,max_num_expected, ...
          default_num_expected,num_header_fields,length_field_offset);
      else
        % The "out_fn" input argument is ignored. These raw data files do
        % not contain UDP packet headers.
        [hdr,last_bytes_len,num_expected,pkt_counter] = arena_packet_strip_tcp_mex(fn,out_fn,last_bytes,last_bytes_len, ...
          num_expected,pkt_counter,min_num_expected,max_num_expected, ...
          default_num_expected,num_header_fields,length_field_offset);
      end
      
      %% Write header output file
      if strcmpi(mat_or_bin_hdr_output,'.mat')
        if strcmpi(configs(config_idx).radar_name,'ku0001')
          offset = mod(hdr(1,:),2^32);
          mode_latch = mod(hdr(3,:),2^8);
          subchannel = mod(bitshift(hdr(3,:),-8),2^8);
          wg_delay_latch = mod(hdr(4,:),2^16);
          rel_time_cntr_latch = double(hdr(5,:));
          profile_cntr_latch = double(hdr(6,:));
          pps_ftime_cntr_latch = double(hdr(7,:));
          pps_cntr_latch = double(hdr(8,:));
          
          save(out_hdr_fn, 'offset','mode_latch','subchannel','wg_delay_latch', ...
            'rel_time_cntr_latch','profile_cntr_latch','pps_ftime_cntr_latch','pps_cntr_latch', ...
            'last_bytes','last_bytes_len','num_expected','pkt_counter');
          
        elseif strcmpi(configs(config_idx).radar_name,'ku0002') 
          offset = hdr{1}.frame_sync;
          mode_latch = hdr{1}.mode;
          subchannel = hdr{1}.subchannel;
          wg_delay_latch = 0;
          rel_time_cntr_latch = hdr{1}.rel_time_cntr_latch;
          profile_cntr_latch = hdr{1}.profile_cntr_latch;
          pps_ftime_cntr_latch =hdr{1}.pps_ftime_cntr_latch;
          pps_cntr_latch = hdr{1}.pps_cntr_latch;
          
          save(out_hdr_fn, 'offset','mode_latch','subchannel','wg_delay_latch', ...
            'rel_time_cntr_latch','profile_cntr_latch','pps_ftime_cntr_latch','pps_cntr_latch', ...
            'last_bytes','last_bytes_len','num_expected','pkt_counter');
          
        elseif strcmpi(configs(config_idx).radar_name,'KUSnow') ...
          || strcmpi(configs(config_idx).radar_name,'SnowRadar2')
          offset = mod(hdr(1,:),2^32);
          mode_latch = mod(hdr(3,:),2^8);
          subchannel = mod(bitshift(hdr(3,:),-8),2^8);
          wg_delay_latch = mod(hdr(4,:),2^16);
          rel_time_cntr_latch = double(hdr(5,:));
          profile_cntr_latch = double(hdr(6,:));
          pps_ftime_cntr_latch = double(hdr(7,:));
          pps_cntr_latch = double(hdr(8,:));
          
          save(out_hdr_fn, 'offset','mode_latch','subchannel','wg_delay_latch', ...
            'rel_time_cntr_latch','profile_cntr_latch','pps_ftime_cntr_latch','pps_cntr_latch', ...
            'last_bytes','last_bytes_len','num_expected','pkt_counter');
          
        elseif strcmpi(configs(config_idx).radar_name,'TOHFSounder')
          offset = mod(hdr(1,:),2^32);
          mode_latch = mod(hdr(3,:),2^8);
          subchannel = mod(bitshift(hdr(3,:),-8),2^8);
          encoder = mod(hdr(4,:),2^32);
          rel_time_cntr_latch = double(hdr(5,:));
          profile_cntr_latch = double(hdr(6,:));
          pps_ftime_cntr_latch = double(hdr(7,:));
          pps_cntr_latch = double(hdr(8,:));
          
          save(out_hdr_fn, 'offset','mode_latch','subchannel','encoder', ...
            'rel_time_cntr_latch','profile_cntr_latch','pps_ftime_cntr_latch','pps_cntr_latch', ...
            'last_bytes','last_bytes_len','num_expected','pkt_counter');
          
        elseif strcmpi(configs(config_idx).radar_name,'DopplerScat')
          offset = mod(hdr(1,:),2^32);
          mode_latch = mod(hdr(3,:),2^8);
          decimation_ratio = mod(bitshift(hdr(3,:),-8),2^8);
          num_pulses_burst = mod(bitshift(hdr(3,:),-16),2^8);
          encoder = mod(hdr(4,:),2^32);
          rel_time_cntr_latch = double(hdr(5,:));
          profile_cntr_latch = double(hdr(6,:));
          pps_ftime_cntr_latch = double(hdr(7,:));
          pps_cntr_latch = double(hdr(8,:));
          
          save(out_hdr_fn, 'offset','mode_latch','decimation_ratio','num_pulses_burst','encoder', ...
            'rel_time_cntr_latch','profile_cntr_latch','pps_ftime_cntr_latch','pps_cntr_latch', ...
            'last_bytes','last_bytes_len','num_expected','pkt_counter');
        end
        
        if 0
          %% Debug outputs
          % load(out_tmp_fn);
          plot(subchannel,'.');
          sum(subchannel)*2
          length(subchannel)
          plot(profile_cntr_latch(subchannel==0))
          plot(diff(profile_cntr_latch(subchannel==0)),'.')
          plot(diff(profile_cntr_latch(subchannel==1)),'.')
          any(diff(profile_cntr_latch(subchannel==0)) ~= 2)
          any(diff(profile_cntr_latch(subchannel==1)) ~= 2)
          plot(mode_latch(subchannel==0))
          plot(wg_delay_latch(subchannel==0))
          plot(rel_time_cntr_latch(subchannel==0))
          plot(diff(rel_time_cntr_latch(subchannel==0)))
          plot(diff(rel_time_cntr_latch(subchannel==1)))
          plot(profile_cntr_latch(subchannel==0))
          plot(diff(profile_cntr_latch(subchannel==0)))
          plot(diff(profile_cntr_latch(subchannel==1)))
          time = pps_cntr_latch + pps_ftime_cntr_latch/10e6;
          plot(time(subchannel==0))
          plot(diff(time(subchannel==0)))
          plot(time(subchannel==1))
          plot(diff(time(subchannel==1)))
        end
        
      else
        out_hdr_fid = fopen(out_hdr_fn,'w');
        fwrite(out_hdr_fid,hdr);
        fclose(out_hdr_fid);
      end
      
    end
  end
end

%% Print out segments
% =========================================================================
oparams = {};
for config_idx = 1:length(configs)
  
  if isempty(configs(config_idx).fns{1})
    % Skip this config file if there are no data files
    continue;
  end
  
  % Determine which default parameters to use
  % =======================================================================
  match_idx = [];
  for default_idx = 1:length(param.config.defaults)
    match = regexpi(configs(config_idx).psc.config_name, param.config.defaults{default_idx}.config_regexp);
    if ~isempty(match)
      match_idx = default_idx;
      break;
    end
  end
  if isempty(match_idx)
    error('No match for psc config name %s.', configs(config_idx).psc.config_name);
  end
  oparams{end+1} = param.config.defaults{match_idx};

  % Create map from wfs to board_idx, mode, subchannel, adc
  % =======================================================================
  data_map = oparams{end}.records.data_map;
  board_idx_map = [];
  mode_map = [];
  subchannel_map = [];
  wfs_map = [];
  adc_map = [];
  for board_idx = 1:length(data_map)
    for data_idx = 1:size(data_map{board_idx},1)
      board_idx_map(end+1) = board_idx;
      mode_map(end+1) = data_map{board_idx}(data_idx,1);
      subchannel_map(end+1) = data_map{board_idx}(data_idx,2);
      wfs_map(end+1) = data_map{board_idx}(data_idx,3);
      adc_map(end+1) = data_map{board_idx}(data_idx,4);
    end
  end
  [wfs,unique_map] = unique(wfs_map);
  board_idx_wfs = board_idx_map(unique_map);
  mode_wfs = mode_map(unique_map);
  subchannel_wfs = subchannel_map(unique_map);
  adc_wfs = adc_map(unique_map);
  
  % Parameter spreadsheet
  % =======================================================================
  [~,config_fn_name] = fileparts(configs(config_idx).config_fn);
  oparams{end}.day_seg = sprintf('%s_%02d',config_fn_name(1:8),numel(oparams));
  oparams{end}.cmd.notes = configs(config_idx).psc.config_name(5:end);
  
  oparams{end}.records.file.version = 103;
  oparams{end}.records.file.boards = param.config.board_map;
  oparams{end}.records.file.prefix = datestr(configs(config_idx).config_fname_info.datenum,'YYYYmmDD_HHMMSS');
  for board_idx = 1:length(param.config.board_map)
    oparams{end}.records.file.start_idx(board_idx) = 1;
    oparams{end}.records.file.stop_idx(board_idx) = length(configs(config_idx).fns{board_idx});
  end
  if strcmpi(configs(config_idx).datastream_type,'udp') ...
      || ~isfield(param.config.arena,'daq') ...
      || ~isfield(param.config.arena.daq,'udp_packet_headers') ...
      || param.config.arena.daq.udp_packet_headers
    % See earlier discussion on udp_packet_headers. The raw data files with
    % the UDP packet headers removed are stored in ct_tmp.
    oparams{end}.records.file.base_dir = ct_filename_ct_tmp(param,'','headers','');
  else
    oparams{end}.records.file.base_dir = base_dir;
  end
  oparams{end}.records.file.board_folder_name = param.config.board_folder_name;
  if ~isnan(str2double(oparams{end}.records.file.board_folder_name))
    oparams{end}.records.file.board_folder_name = ['/' oparams{end}.records.file.board_folder_name];
  end
  oparams{end}.records.file.clk = 10e6;
  oparams{end}.records.gps.time_offset = 0;
  oparams{end}.records.gps.en = 1;
  [~,config_fn_name] = fileparts(configs(config_idx).config_fn);
  oparams{end}.records.config_fn = fullfile(param.config.config_folder_name, [config_fn_name '.xml']);
  
  oparams{end}.radar.fs = configs(config_idx).adc{board_idx_wfs(1),1+mode_wfs(1),1+subchannel_wfs(1)}.sampFreq;
  oparams{end}.radar.prf = configs(config_idx).prf * configs(config_idx).total_presums;
  
  % Usually the default.radar.wfs structure only has one waveform
  % entry which is to be copied to all the waveforms.
  if numel(oparams{end}.radar.wfs) == 1
    oparams{end}.radar.wfs = repmat(oparams{end}.radar.wfs,[1 numel(wfs)]);
  end
  
  for wf_idx = 1:numel(wfs)
    wf = wfs(wf_idx);
    
    board_idx = board_idx_wfs(wf_idx);
    mode_latch = mode_wfs(wf_idx);
    subchannel = subchannel_wfs(wf_idx);
    if strcmpi(param.season_name,'2022_Greenland_P3')
      % Delete This After AITT 2022 done
      configs(config_idx).dac{1,mode_latch+1}.wfs{1}.centerFreq = 10e9/1e6; % TEMP FIX
      configs(config_idx).dac{1,mode_latch+1}.wfs{1}.bandwidth = 15e9/1e6; % TEMP FIX
      configs(config_idx).dac{1,mode_latch+1}.sampFreq = 2400; % TEMP FIX
      configs(config_idx).dac{1,mode_latch+1}.wfs{1}.numPoints = 240e-6*configs(config_idx).dac{1,mode_latch+1}.sampFreq*1e6; % TEMP FIX
      configs(config_idx).dac{1}.delay = 0; % TEMP FIX
      configs(config_idx).dac{1,mode_latch+1}.wfs{1}.alpha = 15e9/240e-6;
      for tx = 1:length(oparams{end}.radar.wfs(wf).tx_paths)
        tx_idx = oparams{end}.radar.wfs(wf).tx_paths(tx);
        if isfinite(oparams{end}.radar.wfs(wf).tx_paths(tx))
          configs(config_idx).dac{tx_idx,mode_latch+1}.wfs{1}.scale = 1;
        else
          scale(tx) = 0;
        end
      end
    end
    
    fc = configs(config_idx).dac{1,mode_latch+1}.wfs{1}.centerFreq*1e6;
    BW = configs(config_idx).dac{1,mode_latch+1}.wfs{1}.bandwidth*1e6;
    if strcmpi(configs(config_idx).radar_name,'ku0002')
      fc = fc*param.config.defaults{1}.radar.wfs(wf).fmult + param.config.defaults{1}.radar.wfs(wf).fLO;
      BW = BW*param.config.defaults{1}.radar.wfs(wfs).fmult;
    end
    Nt = configs(config_idx).dac{1,mode_latch+1}.wfs{1}.numPoints;
    fs = configs(config_idx).dac{1,mode_latch+1}.sampFreq*1e6;
    Tpd = Nt/fs;
    t_dac = (configs(config_idx).dac{1}.delay) * 1e-6;
    
    switch (configs(config_idx).adc{board_idx_wfs(1),1+mode_wfs(1),1+subchannel_wfs(1)}.adcMode)
      case 0
        oparams{end}.radar.wfs(wf).DDC_dec = 1;
      case 1
        oparams{end}.radar.wfs(wf).DDC_dec = 2;
      case 2
        oparams{end}.radar.wfs(wf).DDC_dec = 4;
    end
    
    if subchannel_wfs(1) == 0
      switch (configs(config_idx).adc{board_idx_wfs(1),1+mode_wfs(1),1+subchannel_wfs(1)}.ddc0NcoMode)
        case 0
          oparams{end}.radar.wfs(wf).DDC_freq = 0;
        case 1
          oparams{end}.radar.wfs(wf).DDC_freq = oparams{end}.radar.fs/4;
        case 2
          oparams{end}.radar.wfs(wf).DDC_freq = configs(config_idx).adc{board_idx_wfs(1),1+mode_wfs(1),1+subchannel_wfs(1)}.ddc0NcoFreq;
      end
    elseif subchannel_wfs(1) == 1
      switch (configs(config_idx).adc{board_idx_wfs(1),1+mode_wfs(1),1+subchannel_wfs(1)}.ddc1NcoMode)
        case 0
          oparams{end}.radar.wfs(wf).DDC_freq = 0;
        case 1
          oparams{end}.radar.wfs(wf).DDC_freq = oparams{end}.radar.fs/4;
        case 2
          oparams{end}.radar.wfs(wf).DDC_freq = configs(config_idx).adc{board_idx_wfs(1),1+mode_wfs(1),1+subchannel_wfs(1)}.ddc1NcoFreq;
      end
    end
    
    oparams{end}.radar.wfs(wf).f0 = fc-BW/2;
    oparams{end}.radar.wfs(wf).f1 = fc+BW/2;
    oparams{end}.radar.wfs(wf).tukey = configs(config_idx).dac{1,mode_latch+1}.wfs{1}.alpha;
    oparams{end}.radar.wfs(wf).BW_window = [fc-BW/2 fc+BW/2];
    oparams{end}.radar.wfs(wf).Tpd = Tpd;
    scale = [];
    for tx = 1:length(oparams{end}.radar.wfs(wf).tx_paths)
      if isfinite(oparams{end}.radar.wfs(wf).tx_paths(tx))
        tx_idx = oparams{end}.radar.wfs(wf).tx_paths(tx);
        scale(tx) = configs(config_idx).dac{tx_idx,mode_latch+1}.wfs{1}.scale;
      else
        scale(tx) = 0;
      end
    end
    oparams{end}.radar.wfs(wf).tx_weights = scale;
    oparams{end}.radar.wfs(wf).presums = configs(config_idx).adc{board_idx,mode_latch+1,subchannel+1}.presums;
    adc_idxs = find(wfs_map == wf);
    for adc_idx = adc_idxs
      adc = adc_map(adc_idx);
      adc_board_idx = board_idx_map(adc_idx);
      adc_mode = mode_map(adc_idx);
      adc_subchannel = subchannel_map(adc_idx);
      oparams{end}.radar.wfs(wf).bit_shifts(adc) = configs(config_idx).adc{adc_board_idx,adc_mode+1,adc_subchannel+1}.shiftLSB - 2 - param.config.arena.adc(adc_board_idx).gain_dB(adc_subchannel+1)/6;
    end
    oparams{end}.radar.wfs(wf).Tadc = sscanf(configs(config_idx).adc{board_idx,mode_latch+1,subchannel+1}.rg,'%d') ...
      / oparams{end}.radar.fs*oparams{end}.radar.wfs(wf).DDC_dec ...
      - param.config.arena.param.ADC_time_delay - t_dac;
    
  end
end

%% Print out segments
% =========================================================================
if ~isempty(param.config.param_fn)
  % Print parameter spreadsheet values to stdout and param_txt_fn
  % =========================================================================
  fid = 1;
  fprintf(fid,'<strong>%s\n','='*ones(1,80)); fprintf(fid,'  cmd\n'); fprintf(fid,'%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'cmd',oparams,fid);
  fprintf(fid,'<strong>%s\n','='*ones(1,80)); fprintf(fid,'  records\n'); fprintf(fid,'%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'records',oparams,fid);
  fprintf(fid,'<strong>%s\n','='*ones(1,80)); fprintf(fid,'  qlook\n'); fprintf(fid,'%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'qlook',oparams,fid);
  fprintf(fid,'<strong>%s\n','='*ones(1,80)); fprintf(fid,'  sar\n'); fprintf(fid,'%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'sar',oparams,fid);
  fprintf(fid,'<strong>%s\n','='*ones(1,80)); fprintf(fid,'  array\n'); fprintf(fid,'%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'array',oparams,fid);
  fprintf(fid,'<strong>%s\n','='*ones(1,80)); fprintf(fid,'  radar\n'); fprintf(fid,'%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'radar',oparams,fid);
  fprintf(fid,'<strong>%s\n','='*ones(1,80)); fprintf(fid,'  post\n'); fprintf(fid,'%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'post',oparams,fid);
  fprintf(fid,'<strong>%s\n','='*ones(1,80)); fprintf(fid,'  analysis\n'); fprintf(fid,'%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'analysis',oparams,fid);
  fprintf(fid,'\n');
  
  param_txt_fn = ct_filename_ct_tmp(param,'','param', [param.config.date_str,'.txt']);
  fprintf('Writing %s\n\n', param_txt_fn);
  param_txt_fn_dir = fileparts(param_txt_fn);
  if ~exist(param_txt_fn_dir,'dir')
    mkdir(param_txt_fn_dir);
  end
  [fid,msg] = fopen(param_txt_fn,'wb');
  if fid<0
    error('Could not write to %s: %s\n', param_txt_fn, msg);
  end
  fprintf(fid,'%s\n','='*ones(1,80)); fprintf(fid,'  cmd\n'); fprintf(fid,'%s\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'cmd',oparams,fid);
  fprintf(fid,'\n');
  fprintf(fid,'%s\n','='*ones(1,80)); fprintf(fid,'  records\n'); fprintf(fid,'%s\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'records',oparams,fid);
  fprintf(fid,'\n');
  fprintf(fid,'%s\n','='*ones(1,80)); fprintf(fid,'  qlook\n'); fprintf(fid,'%s\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'qlook',oparams,fid);
  fprintf(fid,'\n');
  fprintf(fid,'%s\n','='*ones(1,80)); fprintf(fid,'  sar\n'); fprintf(fid,'%s\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'sar',oparams,fid);
  fprintf(fid,'\n');
  fprintf(fid,'%s\n','='*ones(1,80)); fprintf(fid,'  array\n'); fprintf(fid,'%s\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'array',oparams,fid);
  fprintf(fid,'\n');
  fprintf(fid,'%s\n','='*ones(1,80)); fprintf(fid,'  radar\n'); fprintf(fid,'%s\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'radar',oparams,fid);
  fprintf(fid,'\n');
  fprintf(fid,'%s\n','='*ones(1,80)); fprintf(fid,'  post\n'); fprintf(fid,'%s\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'post',oparams,fid);
  fprintf(fid,'\n');
  fprintf(fid,'%s\n','='*ones(1,80)); fprintf(fid,'  analysis\n'); fprintf(fid,'%s\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'analysis',oparams,fid);
  fprintf(fid,'\n');
  fclose(fid);
end

%% Exit task
% =========================================================================
fprintf('%s done %s\n', mfilename, datestr(now));

success = true;
