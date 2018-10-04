function success = arena_packet_strip_task(param)
% success = arena_packet_strip_task(param)
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
base_dir = fullfile(param.arena_packet_strip.base_dir);
config_folder_name = fullfile(param.arena_packet_strip.config_folder_name);
reuse_tmp_files = param.arena_packet_strip.reuse_tmp_files;
mat_or_bin_hdr_output = param.arena_packet_strip.mat_or_bin_hdr_output;

%% Read each config/system XML file pair into a configs structure
% =========================================================================
% config_fns = get_filenames(fullfile(base_dir,config_folder_name),'','','config.xml',struct('recursive',true));
config_fns = get_filenames(fullfile(base_dir,config_folder_name),'','','config.xml');

clear configs;
for config_idx = 1:length(config_fns)
  config_fn = config_fns{config_idx};
  
  configs(config_idx) = read_arena_xml(config_fn,'',param.arena_packet_strip.board_map,param.arena_packet_strip.tx_map);
end

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
  if strcmpi(configs(config_idx).radar_name,'KUSnow')
    radar_header_type = snow_radar_header_type;
    min_num_expected = int32(0);
    max_num_expected = int32(configs(config_idx).max_num_bins);
    default_num_expected = int32(512);
    num_header_fields = int32(9);
    length_field_offset = int32(68);
  elseif strcmpi(configs(config_idx).radar_name,'TOHFSounder')
    radar_header_type = hf_sounder_radar_header_type;
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
  elseif strcmpi(configs(config_idx).radar_name,'ku0001')
    radar_header_type = ku0001_radar_header_type;
    min_num_expected = int32(0);
    max_num_expected = int32(configs(config_idx).max_num_bins);
    default_num_expected = int32(512);
    num_header_fields = int32(9);
    length_field_offset = int32(68);
    
  else
    error('Not a supported radar header type %d.', radar_header_type);
  end
  
  %% Get Files for each board
  for board_idx = 1:length(param.arena_packet_strip.board_map)
    board = param.arena_packet_strip.board_map(board_idx);
    
    board_folder_name = fullfile(param.arena_packet_strip.board_folder_name);
    
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
        copyfile(configs(config_idx).config_fn,out_config_fn);
        if ispc
          fileattrib(out_config_fn,'+w');
        else
          fileattrib(out_config_fn,'+w -x');
        end
        [~,out_config_fn_name,out_config_fn_ext] = fileparts(configs(config_idx).system_fn);
        out_config_fn = ct_filename_ct_tmp(param,'','headers', ...
          fullfile(config_folder_name, [out_config_fn_name '.xml']));
        fprintf('Copy %s\n  %s\n', configs(config_idx).system_fn,out_config_fn);
        copyfile(configs(config_idx).system_fn,out_config_fn);
        if ispc
          fileattrib(out_config_fn,'+w');
        else
          fileattrib(out_config_fn,'+w -x');
        end
      end
      
      % Copy Log Files
      if fn_idx == 1 && board_idx == 1
        [out_config_fn_dir] = fileparts(configs(config_idx).config_fn);
        log_files = fullfile(out_config_fn_dir,'logs/*');
        global gRadar
        out_log_dir = fullfile(gRadar.data_support_path, param.season_name, param.arena_packet_strip.config_folder_name);
        if ~exist(out_log_dir,'dir')
          mkdir(out_log_dir)
        end
        fprintf('Copy %s\n  %s\n', log_files, out_log_dir);
        copyfile(log_files, out_log_dir);
      end
      
      % Check to see if outputs already exist
      if reuse_tmp_files && exist(out_fn,'file') && exist(out_hdr_fn,'file')
        continue;
      end
      
      %% Read in headers from data file and create network packet stripped data file
      [hdr,last_bytes_len,num_expected,pkt_counter] = arena_packet_strip_mex(fn,out_fn,last_bytes,last_bytes_len, ...
        num_expected,pkt_counter,min_num_expected,max_num_expected, ...
        default_num_expected,num_header_fields,length_field_offset);
      
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
            'rel_time_cntr_latch','profile_cntr_latch','pps_ftime_cntr_latch','pps_cntr_latch');
          
        elseif strcmpi(configs(config_idx).radar_name,'KUSnow')
          offset = mod(hdr(1,:),2^32);
          mode_latch = mod(hdr(3,:),2^8);
          subchannel = mod(bitshift(hdr(3,:),-8),2^8);
          wg_delay_latch = mod(hdr(4,:),2^16);
          rel_time_cntr_latch = double(hdr(5,:));
          profile_cntr_latch = double(hdr(6,:));
          pps_ftime_cntr_latch = double(hdr(7,:));
          pps_cntr_latch = double(hdr(8,:));
          
          save(out_hdr_fn, 'offset','mode_latch','subchannel','wg_delay_latch', ...
            'rel_time_cntr_latch','profile_cntr_latch','pps_ftime_cntr_latch','pps_cntr_latch');
          
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
            'rel_time_cntr_latch','profile_cntr_latch','pps_ftime_cntr_latch','pps_cntr_latch');
          
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
            'rel_time_cntr_latch','profile_cntr_latch','pps_ftime_cntr_latch','pps_cntr_latch');
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
for config_idx = 1:length(configs)
  
  % Determine which default parameters to use
  % =======================================================================
  default_param = param.arena_packet_strip.default_param;
  defaults = param.arena_packet_strip.defaults;
  match_idx = [];
  for default_idx = 1:length(defaults)
    match = regexpi(configs(config_idx).psc.config_name, defaults{default_idx}.config_regexp);
    if ~isempty(match)
      match_idx = default_idx;
      break;
    end
  end
  if isempty(match_idx)
    error('No match for psc config name %s.', configs(config_idx).psc.config_name);
  end
  if config_idx == 1
    oparams = default_param;
  end

  % Create map from wfs to board_idx, mode, subchannel, adc
  % =======================================================================
  data_map = defaults{match_idx}.records.data_map;
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
  oparams(config_idx).day_seg = sprintf('%s_%02d',config_fn_name(1:8),config_idx);
  oparams(config_idx).cmd.notes = configs(config_idx).psc.config_name(5:end);
  oparams(config_idx).records = defaults{match_idx}.records;
  oparams(config_idx).qlook = defaults{match_idx}.qlook;
  oparams(config_idx).sar = defaults{match_idx}.sar;
  oparams(config_idx).array = defaults{match_idx}.array;
  
  oparams(config_idx).records.file.version = 103;
  oparams(config_idx).records.file.boards = param.arena_packet_strip.board_map;
  oparams(config_idx).records.file.prefix = datestr(configs(config_idx).config_fname_info.datenum,'YYYYmmDD_HHMMSS');
  for board_idx = 1:length(param.arena_packet_strip.board_map)
    oparams(config_idx).records.file.start_idx(board_idx) = 1;
    oparams(config_idx).records.file.stop_idx(board_idx) = length(configs(config_idx).fns{board_idx});
  end
  oparams(config_idx).records.file.base_dir = ct_filename_ct_tmp(param,'','headers','');
  oparams(config_idx).records.file.board_folder_name = param.arena_packet_strip.board_folder_name;
  if ~isnan(str2double(oparams(config_idx).records.file.board_folder_name))
    oparams(config_idx).records.file.board_folder_name = ['/' oparams(config_idx).records.file.board_folder_name];
  end
  oparams(config_idx).records.file.clk = 10e6;
  oparams(config_idx).records.gps.time_offset = 0;
  oparams(config_idx).records.gps.en = 1;
  [~,config_fn_name] = fileparts(configs(config_idx).config_fn);
  oparams(config_idx).records.config_fn = fullfile(param.arena_packet_strip.config_folder_name, [config_fn_name '.xml']);
  
  oparams(config_idx).radar.fs = configs(config_idx).adc{board_idx_wfs(1),1+mode_wfs(1),1+subchannel_wfs(1)}.sampFreq;
  oparams(config_idx).radar.prf = configs(config_idx).prf * configs(config_idx).total_presums;
  oparams(config_idx).radar.adc_bits = defaults{match_idx}.radar.adc_bits;
  oparams(config_idx).radar.Vpp_scale = defaults{match_idx}.radar.Vpp_scale;
  oparams(config_idx).radar.lever_arm_fh = defaults{match_idx}.radar.lever_arm_fh;
  
  for wf_idx = 1:length(wfs)
    wf = wfs(wf_idx);
    
    board_idx = board_idx_wfs(wf_idx);
    mode_latch = mode_wfs(wf_idx);
    subchannel = subchannel_wfs(wf_idx);
    
    fc = configs(config_idx).dac{1,mode_latch+1}.wfs{1}.centerFreq*1e6;
    BW = configs(config_idx).dac{1,mode_latch+1}.wfs{1}.bandwidth*1e6;
    Nt = configs(config_idx).dac{1,mode_latch+1}.wfs{1}.numPoints;
    fs = configs(config_idx).dac{1,mode_latch+1}.sampFreq*1e6;
    Tpd = Nt/fs;
    t_dac = (configs(config_idx).dac{1}.delay) * 1e-6;
    
    switch (configs(config_idx).adc{board_idx_wfs(1),1+mode_wfs(1),1+subchannel_wfs(1)}.adcMode)
      case 0
        oparams(config_idx).radar.wfs(wf).DDC_dec = 1;
      case 1
        oparams(config_idx).radar.wfs(wf).DDC_dec = 2;
      case 2
        oparams(config_idx).radar.wfs(wf).DDC_dec = 4;
    end
    
    if subchannel_wfs(1) == 0
      switch (configs(config_idx).adc{board_idx_wfs(1),1+mode_wfs(1),1+subchannel_wfs(1)}.ddc0NcoMode)
        case 0
          oparams(config_idx).radar.wfs(wf).DDC_freq = 0;
        case 1
          oparams(config_idx).radar.wfs(wf).DDC_freq = oparams(config_idx).radar.fs/4;
        case 2
          oparams(config_idx).radar.wfs(wf).DDC_freq = configs(config_idx).adc{board_idx_wfs(1),1+mode_wfs(1),1+subchannel_wfs(1)}.ddc0NcoFreq;
      end
    elseif subchannel_wfs(1) == 1
      switch (configs(config_idx).adc{board_idx_wfs(1),1+mode_wfs(1),1+subchannel_wfs(1)}.ddc1NcoMode)
        case 0
          oparams(config_idx).radar.wfs(wf).DDC_freq = 0;
        case 1
          oparams(config_idx).radar.wfs(wf).DDC_freq = oparams(config_idx).radar.fs/4;
        case 2
          oparams(config_idx).radar.wfs(wf).DDC_freq = configs(config_idx).adc{board_idx_wfs(1),1+mode_wfs(1),1+subchannel_wfs(1)}.ddc1NcoFreq;
      end
    end
    
    oparams(config_idx).radar.wfs(wf).adcs = defaults{match_idx}.radar.wfs(wf).adcs;
    oparams(config_idx).radar.wfs(wf).f0 = fc-BW/2;
    oparams(config_idx).radar.wfs(wf).f1 = fc+BW/2;
    oparams(config_idx).radar.wfs(wf).tukey = configs(config_idx).dac{1,mode_latch+1}.wfs{1}.alpha;
    oparams(config_idx).radar.wfs(wf).BW_window = [fc-BW/2 fc+BW/2];
    oparams(config_idx).radar.wfs(wf).Tpd = Tpd;
    scale = [];
    for tx_idx = 1:size(configs(config_idx).dac,1)
      scale(tx_idx) = configs(config_idx).dac{tx_idx,mode_latch+1}.wfs{1}.scale;
    end
    oparams(config_idx).radar.wfs(wf).tx_weights = scale;
    oparams(config_idx).radar.wfs(wf).rx_paths = defaults{match_idx}.radar.wfs(wf).rx_paths;
    oparams(config_idx).radar.wfs(wf).adc_gains_dB = round(defaults{match_idx}.radar.wfs(wf).adc_gains_dB*10)/10;
    oparams(config_idx).radar.wfs(wf).chan_equal_dB = round(defaults{match_idx}.radar.wfs(wf).chan_equal_dB*10)/10;
    oparams(config_idx).radar.wfs(wf).chan_equal_deg = round(defaults{match_idx}.radar.wfs(wf).chan_equal_deg*10)/10;
    oparams(config_idx).radar.wfs(wf).Tsys = defaults{match_idx}.radar.wfs(wf).chan_equal_Tsys;
    oparams(config_idx).radar.wfs(wf).presums = configs(config_idx).adc{board_idx,mode_latch+1,subchannel+1}.presums;
    oparams(config_idx).radar.wfs(wf).bit_shifts = configs(config_idx).adc{board_idx,mode_latch+1,subchannel+1}.shiftLSB;
    oparams(config_idx).radar.wfs(wf).Tadc = sscanf(configs(config_idx).adc{board_idx,mode_latch+1,subchannel+1}.rg,'%d') ...
      / oparams(config_idx).radar.fs*oparams(config_idx).radar.wfs(wf).DDC_dec ...
      - param.arena_packet_strip.defaults{1}.arena.param.ADC_time_delay - t_dac;
    
  end
  oparams(config_idx).post = defaults{match_idx}.post;
end

if ~isempty(param.arena_packet_strip.param_fn)
  % Print parameter spreadsheet values
  % =========================================================================
  fprintf('<strong>%s\n','='*ones(1,80)); fprintf('  cmd\n'); fprintf('%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.arena_packet_strip.param_fn,'cmd',oparams);
  fprintf('<strong>%s\n','='*ones(1,80)); fprintf('  records\n'); fprintf('%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.arena_packet_strip.param_fn,'records',oparams);
  fprintf('<strong>%s\n','='*ones(1,80)); fprintf('  qlook\n'); fprintf('%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.arena_packet_strip.param_fn,'qlook',oparams);
  fprintf('<strong>%s\n','='*ones(1,80)); fprintf('  sar\n'); fprintf('%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.arena_packet_strip.param_fn,'sar',oparams);
  fprintf('<strong>%s\n','='*ones(1,80)); fprintf('  array\n'); fprintf('%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.arena_packet_strip.param_fn,'array',oparams);
  fprintf('<strong>%s\n','='*ones(1,80)); fprintf('  radar\n'); fprintf('%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.arena_packet_strip.param_fn,'radar',oparams);
  fprintf('<strong>%s\n','='*ones(1,80)); fprintf('  post\n'); fprintf('%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.arena_packet_strip.param_fn,'post',oparams);
end

%% Exit task
% =========================================================================
fprintf('%s done %s\n', mfilename, datestr(now));

success = true;
