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
% See also: basic_load_arena.m

% Pull in the inputs from param struct
base_dir = fullfile(param.arena_packet_strip.base_dir);
xml_folder_name = fullfile(param.arena_packet_strip.xml_folder_name);
reuse_tmp_files = param.arena_packet_strip.reuse_tmp_files;
mat_or_bin_hdr_output = param.arena_packet_strip.mat_or_bin_hdr_output;

% =========================================================================
%% Read each system XML file into a system structure
system_xml_fns = get_filenames(fullfile(base_dir,xml_folder_name),'','','system.xml',struct('recursive',true));

clear settings;
for xml_idx = 1:length(system_xml_fns)
  xml_fn = system_xml_fns{xml_idx};
  
  settings(xml_idx) = read_arena_xml(xml_fn);
end

% =========================================================================
%% Process each segment of data
for xml_idx = 1:length(settings)
  
  % With each system XML file:
  %   Read in the corresponding config XML file and combine with the system structure
  %   Associate data files with the system structure based on the timestamps
  %   Packet strip using the system/config XML information
  %   Create segment information
  
  %% Initialize variables
  arena_radar_header_type; % Load radar header types
  last_bytes_m = [];
  last_bytes = zeros(64,1,'uint8');
  last_bytes_len = int32(0);
  num_expected = int32(-1);
  pkt_counter = int32(-1);
  if strcmpi(settings(xml_idx).radar_name,'KUSnow')
    radar_header_type = snow_radar_header_type;
    min_num_expected = int32(0);
    max_num_expected = int32(settings(xml_idx).max_num_bins);
    default_num_expected = int32(512);
    num_header_fields = int32(9);
    length_field_offset = int32(68);
  elseif strcmpi(settings(xml_idx).radar_name,'TOHFSounder')
    radar_header_type = hf_sounder_radar_header_type;
    min_num_expected = int32(0);
    max_num_expected = int32(settings(xml_idx).max_num_bins);
    default_num_expected = int32(512);
    num_header_fields = int32(9);
    length_field_offset = int32(68);
  elseif strcmpi(settings(xml_idx).radar_name,'DopplerScat')
    min_num_expected = int32(0);
    max_num_expected = int32(settings(xml_idx).max_num_bins);
    default_num_expected = int32(512);
    num_header_fields = int32(33);
    length_field_offset = int32(260);
  elseif strcmpi(settings(xml_idx).radar_name,'ku0001')
    radar_header_type = ku0001_radar_header_type;
    min_num_expected = int32(0);
    max_num_expected = int32(settings(xml_idx).max_num_bins);
    default_num_expected = int32(512);
    num_header_fields = int32(9);
    length_field_offset = int32(68);
    
  else
    error('Not a supported radar header type %d.', radar_header_type);
  end
  
  %% Get Files for each board
  for board_idx = 1:length(param.arena_packet_strip.boards)
    board = param.arena_packet_strip.boards(board_idx);
    
    board_folder_name = fullfile(param.arena_packet_strip.board_folder_name);
    
    % Replace all "%b" in board_folder_name with the board number
    board_folder_name = regexprep(board_folder_name,'%b',sprintf('%.0f',board));
    
    % =========================================================================
    %% Get Data File list for this board
    fns = get_filenames(fullfile(base_dir,board_folder_name),'','','.dat',struct('recursive',true,'regexp','[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]'));
    fns_datenum = zeros(size(fns));
    for fn_idx = 1:length(fns)
      fname = fname_info_arena(fns{fn_idx});
      fns_datenum(fn_idx) = fname.datenum;
    end
    
    fns_mask = settings(xml_idx).xml_fname.datenum == fns_datenum;
    settings(xml_idx).fns{board_idx} = fns(fns_mask);
    
    
    %% Iterate packet_strip through file list
    old_fn_dir = [];
    for fn_idx = 1:length(settings(xml_idx).fns{board_idx})
      fn = fullfile(settings(xml_idx).fns{board_idx}{fn_idx});
      
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
      fprintf('arena_pkt_strip %d/%d %d/%d %s (%s)\n    %s\n', xml_idx, ...
        length(settings), fn_idx, length(settings(xml_idx).fns{board_idx}), fn, datestr(now), out_fn);
      
      % Check to make sure output directory exists
      if ~exist(out_fn_dir,'dir')
        mkdir(out_fn_dir);
      end
      
      % Copy XML files
      if fn_idx == 1 && board_idx == 1
        [~,out_xml_fn_name,out_xml_fn_ext] = fileparts(settings(xml_idx).xml_fn);
        out_xml_fn = ct_filename_ct_tmp(param,'','headers', ...
          fullfile(xml_folder_name, [out_xml_fn_name '.xml']));
        fprintf('Copy %s\n  %s\n', settings(xml_idx).xml_fn,out_xml_fn);
        copyfile(settings(xml_idx).xml_fn,out_xml_fn);
        if ispc
          fileattrib(out_xml_fn,'+w');
        else
          fileattrib(out_xml_fn,'+w -x');
        end
        [~,out_xml_fn_name,out_xml_fn_ext] = fileparts(settings(xml_idx).config_xml_fn);
        out_xml_fn = ct_filename_ct_tmp(param,'','headers', ...
          fullfile(xml_folder_name, [out_xml_fn_name '.xml']));
        fprintf('Copy %s\n  %s\n', settings(xml_idx).config_xml_fn,out_xml_fn);
        copyfile(settings(xml_idx).config_xml_fn,out_xml_fn);
        if ispc
          fileattrib(out_xml_fn,'+w');
        else
          fileattrib(out_xml_fn,'+w -x');
        end
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
        if strcmpi(settings(xml_idx).radar_name,'ku0001')
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
          
        elseif strcmpi(settings(xml_idx).radar_name,'KUSnow')
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
          
        elseif strcmpi(settings(xml_idx).radar_name,'TOHFSounder')
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
          
        elseif strcmpi(settings(xml_idx).radar_name,'DopplerScat')
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
[default_param,defaults] = default_radar_params_2018_Antarctica_TObas;
for xml_idx = 1:length(settings)
  
  match_idx = [];
  for default_idx = 1:length(defaults)
    match = regexpi(settings(xml_idx).psc_config_name, defaults{default_idx}.xml_regexp);
    if ~isempty(match)
      match_idx = default_idx;
      break;
    end
  end
  if isempty(match_idx)
    error('No match for psc config name %s.', settings(xml_idx).psc_config_name);
  end
  if xml_idx == 1
    oparams = default_param;
  end
  
  [~,xml_fn_name] = fileparts(settings(xml_idx).xml_fn);
  oparams(xml_idx).day_seg = sprintf('%s_%02d',xml_fn_name(1:8),xml_idx);
  oparams(xml_idx).records = defaults{match_idx}.records;
  oparams(xml_idx).qlook = defaults{match_idx}.qlook;
  oparams(xml_idx).sar = defaults{match_idx}.sar;
  oparams(xml_idx).array = defaults{match_idx}.array;
  
  oparams(xml_idx).records.file.version = 103;
  oparams(xml_idx).records.file.boards = param.arena_packet_strip.boards;
  oparams(xml_idx).records.file.prefix = datestr(settings(xml_idx).xml_fname.datenum,'YYYYmmDD_HHMMSS');
  for board_idx = 1:length(param.arena_packet_strip.boards)
    oparams(xml_idx).records.file.start_idx(board_idx) = 1;
    oparams(xml_idx).records.file.stop_idx(board_idx) = length(settings(xml_idx).fns{board_idx});
  end
  oparams(xml_idx).records.file.base_dir = param.arena_packet_strip.base_dir;
  oparams(xml_idx).records.file.board_folder_name = param.arena_packet_strip.board_folder_name;
  oparams(xml_idx).records.file.clk = 10e6;
  oparams(xml_idx).records.gps.time_offset = 0;
  oparams(xml_idx).records.gps.en = 1;
  [~,xml_fn_name] = fileparts(settings(xml_idx).xml_fn);
  oparams(xml_idx).records.xml_fn = fullfile(param.arena_packet_strip.xml_folder_name, [xml_fn_name '.xml']);
  
  if settings(xml_idx).adc{1}.adcMode == 1
    oparams(xml_idx).radar.fs = settings(xml_idx).adc{1}.sampFreq/2;
  end
  oparams(xml_idx).radar.prf = settings(xml_idx).prf;
  oparams(xml_idx).radar.adc_bits = defaults{match_idx}.radar.adc_bits;
  oparams(xml_idx).radar.Vpp_scale = defaults{match_idx}.radar.Vpp_scale;
  oparams(xml_idx).radar.lever_arm_fh = defaults{match_idx}.radar.lever_arm_fh;
  
  % Collect all waveforms
  data_map = defaults{match_idx}.records.arena.data_map;
  board_map = [];
  mode_latch_map = [];
  subchannel_map = [];
  wfs_map = [];
  adc_map = [];
  for board_idx = 1:length(data_map)
    for profile_idx = 1:size(data_map{board_idx},1)
      board_map(end+1) = param.arena_packet_strip.boards(board_idx);
      mode_latch_map(end+1) = data_map{board_idx}(profile_idx,1);
      subchannel_map(end+1) = data_map{board_idx}(profile_idx,2);
      wfs_map(end+1) = data_map{board_idx}(profile_idx,3);
      adc_map(end+1) = data_map{board_idx}(profile_idx,4);
    end
  end
  [wfs,unique_map] = unique(wfs_map);
  board = board_map(unique_map);
  mode_latch = mode_latch_map(unique_map);
  subchannel = subchannel_map(unique_map);
  adc = adc_map(unique_map);
  for wf_idx = 1:length(wfs)
    wf = wfs(wf_idx);
    
    fc = settings(xml_idx).dac{1}.wfs{mode_latch(wf_idx)+1}.centerFreq*1e6;
    BW = settings(xml_idx).dac{1}.wfs{mode_latch(wf_idx)+1}.bandwidth*1e6;
    Nt = settings(xml_idx).dac{1}.wfs{mode_latch(wf_idx)+1}.numPoints;
    fs = settings(xml_idx).dac{1}.sampFreq*1e6;
    Tpd = Nt/fs;
    t_dac = settings(xml_idx).dac{1}.wfs{mode_latch(wf_idx)+1}.initialDelay * 1e-6;
    t_arena = 3.0720e-6;
    
    oparams(xml_idx).radar.wfs(wf).f0 = fc-BW/2;
    oparams(xml_idx).radar.wfs(wf).f1 = fc+BW/2;
    oparams(xml_idx).radar.wfs(wf).tukey = settings(xml_idx).dac{1}.wfs{mode_latch(wf_idx)+1}.alpha;
    oparams(xml_idx).radar.wfs(wf).BW_window = mat2str_generic([fc-BW/2 fc+BW/2]);
    oparams(xml_idx).radar.wfs(wf).Tpd = Tpd;
    oparams(xml_idx).radar.wfs(wf).tx_weights = settings(xml_idx).dac{1}.wfs{mode_latch(wf_idx)+1}.scale;
    oparams(xml_idx).radar.wfs(wf).rx_paths = defaults{match_idx}.radar.rx_paths;
    oparams(xml_idx).radar.wfs(wf).adc_gains = defaults{match_idx}.radar.adc_gains;
    oparams(xml_idx).radar.wfs(wf).chan_equal_dB = defaults{match_idx}.radar.wfs(1).chan_equal_dB;
    oparams(xml_idx).radar.wfs(wf).chan_equal_deg = defaults{match_idx}.radar.wfs(1).chan_equal_deg;
    oparams(xml_idx).radar.wfs(wf).Tsys = defaults{match_idx}.radar.wfs(1).chan_equal_Tsys;
    oparams(xml_idx).radar.wfs(wf).presums = settings(xml_idx).psc.mode_count(mode_latch(wf_idx)+1);
    oparams(xml_idx).radar.wfs(wf).bit_shifts = ceil(max(0,log2( oparams(xml_idx).radar.wfs(wf).presums /4)));
    oparams(xml_idx).radar.wfs(wf).Tadc = sscanf(settings(xml_idx).adc{1}.rg,'%d') / oparams(xml_idx).radar.fs - t_arena - t_dac;
    
  end
%   board
%   wfs
%   mode_latch
%   subchannel
  
%   defaults{match_idx}.records.arena.data_map(:,3
%   wfs = 
  
end

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

fprintf('%s done %s\n', mfilename, datestr(now));

success = true;

return;
