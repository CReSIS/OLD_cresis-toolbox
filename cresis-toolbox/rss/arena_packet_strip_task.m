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
  oparam = default_param;
  oparam.records = defaults{match_idx}.records;
  oparam.qlook = defaults{match_idx}.qlook;
  oparam.sar = defaults{match_idx}.sar;
  oparam.array = defaults{match_idx}.array;
  
  oparam.records.file.version = 103;
  oparam.records.file.boards = param.arena_packet_strip.boards;
  oparam.records.file.prefix = datestr(settings(xml_idx).xml_fname.datenum,'YYYYmmDD_HHMMSS');
  for board_idx = 1:length(param.arena_packet_strip.boards)
    oparam.records.file.start_idx(board_idx) = 1;
    oparam.records.file.stop_idx(board_idx) = length(settings(xml_idx).fns{board_idx});
  end
  oparam.records.file.base_dir = param.arena_packet_strip.base_dir;
  oparam.records.file.board_folder_name = param.arena_packet_strip.board_folder_name;
  oparam.records.file.clk = 10e6;
  oparam.records.gps.time_offset = 0;
  oparam.records.gps.en = 1;
  [~,xml_fn_name] = fileparts(settings(xml_idx).xml_fn);
  oparam.records.xml_fn = fullfile(param.arena_packet_strip.xml_folder_name, xml_fn_name);
  
  % Print results
end

fprintf('%s done %s\n', mfilename, datestr(now));

success = true;

return;
