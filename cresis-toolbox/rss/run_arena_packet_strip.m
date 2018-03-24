% script run_arena_packet_strip.m
%
% Strips data out of Arena packets and creates:
% 1. Files of just radar data
% 2. Files with just header data
%   a. binary flat file .hdr
%   b. Matlab .mat format (if supported)
% NOTE: This script only works on files with fixed header lengths. More
% specifically, every header must have the payload length field in the same spot.
%
% Author: John Paden
%
% See also: basic_load_arena.m

% =========================================================================
%% Get filename list (USER SECTION)
base_dir = 'Z:/HF_Sounder/2016_Greenland_TO/';
adc_folder_name = '20161108A';
reuse_tmp_files = false;
mat_or_bin_hdr_output = '.mat';

% Set param.radar_name and param.season_name
param = default_radar_params_2016_Greenland_TOdtu;

% =========================================================================
%% Read each system XML file into a system structure
system_xml_fns = get_filenames(fullfile(base_dir,adc_folder_name),'','','system.xml',struct('recursive',true));

settings = [];
for xml_idx = 1:length(system_xml_fns)
  xml_fn = system_xml_fns{xml_idx};
  [xml_fn_dir,xml_fn_name] = fileparts(xml_fn);
  xml_type_idx = find(xml_fn_name=='_',1,'last');
  config_xml_fn = fullfile(xml_fn_dir,[xml_fn_name(1:xml_type_idx), 'config.xml']);
  
  settings(xml_idx).xml_fn = xml_fn;
  settings(xml_idx).xml_fname = fname_info_arena(xml_fn);
  settings(xml_idx).config_xml_fn = config_xml_fn;
  settings(xml_idx).max_num_bins = 0;
  
  doc = xmlread(xml_fn);
  %   sys = doc.getElementsByTagName('system'); sys = sys.item(0);
  
  doc_cfg = xmlread(config_xml_fn);
  %   sys = doc.getElementsByTagName('system'); sys = sys.item(0);
  
  %   for idx = 1:sys.getLength
  %     it = sys.item(idx-1);
  %     if it.getNodeType==1
  %       name = it.getNodeName;
  %       type = it.getAttribute('type');
  %       fprintf('%s type %s\n', name.toCharArray, type.toCharArray);
  %     end
  %   end
  
  % Import the XPath classes
  import javax.xml.xpath.*
  
  % Create an XPath expression.
  factory = XPathFactory.newInstance;
  xpath = factory.newXPath;
  
  % Get the system type
  expression = xpath.compile('//system/@type');
  nodeList = expression.evaluate(doc,XPathConstants.NODESET);
  radar_name = nodeList.item(0).getTextContent.toCharArray;
  radar_name = radar_name(:).';
  settings(xml_idx).radar_name = radar_name;
  
  % Get all the ADCs
  expression = xpath.compile('//subSystem[starts-with(@type,"adc")]');
  % expression = xpath.compile('//subSystem[@type="daq"]');
  adcList = expression.evaluate(doc,XPathConstants.NODESET);
  
  for adc_idx = 1:adcList.getLength
    
    % 1. Get the name of the ADC
    expression = xpath.compile('//subSystem[starts-with(@type,"adc")]/name');
    % expression = xpath.compile('//subSystem[@type="daq"]');
    nodeList = expression.evaluate(adcList.item(adc_idx-1),XPathConstants.NODESET);
    name = nodeList.item(0).getTextContent.toCharArray;
    name = name(:).';
    
    % 2. Search for the subSystem in the config XML
    expression = xpath.compile(sprintf('//subSystem/subSystem[name="%s"]',name));
    nodeList = expression.evaluate(doc_cfg,XPathConstants.NODESET);
    match = nodeList.item(0);
    
    % 3. Get the config name and type for this ADC
    expression = xpath.compile('//config/@type');
    nodeList = expression.evaluate(match,XPathConstants.NODESET);
    config_type = nodeList.item(0).getTextContent.toCharArray;
    config_type = config_type(:).';
    expression = xpath.compile('//config/name');
    nodeList = expression.evaluate(match,XPathConstants.NODESET);
    config_name = nodeList.item(0).getTextContent.toCharArray;
    config_name = config_name(:).';
    
    % Get the config associated with this ADC
    expression = xpath.compile(sprintf('//configs/config[(@type="%s" and name="%s")]',config_type,config_name));
    nodeList = expression.evaluate(doc_cfg,XPathConstants.NODESET);
    adc_cfg = nodeList.item(0);

    % Find the longest possible record size
    if strcmpi(config_type,'adc-ads42lb69_0010')
      % TOHFSounder
      expression = xpath.compile('//processing/subChannel/mode/digRx_RG');
      nodeList = expression.evaluate(adc_cfg,XPathConstants.NODESET);
      for mode_idx = 1:nodeList.getLength
        modes = nodeList.item(mode_idx-1);
        range_gates = modes.getTextContent.toCharArray;
        range_gates = range_gates(:).';
        % Assumes simple range gate format "start:stop"
        [start,stop] = strtok(range_gates,':'); stop=stop(2:end);
        num_bins = str2double(stop) - str2double(start) + 1;
        if num_bins > settings(xml_idx).max_num_bins
          settings(xml_idx).max_num_bins = num_bins;
        end
      end
      
    elseif strcmpi(config_type,'adc-isla214p50_0005')
      % KUSnow
      expression = xpath.compile('//subChannels/subChannel/mode/rg');
      nodeList = expression.evaluate(adc_cfg,XPathConstants.NODESET);
      for mode_idx = 1:nodeList.getLength
        modes = nodeList.item(mode_idx-1);
        range_gates = modes.getTextContent.toCharArray;
        range_gates = range_gates(:).';
        % Assumes simple range gate format "start:stop"
        [start,stop] = strtok(range_gates,':'); stop=stop(2:end);
        num_bins = str2double(stop) - str2double(start) + 1;
        if num_bins > settings(xml_idx).max_num_bins
          settings(xml_idx).max_num_bins = num_bins;
        end
      end
      
    elseif strcmpi(config_type,'adc-ads42lb69_0010')
      % DopplerScat
      expression = xpath.compile('//processing/subChannel/mode/digRx_RG');
      nodeList = expression.evaluate(adc_cfg,XPathConstants.NODESET);
      for mode_idx = 1:nodeList.getLength
        modes = nodeList.item(mode_idx-1);
        range_gates = modes.getTextContent.toCharArray;
        range_gates = range_gates(:).';
        % Assumes simple range gate format "start:stop"
        [start,stop] = strtok(range_gates,':'); stop=stop(2:end);
        num_bins = str2double(stop) - str2double(start) + 1;
        if num_bins > settings(xml_idx).max_num_bins
          settings(xml_idx).max_num_bins = num_bins;
        end
      end
      
    else
      error('ADC type %s not supported.', config_type);
    end
  end
  
end

% =========================================================================
%% Get File list
fns = get_filenames(fullfile(base_dir,adc_folder_name),'','','.dat',struct('recursive',true));
fns_datenum = zeros(size(fns));
for fn_idx = 1:length(fns)
  fname = fname_info_arena(fns{fn_idx});
  fns_datenum(fn_idx) = fname.datenum;
end

% =========================================================================
%% Process each segment of data
for xml_idx = 1:length(settings)

  fns_mask = settings(xml_idx).xml_fname.datenum == fns_datenum;
  settings(xml_idx).fns = fns(fns_mask);
  
  % With each system XML file:
  %   Read in the corresponding config XML file and combine with the system structure
  %   Associate data files with the system structure based on the timestamps
  %   Packet strip using the system/config XML information
  %   Create segment information
  
  %% Initialize variables
  doppler_radar_header_type = 2;
  snow_radar_header_type = 5;
  hf_sounder_radar_header_type = 16;
  last_bytes_m = [];
  last_bytes = zeros(64,1,'uint8');
  last_bytes_len = int32(0);
  num_expected = int32(-1);
  pkt_counter = int32(-1);
  if strcmpi(radar_name,'KUSnow')
    radar_header_type = snow_radar_header_type;
    min_num_expected = int32(0);
    max_num_expected = int32(settings(xml_idx).max_num_bins);
    default_num_expected = int32(512);
    num_header_fields = int32(9);
    length_field_offset = int32(68);
  elseif strcmpi(radar_name,'TOHFSounder')
    radar_header_type = hf_sounder_radar_header_type;
    min_num_expected = int32(0);
    max_num_expected = int32(settings(xml_idx).max_num_bins);
    default_num_expected = int32(512);
    num_header_fields = int32(9);
    length_field_offset = int32(68);
  elseif strcmpi(radar_name,'DopplerScat')
    min_num_expected = int32(0);
    max_num_expected = int32(settings(xml_idx).max_num_bins);
    default_num_expected = int32(512);
    num_header_fields = int32(33);
    length_field_offset = int32(260);
  else
    error('Not a supported radar header type %d.', radar_header_type);
  end
  
  %% Iterate packet_strip through file list
  old_fn_dir = [];
  for fn_idx = 1:length(settings(xml_idx).fns)
    fn = settings(xml_idx).fns{fn_idx};
    fprintf('arena_pkt_strip %d/%d %d/%d %s (%s)\n', xml_idx, ...
      length(settings), fn_idx, length(settings(xml_idx).fns), fn, datestr(now));
    start_time = tic;
    
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
      fullfile(adc_folder_name, fn_name));
    [out_fn_dir,out_fn_name] = fileparts(out_fn);
    out_fn = fullfile(out_fn_dir,[out_fn_name,'.bin']);
    if strcmpi(mat_or_bin_hdr_output,'.mat')
      out_hdr_fn = fullfile(out_fn_dir,[out_fn_name,'.mat']);
    else
      out_hdr_fn = fullfile(out_fn_dir,[out_fn_name,'.hdr']);
    end

    % Check to make sure output directory exists
    if ~exist(out_fn_dir,'dir')
      mkdir(out_fn_dir);
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
      if strcmpi(radar_name,'KUSnow')
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
        
      elseif strcmpi(radar_name,'TOHFSounder')
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
        
      elseif strcmpi(radar_name,'DopplerScat')
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

return
