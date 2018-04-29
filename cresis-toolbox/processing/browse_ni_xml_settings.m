% script browse_ni_xml_settings
%
% Reads in all the XML files from a season (probably will need to modify
% script some to get paths right) and then prints out information on each
% and which data files go with each XML files.
%
% Also has an option to write settings to the parameter spreadsheet.
%
% Author: John Paden

% dirs: cell array of strings containing directory paths that you would like
%   to process settings files for
% sub_directory: sub-directory to append to each entry in dirs
% header_load_date: Matlab datenum type, segments will use this date
% param_file: Optional struct for writing settings to a parameter file
%  .write_en: logical, if true settings will be written to parameter xls file
%  .path: optional string to override the default parameter file location
%
% Example:
%   run_browse_ni_xml_settings.m

[output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name);

%% Prepare inputs
xml_version = defaults{1}.xml_version;
cresis_xml_mapping;
if ~isempty(gps_fn)
  try
    gps = load(gps_fn);
  catch ME
    ME.getReport
    gps_fn = '';
  end
end

%% Prepare GPS plots if enabled
if ~isempty(gps_fn)
  try; close([1 2]); end;
  h_roll_fig = figure(1); clf;
  h_roll_axes = axes('Parent',h_roll_fig);
  hold(h_roll_axes,'on')
  xlabel('UTC time (seconds of day)','Parent',h_roll_axes);
  ylabel('Roll (deg)','Parent',h_roll_axes);
  grid(h_roll_axes,'on');
  h_roll_plots = [];
  h_roll_texts = [];
    
  h_geotiff = geotiff(geotiff_fn,2);
end

%% Process XML files for each directory in dirs
% =========================================================================
max_wfs = 0;
segment_offset = 0;
segment_increment = 0;
for adc_folder_name_idx = 1:length(adc_folder_names);
  segment_offset = segment_offset + segment_increment;
  
  adc_folder_name = adc_folder_names{adc_folder_name_idx};
  if exist('settings_folder_name','var')
    settings_fn_dir = fullfile(base_dir, settings_folder_name);
  else
    % Default is for the settings to be in the parent directory
    settings_fn_dir = fullfile(base_dir, adc_folder_name);
    settings_fn_dir = fileparts(settings_fn_dir);
  end
  fprintf('Processing "%s" (%d of %d)\n', settings_fn_dir, adc_folder_name_idx, length(adc_folder_names));
  
  % Read XML files in this directory
  [settings,settings_enc] = read_ni_xml_directory(settings_fn_dir,xml_file_prefix,false);
  if adc_folder_name_idx < length(adc_folder_names)
    if isempty(settings)
      fprintf('No segments (settings files) were found in %s.\n', settings_fn_dir);
      segment_increment = [];
      while length(segment_increment) ~= 1
        segment_increment = input('Input how many segments there were? ');
      end
    else
      segment_increment = segment_increment + length(settings);
    end
    proc_settings = [];
    while length(proc_settings) ~= 1
      proc_settings = input('Processing these settings now? [false]: ');
      if isempty(proc_settings)
        proc_settings = false;
      end
    end
  else
    proc_settings = true;
  end
  if ~proc_settings
    continue
  end
  
  % Get raw data files associated with this directory
  adc = defaults{1}.records.file.adcs(1);
  board = adc_to_board(param.radar_name,adc);
  adc_folder_name = regexprep(adc_folder_name,'%02d',sprintf('%02.0f',adc));
  adc_folder_name = regexprep(adc_folder_name,'%d',sprintf('%.0f',adc));
  adc_folder_name = regexprep(adc_folder_name,'%b',sprintf('%.0f',board));
  if ~isfield(defaults{1},'data_file_regexp')
    get_fns_param = [];
  else
    get_fns_param.regexp = defaults{1}.data_file_regexp;
  end
  if ~isfield(defaults{1},'data_file_midfix')
    defaults{1}.data_file_midfix = '';
  end
  data_fns = get_filenames(fullfile(base_dir,adc_folder_name),defaults{1}.data_file_prefix,defaults{1}.data_file_midfix,'.bin',get_fns_param);
  fn_datenums = [];
  
  % Prepare parameter spreadsheet for writing if enabled
  if param_file.write_en
    new_date_flag = 0;
    param_fn = ct_filename_param(sprintf('%s_param_%s.xls', ct_output_dir(param.radar_name), param.season_name));
    
    if ~exist(param_fn,'file')
      error('Parameter file %s does not exist', param_fn);
    end
    params = read_param_xls(param_fn);
    if str2double(params(1).param_file_version) < 4
      error('Cannot write to param file version less than 4');
    end
    
  end
  
  % Get the date information out of the filename
  for data_fn_idx = 1:length(data_fns)
    fname = fname_info_mcords2(data_fns{data_fn_idx});
    fn_datenums(end+1) = fname.datenum;
  end
  
  %% Print out settings from each XML file (and plot if enabled)
  for set_idx = 1:length(settings)
    settings(set_idx).enabled = true;

    % Print out settings
    [~,settings_fn_name] = fileparts(settings(set_idx).fn);
    fprintf('===================== Setting %d =================\n', set_idx);
    fprintf('%s: %d waveforms\n', settings_fn_name, length(settings(set_idx).(config_var).Waveforms));
    if isfield(settings(set_idx),'XML_File_Path')
      fprintf('  %s\n', settings(set_idx).XML_File_Path{1}.values{1});
    end
    if max_wfs < length(settings(set_idx).(config_var).Waveforms)
      max_wfs = length(settings(set_idx).(config_var).Waveforms);
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
      settings(set_idx).match_idxs = find(fn_datenums >= settings(set_idx).datenum & fn_datenums < settings(set_idx+1).datenum);
    else
      settings(set_idx).match_idxs = find(fn_datenums >= settings(set_idx).datenum);
    end
    
    % Load create_segment results and remove bad files
    settings(set_idx).day_wrap_offset = 0;
    try
      tmp_hdr_fn = ct_filename_ct_tmp(param,'','headers', ...
        fullfile(adc_folder_name, 'create_segment_raw_file_list_v2.mat'));
      create_segment = load(tmp_hdr_fn,'day_string','base_dir','param','file_prefix','file_midfix','file_regexp','segments');

      num_files = [];
      for seg_idx = 1:length(create_segment.segments)
        num_files(seg_idx) = min(create_segment.segments(seg_idx).stop_idx,settings(set_idx).match_idxs(end)) ... 
          - max(create_segment.segments(seg_idx).start_idx,settings(set_idx).match_idxs(1)) + 1;
      end
      [num_files,seg_idx] = max(num_files);
      if num_files < 1
        warning('Did not find matching segment in create_segments.');
        settings(set_idx).enabled = false;
      else
        settings(set_idx).match_idxs = max(create_segment.segments(seg_idx).start_idx,settings(set_idx).match_idxs(1)) ...
          : min(create_segment.segments(seg_idx).stop_idx,settings(set_idx).match_idxs(end));
        settings(set_idx).day_wrap_offset = create_segment.segments(seg_idx).day_wrap_offset;
      end
    end
    
    if isempty(settings(set_idx).match_idxs)
      fprintf('  No data files\n');
      settings(set_idx).match_idxs = -1;
    else
      fprintf('%5.0f files %s\n', length(settings(set_idx).match_idxs), data_fns{settings(set_idx).match_idxs(1)});
      fprintf('%5s       %s\n', '', data_fns{settings(set_idx).match_idxs(end)});
    end
    
    % Plot GPS information if enabled
    if ~isempty(gps_fn)
      % Load in first header from each file, get UTC time SOD, plot
      % position on a map
      hdr_gps_time = [];
      for match_idx = settings(set_idx).match_idxs(1:end-1)
        try
          hdr = defaults{1}.header_load_func(data_fns{match_idx},defaults{1}.header_load_params);
        end
        hdr_gps_time(end+1) = datenum_to_epoch(header_load_date + hdr.utc_time_sod/86400);
      end
      hdr_gps_time = utc_to_gps(hdr_gps_time);
      hdr_roll = interp1(gps.gps_time,gps.roll,hdr_gps_time);
      
      segment = [];
      segment.name = sprintf('%d:%d',adc_folder_name_idx,set_idx);
      segment.lat = interp1(gps.gps_time,gps.lat,hdr_gps_time);
      segment.lat = segment.lat(:);
      segment.lon = interp1(gps.gps_time,gps.lon,hdr_gps_time);
      segment.lon = segment.lon(:);
      segment.value_name = {'Filename','Roll'};
      segment.value = data_fns(settings(set_idx).match_idxs(1:end-1));
      segment.value = segment.value(:);
      segment.value = cat(2,segment.value,mat2cell(hdr_roll(:)*180/pi,ones(size(hdr_roll(:))),1));
      color = h_geotiff.insert_segment(segment);
      
      hdr_utc_time_sod = epoch_to_sod(gps_to_utc(hdr_gps_time),datenum_to_epoch(header_load_date));
      if isempty(hdr_utc_time_sod)
        h_roll_plots(adc_folder_name_idx,set_idx) = plot(NaN,NaN,'.-','Parent',h_roll_axes,'Color',color);
        h_roll_texts(adc_folder_name_idx,set_idx) = text(NaN,NaN,sprintf('%d:%d',adc_folder_name_idx,set_idx),'Color',color,'Parent',h_roll_axes);
      else
        h_roll_plots(adc_folder_name_idx,set_idx) = plot(hdr_utc_time_sod,hdr_roll*180/pi,'.-','Parent',h_roll_axes,'Color',color);
        h_roll_texts(adc_folder_name_idx,set_idx) = text(hdr_utc_time_sod(1),hdr_roll(1)*180/pi,sprintf('%d:%d',adc_folder_name_idx,set_idx),'Color',color,'Parent',h_roll_axes);
      end
    end
    
    % Associate default parameters with each settings
    settings(set_idx).default = default_radar_params_settings_match(defaults,settings(set_idx));
    default = settings(set_idx).default;

    % Determine if using these settings
    if settings(set_idx).match_idxs(end) - settings(set_idx).match_idxs(1) ...
        < MIN_FILES_IN_SEGMENT-1+SKIP_FILES_START+SKIP_FILES_END || ~isempty(regexpi(default.name,'other'))
      settings(set_idx).enabled = false;
    end
    if manual_enable
      user_input = input(sprintf('Keep segment? [%s]: ', mat2str(settings(set_idx).enabled)));
      if ~isempty(user_input)
        settings(set_idx).enabled = logical(user_input);
      end
    end

  end
  
  %% Vectors Worksheet
  fprintf('\nVectors WorkSheet\n');
  fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'Date', ...
    '1', 'file.start_idx', 'file.stop_idx', 'file.base_dir', 'file.adc_folder_name', 'file.file_prefix', 'file.file_midfix','file.regexp','gps.time_offset');
  fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
    'YYYYMMDD','Segment','r','r','t','t','t','t','r','r');
  for set_idx = 1:length(settings)
    if settings(set_idx).enabled
      default = settings(set_idx).default;
      
      fprintf('%s\t%02d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%g\n', datestr(header_load_date,'YYYYmmDD'), ...
        segment_offset+set_idx, settings(set_idx).match_idxs(1)+SKIP_FILES_START, ...
        settings(set_idx).match_idxs(end)-SKIP_FILES_END, base_dir,  ...
        adc_folder_names{adc_folder_name_idx}, defaults{1}.data_file_prefix, '', ...
        '', default.vectors.gps.time_offset+settings(set_idx).day_wrap_offset);
      
      % Write to parameter file if enabled
      if param_file.write_en
        day_seg = sprintf('%s_%02.0f', datestr(header_load_date,'YYYYmmDD'), set_idx);
        segment_found = false;
        if isempty(params(1).day_seg)
          param_idx = 1;
        else
          for param_idx = 1:length(params)
            if strcmpi(params(param_idx).day_seg,day_seg)
              segment_found = true;
              break;
            end
          end
          
          if ~segment_found
            param_idx = length(params)+1;
          end
        end
        settings(set_idx).param_idx = param_idx;
        if ~segment_found
          params(param_idx).day_seg = day_seg;
        end
        
        params(param_idx).vectors.file.radar_num = [];
        params(param_idx).vectors.file.adc = 1;
        start_idx = settings(set_idx).match_idxs(1)+SKIP_FILES_START;
        stop_idx = settings(set_idx).match_idxs(end)-SKIP_FILES_END;
        if segment_found
          if start_idx ~= params(param_idx).vectors.file.start_idx
            warning('params(%d).vectors.file.start_idx %d does not match %d. Verify that you would like to overwrite by typing dbcont.', ...
              param_idx,params(param_idx).vectors.file.start_idx,start_idx);
            keyboard
          end
          if stop_idx ~= params(param_idx).vectors.file.stop_idx
            warning('params(%d).vectors.file.stop_idx %d does not match %d. Verify that you would like to overwrite by typing dbcont.', ...
              param_idx,params(param_idx).vectors.file.stop_idx,stop_idx);
            keyboard
          end
        end
        params(param_idx).vectors.file.start_idx = start_idx;
        params(param_idx).vectors.file.stop_idx = stop_idx;
        params(param_idx).vectors.file.base_dir = base_dir_in_param;
        params(param_idx).vectors.file.adc_folder_name = adc_folder_names{adc_folder_name_idx};
        params(param_idx).vectors.file.file_prefix = data_file_prefix;
        params(param_idx).vectors.out_fn = '';
        params(param_idx).vectors.gps.fn = '';
        params(param_idx).vectors.gps.time_offset = 1;
        params(param_idx).vectors.gps.verification = '';
        params(param_idx).vectors.gps.utc_time_halved = 0;
        
        params(param_idx).cmd = merge_structs(params(param_idx).cmd,default.cmd);
        if isempty(params(param_idx).cmd.notes)
          [~,xml_fn_name] = fileparts(settings(set_idx).XML_File_Path{1}.values{1});
          params(param_idx).cmd.notes = xml_fn_name;
        end
        params(param_idx).records = merge_structs(params(param_idx).records,default.records);
      end

    end
  end
  
  %% Command Worksheet
  fprintf('\nCommand (cmd) WorkSheet\n');
  fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
    'Date','1','frms','create_records','create_frames','get_heights','csarp','combine_wf_chan','generic','mission_names','notes');
  fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
    'YYYYMMDD','Segment','r','b','b','b','b','b','r','t','t');
  for set_idx = 1:length(settings)
    if settings(set_idx).enabled
      default = settings(set_idx).default;

      % Print out data row
      if isfield(settings(set_idx),'XML_File_Path')
        settings_fn = settings(set_idx).XML_File_Path{1}.values{1};
      else
        settings_fn = '';
      end
      fprintf('%s\t%02d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
        datestr(header_load_date,'YYYYmmDD'),segment_offset+set_idx,'','','','','','','',mission_name,settings_fn);

    end
  end
  
  %% Records Worksheet
  fprintf('\nRecords WorkSheet\n');
  fprintf('%s\t%s\t%s\n', ...
    'Date','1','file_version');
  fprintf('%s\t%s\t%s\n', ...
    'YYYYMMDD','Segment','r');
  for set_idx = 1:length(settings)
    if settings(set_idx).enabled
      default = settings(set_idx).default;
            
      % Determine file version
      if strcmpi(radar_name,'mcords3')
        file_version = 403;
      elseif any(strcmpi(radar_name,{'mcords4'}))
        file_version = 404;
      elseif any(strcmpi(radar_name,{'mcords5'}))
        if settings(set_idx).DDC_Ctrl.DDC_sel.Val == 0
          file_version = 408;
        elseif any(settings(set_idx).DDC_Ctrl.DDC_sel.Val == [1 2])
          file_version = 407;
        else
          error('Invalid DDC setting.');
        end
      end
      
      % Print out data row
      fprintf('%s\t%02d\t%d\n',datestr(header_load_date,'YYYYmmDD'),segment_offset+set_idx, file_version);

    end
  end
    
  %% Get Heights Worksheet
  fprintf('\nGet Heights (Quick Look) WorkSheet\n');
  fprintf('%s\t%s\t%s\t%s\t%s\n', ...
    'Date','1','qlook.out_path','qlook.img_comb','imgs');
  fprintf('%s\t%s\t%s\t%s\t%s\n', ...
    'YYYYMMDD','Segment','t','r','r');
  
  for set_idx = 1:length(settings)
    default = settings(set_idx).default;
    if settings(set_idx).enabled
      % Print out data row
      fprintf('%s\t%02d\t%s\t%s\t%s\n',datestr(header_load_date,'YYYYmmDD'),segment_offset+set_idx, ...
        default.get_heights.qlook.out_path, mat2str(default.get_heights.qlook.img_comb), mat2str_generic(default.get_heights.imgs));
      
      % Write to parameter file if enabled
      if param_file.write_en
        param_idx = settings(set_idx).param_idx;
        params(param_idx).get_heights = merge_structs(params(param_idx).get_heights,default.get_heights);
      end
    end
  end
  
  %% CSARP Worksheet
  fprintf('\nCSARP (SAR Processor) WorkSheet\n');
  fprintf('%s\t%s\t%s\t%s\n', ...
    'Date','1','out_path','imgs');
  fprintf('%s\t%s\t%s\t%s\n', ...
    'YYYYMMDD','Segment','t','r');
  
  for set_idx = 1:length(settings)
    default = settings(set_idx).default;
    if settings(set_idx).enabled
      % Print out data row
      fprintf('%s\t%02d\t%s\t%s\n',datestr(header_load_date,'YYYYmmDD'),segment_offset+set_idx, ...
        default.csarp.out_path, mat2str_generic(default.csarp.imgs));
      
      % Write to parameter file if enabled
      if param_file.write_en
        param_idx = settings(set_idx).param_idx;
        
        params(param_idx).csarp = merge_structs(params(param_idx).csarp,default.csarp);
      end
    end
  end
    
  %% Combine Wf Chan Worksheet
  fprintf('\nCombine (Combine Waveforms and Images) WorkSheet\n');
  fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
    'Date','1','in_path','array_path','out_path','img_comb','imgs');
  fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
    'YYYYMMDD','Segment','t','t','t','r','r');
  
  for set_idx = 1:length(settings)
    default = settings(set_idx).default;
    if settings(set_idx).enabled
      % Print out data row
      fprintf('%s\t%02d\t%s\t%s\t%s\t%s\t%s\n',datestr(header_load_date,'YYYYmmDD'),segment_offset+set_idx, ...
        default.combine.in_path, default.combine.array_path, default.combine.out_path, ...
        mat2str(default.combine.img_comb), mat2str_generic(default.combine.imgs));
      
      % Write to parameter file if enabled
      if param_file.write_en
        param_idx = settings(set_idx).param_idx;
        
        params(param_idx).combine = merge_structs(params(param_idx).combine,default.combine);
      end
    end
  end
  
  %% Radar Worksheet
  fprintf('\nRadar WorkSheet\n');
  fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s', ...
    'Date','1','fs','prf','adc_bits','Vpp_scale','size');
  for wf = 1:max_wfs
    for wf_hdr_idx = 1:length(default.radar_worksheet_headers)
      fprintf('\t%s', default.radar_worksheet_headers{wf_hdr_idx});
    end
  end
  fprintf('\n');
  fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s', ...
    'YYYYMMDD','Segment','r','r','r','r','ar:wfs');
  for wf = 1:max_wfs
    for wf_hdr_idx = 1:length(default.radar_worksheet_headers)
      fprintf('\ta%s:wfs(%d)', default.radar_worksheet_headers_type{wf_hdr_idx}, wf);
    end
  end
  fprintf('\n');
  
  for set_idx = 1:length(settings)
    default = settings(set_idx).default;
      
    if settings(set_idx).enabled
      
      if param_file.write_en
        param_idx = settings(set_idx).param_idx;
        params(param_idx).radar.fs = default.radar.fs;
        params(param_idx).radar.prf = settings(set_idx).(config_var).(prf_var);
        params(param_idx).radar.adc_bits = default.radar.adc_bits;
        params(param_idx).radar.Vpp_scale = default.radar.adc_full_scale;
      end
      
      row_str = sprintf('%s\t%02d\t%.16g\t%g\t%d\t%g\t%d',datestr(header_load_date,'YYYYmmDD'),segment_offset+set_idx, ...
        default.radar.fs, settings(set_idx).(config_var).(prf_var), default.radar.adc_bits, default.radar.adc_full_scale, length(settings(set_idx).(config_var).Waveforms));
      
      for wf = 1:length(settings(set_idx).(config_var).Waveforms)
        
        if any(strcmpi('Tpd',default.radar_worksheet_headers))
          if param_file.write_en
            params(param_idx).radar.wfs(wf).Tpd = double(settings(set_idx).(config_var).Waveforms(wf).Len_Mult)*settings(set_idx).(config_var).Base_Len;
          end
          row_str = cat(2,row_str, sprintf('\t%g',double(settings(set_idx).(config_var).Waveforms(wf).Len_Mult)*settings(set_idx).(config_var).Base_Len));
        end
        if any(strcmpi('Tadc',default.radar_worksheet_headers))
          if param_file.write_en
            params(param_idx).radar.wfs(wf).Tadc = default.radar.Tadc;
          end
          row_str = cat(2,row_str, sprintf('\t%g',default.radar.Tadc));
        end
        if any(strcmpi('Tadc_adjust',default.radar_worksheet_headers))
          if param_file.write_en
            params(param_idx).radar.wfs(wf).Tadc_adjust = default.radar.Tadc_adjust;
          end
          str = char(format(java.text.DecimalFormat('0.0000000000'),default.radar.Tadc_adjust));
%          row_str = cat(2,row_str, sprintf('\t%g',default.radar.Tadc_adjust));
          row_str = cat(2,row_str, sprintf('\t%s',str));
        end
        if any(strcmpi('f0',default.radar_worksheet_headers))
          if param_file.write_en
            params(param_idx).radar.wfs(wf).f0 = settings(set_idx).(config_var).Waveforms(wf).Start_Freq(1);
          end
          row_str = cat(2,row_str, sprintf('\t%.16g',settings(set_idx).(config_var).Waveforms(wf).Start_Freq(1)));
        end
        if any(strcmpi('f1',default.radar_worksheet_headers))
          if param_file.write_en
            params(param_idx).radar.wfs(wf).f1 = settings(set_idx).(config_var).Waveforms(wf).Stop_Freq(1);
          end
          row_str = cat(2,row_str, sprintf('\t%.16g',settings(set_idx).(config_var).Waveforms(wf).Stop_Freq(1)));
        end
        if any(strcmpi('ft_dec',default.radar_worksheet_headers))
          if param_file.write_en
            params(param_idx).radar.wfs(wf).ft_dec = default.radar.ft_dec;
          end
          row_str = cat(2,row_str, ...
              sprintf('\t[%d',default.radar.ft_dec(1)));
          row_str = cat(2,row_str, ...
              sprintf(' %d',default.radar.ft_dec(2)));
         row_str = cat(2,row_str, ...
              sprintf(']'));         
        end
        if any(strcmpi('ref_fn',default.radar_worksheet_headers))
          if param_file.write_en
            params(param_idx).radar.wfs(wf).ref_fn = default.radar.ref_fn;
          end
          row_str = cat(2,row_str, sprintf('\t''%s''',default.radar.ref_fn));
        end
        if any(strcmpi('tukey',default.radar_worksheet_headers))
          if param_file.write_en
            params(param_idx).radar.wfs(wf).tukey = settings(set_idx).(config_var).RAM_Taper;
          end
          row_str = cat(2,row_str, sprintf('\t%g',settings(set_idx).(config_var).RAM_Taper));
        end
        
        % Transmit weights
        if any(strcmpi(radar_name,{'mcords3','mcords5'}))
          tx_mask_inv = fliplr(~(dec2bin(double(settings(set_idx).(config_var).Waveforms(wf).TX_Mask),8) - '0'));
          tx_weights = double(settings(set_idx).(config_var).(ram_var)) .* tx_mask_inv / default.max_DDS_RAM*default.tx_voltage;
        else
          tx_mask_inv = ~(dec2bin(double(settings(set_idx).(config_var).Waveforms(wf).TX_Mask),8) - '0');
          tx_weights = double(settings(set_idx).(config_var).(ram_var)) .* tx_mask_inv / default.max_DDS_RAM*default.tx_voltage;
        end
        tx_weights = tx_weights(logical(default.tx_DDS_mask));
        if param_file.write_en
          params(param_idx).radar.wfs(wf).tx_weights = tx_weights;
        end
        row_str = cat(2,row_str, ...
          sprintf('\t[%g', tx_weights(1)));
        row_str = cat(2,row_str, ...
          sprintf(' %g', tx_weights(2:end)));
        row_str = cat(2,row_str, ...
          sprintf(']'));
        
        % Rx paths
        if param_file.write_en
          params(param_idx).radar.wfs(wf).rx_paths = default.radar.rx_paths;
        end
        row_str = cat(2,row_str, ...
          sprintf('\t[%d', default.radar.rx_paths(1)));
        row_str = cat(2,row_str, ...
          sprintf(' %d', default.radar.rx_paths(2:end)));
        row_str = cat(2,row_str, ...
          sprintf(']'));
        
        % ADC Gains
        atten = double(settings(set_idx).(config_var).Waveforms(wf).Attenuator_1(1)) ...
          + double(settings(set_idx).(config_var).Waveforms(wf).Attenuator_2(1));
        if param_file.write_en
          params(param_idx).radar.wfs(wf).adc_gains = 10.^((default.radar.rx_gain - atten(1)*ones(1,length(default.radar.rx_paths)))/20);
        end
        row_str = cat(2,row_str, ...
          sprintf('\t10.^((%g - %g*ones(1,%d))/20)', default.radar.rx_gain, atten(1), length(default.radar.rx_paths)));
        
        default_wf = min(wf,length(default.radar.wfs));
        % Chan Equal DB
        if param_file.write_en
          params(param_idx).radar.wfs(wf).chan_equal_dB = default.radar.wfs(default_wf).chan_equal_dB;
        end
        row_str = cat(2,row_str, sprintf('\t%s',mat2str(default.radar.wfs(default_wf).chan_equal_dB)));
        
        % Chan Equal deg
        if param_file.write_en
          params(param_idx).radar.wfs(wf).chan_equal_deg = default.radar.wfs(default_wf).chan_equal_deg;
        end
        row_str = cat(2,row_str, sprintf('\t%s',mat2str(default.radar.wfs(default_wf).chan_equal_deg)));
        
        % Tsys
        if param_file.write_en
          params(param_idx).radar.wfs(wf).Tsys = default.radar.wfs(default_wf).chan_equal_Tsys;
        end
        row_str = cat(2,row_str, sprintf('\t%s',mat2str(default.radar.wfs(default_wf).chan_equal_Tsys)));
        
        % DC Adjust
        if param_file.write_en
          if wf <= length(default.radar.DC_adjust)
            params(param_idx).radar.wfs(wf).DC_adjust = '';
          else
            params(param_idx).radar.wfs(wf).DC_adjust = default.radar.DC_adjust{wf};
          end
        end
        if wf <= length(default.radar.DC_adjust)
          row_str = cat(2,row_str, sprintf('\t''%s''',default.radar.DC_adjust{wf}));
        else
          row_str = cat(2,row_str, sprintf('\t%s',''));
        end
        
        % DDC mode and frequency
        if any(strcmpi('DDC_mode',default.radar_worksheet_headers))
          if param_file.write_en
            params(param_idx).radar.wfs(wf).DDC_mode = settings(set_idx).DDC_Ctrl.DDC_sel.Val;
          end
          row_str = cat(2,row_str, sprintf('\t%d',settings(set_idx).DDC_Ctrl.DDC_sel.Val));
        end
        if any(strcmpi('DDC_freq',default.radar_worksheet_headers))
          if param_file.write_en
            params(param_idx).radar.wfs(wf).DDC_freq = settings(set_idx).DDC_Ctrl.(NCO_freq)*1e6;
          end
          row_str = cat(2,row_str, sprintf('\t%g',settings(set_idx).DDC_Ctrl.(NCO_freq)*1e6));
        end
        
      end
      fprintf('%s\n',row_str)
    end
  end

  %% Write to parameter spreadsheet
  if param_file.write_en
    insert_param_xls(param_fn,params,[],'wfs');
    try
    catch ME
      ME.getReport
      keyboard
    end
  end
  
end

return;
