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

%% Prepare inputs
cresis_xml_mapping;
if ~isempty(gps_fn)
  gps = load(gps_fn);
end

%% Prepare GPS plots if enabled
if ~isempty(gps_fn)
  h_roll_fig = figure(1); clf;
  h_roll_axes = axes('Parent',h_roll_fig);
  h_fig = figure(2); clf;
  if ~isempty(geotiff_fn)
    [proj] = plot_geotiff(geotiff_fn,[],[],h_fig);
    h_axes = get(h_fig,'Children');
  else
    h_axes = axes('Parent',h_fig);
  end
  hold(h_axes,'on')
  hold(h_roll_axes,'on')
  h_plots = [];
  h_texts = [];
  h_roll_plots = [];
  h_roll_texts = [];
end

%% Process XML files for each directory in dirs
max_wfs = 0;
for adc_folder_name_idx = 1:length(adc_folder_names);
  
  adc_folder_name = adc_folder_names{adc_folder_name_idx};
  settings_fn_dir = fullfile(base_dir, adc_folder_name);
  fprintf('Processing "%s" (%d of %d)\n', settings_fn_dir, adc_folder_name_idx, length(adc_folder_names));
  
  %% Read XML files in this directory
  [settings,settings_enc] = read_ni_xml_directory(settings_fn_dir,xml_file_prefix,false);
  
  %% Get raw data files associated with this directory
  data_fns = get_filenames(fullfile(base_dir,adc_folder_name,board_sub_directory),data_file_prefix,'','.bin');
  fn_datenums = [];
  
  %% Prepare parameter spreadsheet for writing if enabled
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
  
  %% Get the date information out of the filename
  for data_fn_idx = 1:length(data_fns)
    fname = fname_info_mcords2(data_fns{data_fn_idx});
    fn_datenums(end+1) = fname.datenum;
  end
  
  %% Print out settings from each XML file (and plot if enabled)
  for set_idx = 1:length(settings)
    %% Print out settings
    [~,settings_fn_name] = fileparts(settings(set_idx).fn);
    fprintf('===================== Setting %d =================\n', set_idx);
    fprintf('%s: %d waveforms\n', settings_fn_name, length(settings(set_idx).(config_var).Waveforms));
    if max_wfs < length(settings(set_idx).(config_var).Waveforms)
      max_wfs = length(settings(set_idx).(config_var).Waveforms);
    end
    fprintf('  PRF:'); fprintf(' %g', settings(set_idx).(config_var).(prf_var)); fprintf('\n');
    fprintf('  Amp:'); fprintf(' %g', settings(set_idx).(config_var).(ram_amp_var)); fprintf('\n');
    fprintf('  Tukey:'); fprintf(' %g', settings(set_idx).(config_var).RAM_Taper); fprintf('\n');
    fprintf('  Stop:'); fprintf(' %g', settings(set_idx).(config_var).Waveforms(1).Stop_Freq); fprintf('\n');
    fprintf('  Tx Mask:'); fprintf(' %g', settings(set_idx).(config_var).Waveforms(1).TX_Mask); fprintf('\n');
    for wf = 1:length(settings(set_idx).(config_var).Waveforms)
      fprintf('  WF %d Atten:', wf); fprintf(' %g', settings(set_idx).(config_var).Waveforms(wf).Attenuator_2); fprintf('\n');
      fprintf('  WF %d Len:', wf); fprintf(' %.1f us', 1e6*settings(set_idx).(config_var).Base_Len*settings(set_idx).(config_var).Waveforms(wf).Len_Mult); fprintf('\n');
    end
    if set_idx < length(settings)
      settings(set_idx).match_idxs = find(fn_datenums >= settings(set_idx).datenum & fn_datenums < settings(set_idx+1).datenum);
    else
      settings(set_idx).match_idxs = find(fn_datenums >= settings(set_idx).datenum);
    end
    if isempty(settings(set_idx).match_idxs)
      fprintf('    No data files\n');
      settings(set_idx).match_idxs = -1;
    else
      fprintf('    %s\n', data_fns{settings(set_idx).match_idxs(1)});
      fprintf('    %s\n', data_fns{settings(set_idx).match_idxs(end)});
    end
    
    %% Plot GPS information if enabled
    if ~isempty(gps_fn)
      % Load in first header from each file, get UTC time SOD, plot
      % position on a map
      for match_idx = settings(set_idx).match_idxs(1:end-1)
        try
          hdr = header_load_func(data_fns{match_idx},header_load_params);
        end
        hdr_gps_time = datenum_to_epoch(header_load_date + hdr.utc_time_sod/86400);
        hdr_gps_time = hdr_gps_time + utc_leap_seconds(hdr_gps_time);
        lat = interp1(gps.gps_time,gps.lat,hdr_gps_time);
        lon = interp1(gps.gps_time,gps.lon,hdr_gps_time);
        if ~isempty(geotiff_fn)
          [x,y] = projfwd(proj,lat,lon);
        else
          x = lon;
          y = lat;
        end
        h_plots(dir_idx,set_idx) = plot(x/1e3,y/1e3,'.','Parent',h_axes);
        hdr_roll = interp1(gps.gps_time,gps.roll,hdr_gps_time);
        h_roll_plots(dir_idx,set_idx) = plot(hdr_gps_time,hdr_roll,'Parent',h_roll_axes);
      end
      h_texts(dir_idx,set_idx) = text(x(1)/1e3,y(1)/1e3,sprintf('%d:%d',dir_idx,set_idx),'Parent',h_axes);
      h_roll_texts(dir_idx,set_idx) = text(hdr_gps_time(1),hdr_roll(1),sprintf('%d:%d',dir_idx,set_idx),'Parent',h_roll_axes);
    end
  end
  
  %% Vectors Worksheet
  fprintf('Vectors and Records WorkSheet\n');
  fprintf('%s\t%s\t%s\t%s\t%s\n', ...
    'Date','1','file.start_idx','file.stop_idx','file_version');
  fprintf('%s\t%s\t%s\t%s\t%s\n', ...
    'YYYYMMDD','Segment','r','r','r');
  for set_idx = 1:length(settings)
    if settings(set_idx).match_idxs(end) - settings(set_idx).match_idxs(1) >= MIN_FILES_IN_SEGMENT-1
      
      %% Determine file version
      if strcmpi(param.radar_name,'mcords3')
        file_version = 403;
      elseif any(strcmpi(param.radar_name,{'mcords4'}))
        file_version = 404;
      elseif any(strcmpi(param.radar_name,{'mcords5'}))
        if settings(set_idx).DDC_Ctrl.DDC_sel.Val == 0
          file_version = 408;
        elseif any(settings(set_idx).DDC_Ctrl.DDC_sel.Val == [1 2])
          file_version = 407;
        else
          error('Invalid DDC setting.');
        end
      end
      
      %% Print out data row
      fprintf('%s\t%02d\t%d\t%d\t%d\n',datestr(header_load_date,'YYYYmmDD'),set_idx, ...
        settings(set_idx).match_idxs(1)+SKIP_FILES_START, settings(set_idx).match_idxs(end)-SKIP_FILES_END, ...
        file_version);
      
      %% Write to parameter file if enabled
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
        params(param_idx).vectors.file.adc_folder_name = sprintf('/%s/chan%%d/',adc_folder_name);
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
  
  %% Radar Worksheet
  fprintf('Radar WorkSheet\n');
  fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s', ...
    'Date','1','fs','prf','adc_bits','Vpp_scale','size');
  for wf = 1:max_wfs
    for wf_hdr_idx = 1:length(radar_worksheet_headers)
      fprintf('\t%s', radar_worksheet_headers{wf_hdr_idx});
    end
  end
  fprintf('\n');
  fprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s', ...
    'YYYYMMDD','Segment','r','r','r','r','ar:wfs');
  for wf = 1:max_wfs
    for wf_hdr_idx = 1:length(radar_worksheet_headers)
      fprintf('\ta%s:wfs(%d)', radar_worksheet_headers_type{wf_hdr_idx}, wf);
    end
  end
  fprintf('\n');
  
  for set_idx = 1:length(settings)
    if settings(set_idx).match_idxs(end) - settings(set_idx).match_idxs(1) >= MIN_FILES_IN_SEGMENT-1
      
      if param_file.write_en
        param_idx = settings(set_idx).param_idx;
        params(param_idx).radar.fs = fs;
        params(param_idx).radar.prf = settings(set_idx).(config_var).(prf_var);
        params(param_idx).radar.adc_bits = adc_bits;
        params(param_idx).radar.Vpp_scale = adc_full_scale;
      end
      
      row_str = sprintf('%s\t%02d\t%.16g\t%g\t%d\t%g\t%d',datestr(header_load_date,'YYYYmmDD'),set_idx, ...
        fs, settings(set_idx).(config_var).(prf_var), adc_bits, adc_full_scale, length(settings(set_idx).(config_var).Waveforms));
      
      for wf = 1:length(settings(set_idx).(config_var).Waveforms)
        
        if any(strcmpi('Tpd',radar_worksheet_headers))
          if param_file.write_en
            params(param_idx).radar.wfs(wf).Tpd = double(settings(set_idx).(config_var).Waveforms(wf).Len_Mult)*settings(set_idx).(config_var).Base_Len;
          end
          row_str = cat(2,row_str, sprintf('\t%g',double(settings(set_idx).(config_var).Waveforms(wf).Len_Mult)*settings(set_idx).(config_var).Base_Len));
        end
        if any(strcmpi('Tadc',radar_worksheet_headers))
          if param_file.write_en
            params(param_idx).radar.wfs(wf).Tadc = Tadc;
          end
          row_str = cat(2,row_str, sprintf('\t%g',Tadc));
        end
        if any(strcmpi('Tadc_adjust',radar_worksheet_headers))
          if param_file.write_en
            params(param_idx).radar.wfs(wf).Tadc_adjust = Tadc_adjust;
          end
          row_str = cat(2,row_str, sprintf('\t%g',Tadc_adjust));
        end
        if any(strcmpi('f0',radar_worksheet_headers))
          if param_file.write_en
            params(param_idx).radar.wfs(wf).f0 = settings(set_idx).(config_var).Waveforms(wf).Start_Freq(1);
          end
          row_str = cat(2,row_str, sprintf('\t%.16g',settings(set_idx).(config_var).Waveforms(wf).Start_Freq(1)));
        end
        if any(strcmpi('f1',radar_worksheet_headers))
          if param_file.write_en
            params(param_idx).radar.wfs(wf).f1 = settings(set_idx).(config_var).Waveforms(wf).Stop_Freq(1);
          end
          row_str = cat(2,row_str, sprintf('\t%.16g',settings(set_idx).(config_var).Waveforms(wf).Stop_Freq(1)));
        end
        if any(strcmpi('ref_fn',radar_worksheet_headers))
          if param_file.write_en
            params(param_idx).radar.wfs(wf).ref_fn = '';
          end
          row_str = cat(2,row_str, sprintf('\t'));
        end
        if any(strcmpi('tukey',radar_worksheet_headers))
          if param_file.write_en
            params(param_idx).radar.wfs(wf).tukey = settings(set_idx).(config_var).RAM_Taper;
          end
          row_str = cat(2,row_str, sprintf('\t%g',settings(set_idx).(config_var).RAM_Taper));
        end
        
        % Transmit weights
        if any(strcmpi(param.radar_name,{'mcords3','mcords5'}))
          tx_mask_inv = fliplr(~(dec2bin(double(settings(set_idx).(config_var).Waveforms(wf).TX_Mask),8) - '0'));
          tx_weights = double(settings(set_idx).(config_var).(ram_var)) .* tx_mask_inv / max_DDS_RAM*tx_voltage;
        else
          tx_mask_inv = ~(dec2bin(double(settings(set_idx).(config_var).Waveforms(wf).TX_Mask),8) - '0');
          tx_weights = double(settings(set_idx).(config_var).(ram_var)) .* tx_mask_inv / max_DDS_RAM*tx_voltage;
        end
        tx_weights = tx_weights(logical(tx_DDS_mask));
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
          params(param_idx).radar.wfs(wf).rx_paths = rx_paths;
        end
        row_str = cat(2,row_str, ...
          sprintf('\t[%d', rx_paths(1)));
        row_str = cat(2,row_str, ...
          sprintf(' %d', rx_paths(2:end)));
        row_str = cat(2,row_str, ...
          sprintf(']'));
        
        % ADC Gains
        atten = double(settings(set_idx).(config_var).Waveforms(wf).Attenuator_1(1)) ...
          + double(settings(set_idx).(config_var).Waveforms(wf).Attenuator_2(1));
        if param_file.write_en
          params(param_idx).radar.wfs(wf).adc_gains = 10.^((max_adc_gain_dB - atten(1)*ones(1,length(rx_paths)))/20);
        end
        row_str = cat(2,row_str, ...
          sprintf('\t10.^((%g - %g*ones(1,%d))/20)', max_adc_gain_dB, atten(1), length(rx_paths)));
        
        % Chan Equal DB
        if param_file.write_en
          params(param_idx).radar.wfs(wf).chan_equal_dB = chan_equal_dB;
        end
        row_str = cat(2,row_str, sprintf('\t%s',chan_equal_dB));
        
        % Chan Equal deg
        if param_file.write_en
          params(param_idx).radar.wfs(wf).chan_equal_deg = chan_equal_deg;
        end
        row_str = cat(2,row_str, sprintf('\t%s',chan_equal_deg));
        
        % Tsys
        if param_file.write_en
          params(param_idx).radar.wfs(wf).Tsys = chan_equal_Tsys;
        end
        row_str = cat(2,row_str, sprintf('\t%s',chan_equal_Tsys));
        if any(strcmpi('DDC_mode',radar_worksheet_headers))
          if param_file.write_en
            params(param_idx).radar.wfs(wf).DDC_mode = settings(set_idx).DDC_Ctrl.DDC_sel.Val;
          end
          row_str = cat(2,row_str, sprintf('\t%d',settings(set_idx).DDC_Ctrl.DDC_sel.Val));
        end
        if any(strcmpi('DDC_freq',radar_worksheet_headers))
          if param_file.write_en
            params(param_idx).radar.wfs(wf).DDC_freq = settings(set_idx).DDC_Ctrl.(NCO_freq)*1e6;
          end
          row_str = cat(2,row_str, sprintf('\t%g',settings(set_idx).DDC_Ctrl.(NCO_freq)*1e6));
        end
        
      end
      fprintf('%s\n',row_str)
    end
  end
  
  %% Processing Worksheets: Get Heights, CSARP, and Combine Wf Chan
  fprintf('Processing WorkSheets (Get Heights, CSARP, and Combine Wf Chan)\n');
  for set_idx = 1:length(settings)
    if settings(set_idx).match_idxs(end) - settings(set_idx).match_idxs(1) >= MIN_FILES_IN_SEGMENT-1
      if abs(iq_mode)
        wfs = 1:2:length(settings(set_idx).(config_var).Waveforms);
      else
        wfs = 1:length(settings(set_idx).(config_var).Waveforms);
      end
      Tpd= double(cell2mat({settings(set_idx).(config_var).Waveforms(wfs).Len_Mult}))*settings(set_idx).(config_var).Base_Len;
      [~,sort_idx] = sort(Tpd);
      
      if length(settings(set_idx).(config_var).Waveforms) > 1+abs(iq_mode)
        wf_idx = 2;
        Tpd_prev = double(settings(set_idx).(config_var).Waveforms(wfs(wf_idx-1)).Len_Mult)*settings(set_idx).(config_var).Base_Len;
        Tpd = double(settings(set_idx).(config_var).Waveforms(wfs(wf_idx)).Len_Mult)*settings(set_idx).(config_var).Base_Len;
        row_str = sprintf('[%g -inf %g', Tpd, Tpd_prev);
        for wf_idx = 3:length(wfs)
          Tpd_prev = double(settings(set_idx).(config_var).Waveforms(wfs(wf_idx-1)).Len_Mult)*settings(set_idx).(config_var).Base_Len;
          Tpd = double(settings(set_idx).(config_var).Waveforms(wfs(wf_idx)).Len_Mult)*settings(set_idx).(config_var).Base_Len;
          row_str = cat(2,row_str,sprintf(' %g -inf %g', Tpd, Tpd_prev));
        end
        row_str = cat(2,row_str,sprintf(']\t'));
      else
        row_str = sprintf('\t');
      end
      
      row_str = cat(2,row_str,sprintf('{'));
      for wf_idx = 1:length(wfs)
        wf = wfs(wf_idx);
        if wf_idx == 1
          row_str = cat(2,row_str,sprintf('['));
        else
          row_str = cat(2,row_str,sprintf(',['));
        end
        imgs_adcs_idx = 1;
        if abs(iq_mode)
          row_str = cat(2,row_str,sprintf('%d %d', iq_mode*wf, imgs_adcs(imgs_adcs_idx)));
          for imgs_adcs_idx = 2:length(imgs_adcs)
            row_str = cat(2,row_str,sprintf('; %d %d', iq_mode*wf, imgs_adcs(imgs_adcs_idx)));
          end
        else
          row_str = cat(2,row_str,sprintf('%d %d', wf, imgs_adcs(imgs_adcs_idx)));
          for imgs_adcs_idx = 2:length(imgs_adcs)
            row_str = cat(2,row_str,sprintf('; %d %d', wf, imgs_adcs(imgs_adcs_idx)));
          end
        end
        row_str = cat(2,row_str,sprintf(']'));
      end
      row_str = cat(2,row_str,sprintf('}'));
      
      fprintf('%s\n',row_str)
      
      %% Write to parameter file if enabled
      if param_file.write_en
        param_idx = settings(set_idx).param_idx;
        
        params(param_idx).get_heights = merge_structs(params(param_idx).get_heights,default.get_heights);
        params(param_idx).csarp = merge_structs(params(param_idx).csarp,default.csarp);
        params(param_idx).combine = merge_structs(params(param_idx).combine,default.combine);
      end
      
    end
  end
  
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
