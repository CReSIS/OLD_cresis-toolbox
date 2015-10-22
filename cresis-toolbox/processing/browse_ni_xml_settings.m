% script browse_ni_xml_settings
%
% Reads in all the XML files from a season (probably will need to modify
% script some to get paths right) and then prints out information on each
% and which data files go with each XML files.
% Skeleton script for debugging purposes.
%
% Author: John Paden

fns = {'D:\20150913\awi\'};
sub_directory = '';

header_load_date = datenum(2015,9,13);

%% MCORDS 5: 2015 Greenland C130
param.radar_name = 'mcords5';
xml_file_prefix = 'mcords5_20150913';
data_file_prefix = 'mcords5';
board_sub_directory = 'chan1';
header_load_func = @basic_load_mcords5;
header_load_params.clk = 1600e6;
adc_bits = 12;
adc_full_scale = 2;
rx_paths = [1:24];
chan_equal_dB = '[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
chan_equal_deg = '[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
chan_equal_Tsys = '[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]/1e9';
Tadc = []; % normally leave empty to use value in file header
Tadc_adjust = 9.92e-6; % leave this empty or set it to zero at first, determine this value later using data over surface with known height or from surface multiple
max_DDS_RAM = 3000;
tx_voltage = sqrt(1000*50)*10^(-2/20);
iq_mode = 0;
fs = 1600e6;
fs_sync = 1.6e9/8;
tx_DDS_mask = [1 1 1 1 1 1 1 1];
imgs_adcs = 1:8;
radar_worksheet_headers = {'Tpd','Tadc','Tadc_adjust','f0','f1','ref_fn','tukey','tx_weights','rx_paths','adc_gains','chan_equal_dB','chan_equal_deg','Tsys','DDC_mode','DDC_freq'};
radar_worksheet_headers_type = {'r','r','r','r','r','r','r','r','r','r','r','r','r','r','r'};

%% MCORDS 4: 2013 Antarctica Basler
% param.radar_name = 'mcords4';
% xml_file_prefix = 'mcords4';
% data_file_prefix = 'mcords4';
% board_sub_directory = 'chan2';
% header_load_func = @basic_load_mcords4;
% header_load_params.clk = 1e9/2;
% adc_bits = 12;
% adc_full_scale = 2;
% rx_paths = [1 3 5 7 2 4 6 8];
% chan_equal_dB = '[-0.9 0 0.3 0.4 0.2 0.4 0.2 -1]';
% chan_equal_deg = '[12.6 19.3 27.3 35.2 28.5 28 9.4 6.3]';
% chan_equal_Tsys = '[82.06 82.00 81.88 82.10 82.08 81.97 81.90 81.77]/1e9';
% Tadc_adjust = 0; % leave this empty or set it to zero first, you can later determine this value from surface multiple.
% max_DDS_RAM = 60000;
% tx_voltage = sqrt(250*50)*10^(-2/20);
% iq_mode = -j;
% fs = 1e9/2;
% fs_sync = 1e9/8;
% tx_DDS_mask = [1 1 1 1 1 1 1 1];
% imgs_adcs = 1:8;
% % radar_worksheet_headers = {'Tpd','Tadc','Tadc_adjust','f0','f1','ref_fn','tukey','tx_weights','rx_paths','adc_gains','chan_equal_dB','chan_equal_deg','Tsys'};
% % radar_worksheet_headers_type = {'r','r','r','r','r','r','r','r','r','r','r','r','r'};
% radar_worksheet_headers = {'Tpd','Tadc','f0','f1','ref_fn','tukey','tx_weights','rx_paths','adc_gains','chan_equal_dB','chan_equal_deg','Tsys'};
% radar_worksheet_headers_type = {'r','r','r','r','r','r','r','r','r','r','r','r'};

%% MCORDS 3: 2014 Greenland P3
% param.radar_name = 'mcords3';
% data_file_prefix = 'mcords3';
% board_sub_directory = 'board0';
% header_load_func = @basic_load_mcords3;
% header_load_params.clk = 1e9/9;
% adc_bits = 14;
% adc_full_scale = 2;
% rx_paths = [1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
% chan_equal_dB = '[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
% chan_equal_deg = '[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
% chan_equal_Tsys = '([0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]+1640)/1e9';
% max_DDS_RAM = 30000;
% tx_voltage = sqrt(300*50)*10^(-2/20);
% iq_mode = 0;
% fs = 1e9/9;
% fs_sync = 1e9/18;
% tx_DDS_mask = [1 1 1 1 1 1 1 0];
% imgs_adcs = 2:15;
% radar_worksheet_headers = {'Tpd','Tadc','Tadc_adjust','f0','f1','ref_fn','tukey','tx_weights','rx_paths','adc_gains','chan_equal_dB','chan_equal_deg','Tsys'};
% radar_worksheet_headers_type = {'r','r','r','r','r','r','r','r','r','r','r','r','r'};

%% MCORDS 3: 2014 Antarctica DC8
% param.radar_name = 'mcords3';
% data_file_prefix = 'mcords3';
% board_sub_directory = 'board0';
% header_load_func = @basic_load_mcords3;
% header_load_params.clk = 300e6/2;
% adc_bits = 14;
% adc_full_scale = 2;
% %rx_paths = [1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
% rx_paths = [1 2 3 4 5 6 1 1];
% chan_equal_dB = '[-0.3 -3.7 -1.9  0 -4.7  0.30]';
% chan_equal_deg = '[-11.4 -33.9 -44.6 0 -39.8 -48.5]';
% chan_equal_Tsys = '[0 0 0 0 0 0]/1e9';
% Tadc_adjust = 1.6e-6; % leave this empty or set it to zero first, you can later determine this value from surface multiple.
% max_DDS_RAM = 30000;
% tx_voltage = sqrt(300*50)*10^(-2/20);
% iq_mode = 0;
% %fs = 1e9/9;
% %fs_sync = 1e9/18;
% fs = 300e6/2;
% fs_sync = 300e6/4;
% tx_DDS_mask = [1 1 1 1 1 1];
% imgs_adcs = [1 2 3 4 5 6];
% radar_worksheet_headers = {'Tpd','Tadc','Tadc_adjust','f0','f1','ref_fn','tukey','tx_weights','rx_paths','adc_gains','chan_equal_dB','chan_equal_deg','Tsys'};
% radar_worksheet_headers_type = {'r','r','r','r','r','r','r','r','r','r','r','r','r'};

%% MCORDS 2
% param.radar_name = 'mcords2';
% data_file_prefix = 'mcords2';
% board_sub_directory = 'board0';
% header_load_func = @basic_load_mcords2;
% header_load_params.clk = 1e9/9;
% adc_bits = 14;
% adc_full_scale = 2;
% rx_paths = [1 1 2 3 4 5 6 7];
% chan_equal_dB = '[0 0 0 0 0 0 0 0]';
% chan_equal_deg = '[0 0 0 0 0 0 0 0]';
% chan_equal_Tsys = '[790 790 790 790 790 790 790 790]/1e9';
% max_DDS_RAM = 30000;
% tx_voltage = sqrt(300*50)*10^(-2/20);
% iq_mode = 0;
% fs = 1e9/9;
% tx_DDS_mask = [1 1 1 1 1 1 1 0];
% imgs_adcs = 2:15;
% radar_worksheet_headers = {'Tpd','Tadc','Tadc_adjust','f0','f1','ref_fn','tukey','tx_weights','rx_paths','adc_gains','chan_equal_dB','chan_equal_deg','Tsys'};
% radar_worksheet_headers_type = {'r','r','r','r','r','r','r','r','r','r','r','r','r'};

%% 
% xml_version = There are three versions of XML files
%   mcords3: pre-2014 Greenland P3 is 1
%   mcords3: post-2014 Greenland P3 is 3
%   mcords4: 2
%   mcords5: 4 (e.g. 2015 LC130, 2015 Polar6)
%   (THIS OFTEN NEEDS TO BE SET)
xml_version = 4;

gps_fn = '/cresis/snfs1/dataproducts/csarp_support/gps/2013_Antarctica_Basler/gps_20131216.mat';
gps_fn = '';

if ~isempty(gps_fn)
  gps = load(gps_fn);
end

MIN_FILES_IN_SEGMENT = 2;
SKIP_FILES_START = 0;
SKIP_FILES_END = 0;

geotiff_fn = '/mnt/products/GIS_data/antarctica/Landsat-7/Antarctica_LIMA_480m.tif';
geotiff_fn = '';

% =======================================================================
% =======================================================================
%% Automated Section
% =======================================================================
% =======================================================================

h_roll_fig = 2;
figure(h_roll_fig); clf;
h_roll_axes = axes('Parent',h_roll_fig);
h_fig = 1;
figure(h_fig); clf;
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

if strcmpi(param.radar_name,'mcords2')
  if xml_version == 1
    config_var = 'Configuration';
    config_var_enc = 'Configuration';
    prf_var = 'PRF_Hz';
    ram_var = 'RAM Taper';
    ram_var_enc = 'RAMZ20Taper';
    xml_file_prefix = 'DDS';
    phase_var = 'Phase_Offset_deg';
    phase_var_enc = 'PhaseZ20OffsetZ20Z28degZ29';
    wave_var_enc = 'Z23Waveforms';
    ttl_start_var_enc = 'TTLZ20StartZ20Z28TTLZ30Z2DZ3ETTLZ37Z29';
    ttl_length_var_enc = 'TTLZ20LengthZ20Z28TTLZ30Z2DZ3ETTLZ37Z29';
    TTL_prog_delay = 650;
  end
  fs = 1e9/9;
  fs_sync = 1e9/18;
elseif strcmpi(param.radar_name,'mcords3')
  if xml_version == 1
    config_var = 'Configuration';
    config_var_enc = 'Configuration';
    prf_var = 'PRF_Hz';
    ram_var = 'RAM Taper';
    ram_var_enc = 'RAMZ20Taper';
    xml_file_prefix = 'DDS';
    phase_var = 'Phase_Offset_deg';
    phase_var_enc = 'PhaseZ20OffsetZ20Z28degZ29';
    wave_var_enc = 'Z23Waveforms';
    ttl_start_var_enc = 'TTLZ20StartZ20Z28TTLZ30Z2DZ3ETTLZ37Z29';
    ttl_length_var_enc = 'TTLZ20LengthZ20Z28TTLZ30Z2DZ3ETTLZ37Z29';
    TTL_prog_delay = 650;
  elseif xml_version == 3
    config_var = 'DDS_Setup';
    config_var_enc = 'DDSZ20Setup';
    prf_var = 'PRF';
    ram_var = 'Ram_Amplitude';
    ram_var_enc = 'RamZ20Amplitude';
    xml_file_prefix = 'radar';
    phase_var = 'Phase_Offset';
    phase_var_enc = 'PhaseZ20Offset';
    wave_var_enc = 'Z23Wave';
    ttl_start_var_enc = 'TTLZ20Start';
    ttl_length_var_enc = 'TTLZ20Length';
    TTL_prog_delay = 650;
  else
    error('Unsupported version');
  end
elseif strcmpi(param.radar_name,'mcords4')
  if xml_version == 2
    config_var = 'DDS_Setup';
    config_var_enc = 'DDSZ20Setup';
    prf_var = 'PRF';
    ram_var = 'Ram_Amplitude';
    ram_var_enc = 'RamZ20Amplitude';
    xml_file_prefix = 'mcords4';
    phase_var = 'Phase_Offset';
    phase_var_enc = 'PhaseZ20Offset';
    wave_var_enc = 'Z23Wave';
    ttl_start_var_enc = 'TTLZ20Start';
    ttl_length_var_enc = 'TTLZ20Length';
    TTL_prog_delay = 650;
  else
    error('Unsupported version');
  end
elseif strcmpi(param.radar_name,'mcords5')
  if any(xml_version == [4])
    config_var = 'DDS_Setup';
    NCO_freq = 'NCO_freq';
    config_var_enc = 'DDSZ20Setup';
    prf_var = 'PRF';
    ram_var = 'Ram_Amplitude';
    ram_var_enc = 'RamZ20Amplitude';
    phase_var = 'Phase_Offset';
    phase_var_enc = 'PhaseZ20Offset';
    wave_var_enc = 'Z23Wave';
    ttl_start_var_enc = 'TTLZ20Start';
    ttl_length_var_enc = 'TTLZ20Length';
    TTL_prog_delay = 650;
  else
    error('Unsupported version');
  end
end

max_wfs = 0;
for fn_idx = 1:length(fns);
  fn = fns{fn_idx};
  [~,fn_name] = fileparts(fn);
  fprintf('Processing %s (%d of %d)\n', fn_name, fn_idx, length(fns));
  settings = read_ni_xml_directory(fullfile(fn,sub_directory),xml_file_prefix,false);
  data_fns = get_filenames(fullfile(fn,sub_directory,board_sub_directory),data_file_prefix,'','.bin',struct('recursive',1));
  fn_datenums = [];
  for data_fn_idx = 1:length(data_fns)
    fname = fname_info_mcords2(data_fns{data_fn_idx});
    fn_datenums(end+1) = fname.datenum;
  end
  for set_idx = 1:length(settings)
    fprintf('  ===================== Setting %d =================', set_idx);
    if isfield(settings(set_idx),'Configuration')
      fprintf('  %s: %d waveforms\n', settings(set_idx).fn, length(settings(set_idx).Configuration.Waveforms));
      if max_wfs < length(settings(set_idx).Configuration.Waveforms)
        max_wfs = length(settings(set_idx).Configuration.Waveforms);
      end
    else
      fprintf('  %s: %d waveforms\n', settings(set_idx).fn, length(settings(set_idx).(config_var).Waveforms));
      if max_wfs < length(settings(set_idx).(config_var).Waveforms)
        max_wfs = length(settings(set_idx).(config_var).Waveforms);
      end
      fprintf('PRF:'); fprintf(' %g', settings(set_idx).(config_var).(prf_var)); fprintf('\n');
      fprintf('Amp:'); fprintf(' %g', settings(set_idx).(config_var).Ram_Amplitude); fprintf('\n');
      fprintf('Tukey:'); fprintf(' %g', settings(set_idx).(config_var).RAM_Taper); fprintf('\n');
      fprintf('Stop:'); fprintf(' %g', settings(set_idx).(config_var).Waveforms(1).Stop_Freq); fprintf('\n');
      fprintf('Tx Mask:'); fprintf(' %g', settings(set_idx).(config_var).Waveforms(1).TX_Mask); fprintf('\n');
      for wf = 1:2:length(settings(set_idx).(config_var).Waveforms)
        fprintf('WF %d Atten:', wf); fprintf(' %g', settings(set_idx).(config_var).Waveforms(wf).Attenuator_2); fprintf('\n');
        fprintf('WF %d Len:', wf); fprintf(' %g', settings(set_idx).(config_var).Waveforms(wf).Len_Mult); fprintf('\n');
      end
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
        h_plots(fn_idx,set_idx) = plot(x/1e3,y/1e3,'.','Parent',h_axes);
        hdr_roll = interp1(gps.gps_time,gps.roll,hdr_gps_time);
        h_roll_plots(fn_idx,set_idx) = plot(hdr_gps_time,hdr_roll,'Parent',h_roll_axes);
      end
      h_texts(fn_idx,set_idx) = text(x(1)/1e3,y(1)/1e3,sprintf('%d:%d',fn_idx,set_idx),'Parent',h_axes);
      h_roll_texts(fn_idx,set_idx) = text(hdr_gps_time(1),hdr_roll(1),sprintf('%d:%d',fn_idx,set_idx),'Parent',h_roll_axes);
    end
  end
end

%% Vectors Worksheet
fprintf('Vectors and Records WorkSheet\n');
if any(strcmpi(param.radar_name,{'mcords3','mcords4'}))
  fprintf('%s\t%s\t%s\t%s\n', ...
    'Date','1','file.start_idx','file.stop_idx');
  fprintf('%s\t%s\t%s\t%s\n', ...
    'YYYYMMDD','Segment','r','r');
elseif any(strcmpi(param.radar_name,{'mcords5'}))
  fprintf('%s\t%s\t%s\t%s\t%s\n', ...
    'Date','1','file.start_idx','file.stop_idx','file_version');
  fprintf('%s\t%s\t%s\t%s\t%s\n', ...
    'YYYYMMDD','Segment','r','r','r');
end
for set_idx = 1:length(settings)
  if settings(set_idx).match_idxs(end) - settings(set_idx).match_idxs(1) >= MIN_FILES_IN_SEGMENT-1
    if strcmpi(param.radar_name,'mcords3')
      fprintf('%s\t%02d\t%d\t%d\n',datestr(header_load_date,'YYYYmmDD'),set_idx, ...
        settings(set_idx).match_idxs(1)+1, settings(set_idx).match_idxs(end))
    elseif any(strcmpi(param.radar_name,{'mcords4'}))
      fprintf('%s\t%02d\t%d\t%d\n',datestr(header_load_date,'YYYYmmDD'),set_idx, ...
        settings(set_idx).match_idxs(1), settings(set_idx).match_idxs(end)-1)
    elseif any(strcmpi(param.radar_name,{'mcords5'}))
      if settings(set_idx).DDC_Ctrl.DDC_sel.Val == 0
        file_version = 408;
      elseif any(settings(set_idx).DDC_Ctrl.DDC_sel.Val == [1 2])
        file_version = 407;
      else
        error('Invalid DDC setting.');
      end
      fprintf('%s\t%02d\t%d\t%d\t%d\n',datestr(header_load_date,'YYYYmmDD'),set_idx, ...
        settings(set_idx).match_idxs(1)+SKIP_FILES_START, settings(set_idx).match_idxs(end)-SKIP_FILES_END, ...
        file_version);
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
    row_str = sprintf('%s\t%02d\t%.16g\t%g\t%d\t%g\t%d',datestr(header_load_date,'YYYYmmDD'),set_idx, ...
      fs, settings(set_idx).(config_var).(prf_var), adc_bits, adc_full_scale, length(settings(set_idx).(config_var).Waveforms));
    for wf = 1:length(settings(set_idx).(config_var).Waveforms)

      if any(strcmpi('Tpd',radar_worksheet_headers))
        row_str = cat(2,row_str, sprintf('\t%g',double(settings(set_idx).(config_var).Waveforms(wf).Len_Mult)*settings(set_idx).(config_var).Base_Len));
      end
      if any(strcmpi('Tadc',radar_worksheet_headers))
        row_str = cat(2,row_str, sprintf('\t%g',Tadc));
      end
      if any(strcmpi('Tadc_adjust',radar_worksheet_headers))
        row_str = cat(2,row_str, sprintf('\t%g',Tadc_adjust));
      end
      if any(strcmpi('f0',radar_worksheet_headers))
        row_str = cat(2,row_str, sprintf('\t%.16g',settings(set_idx).(config_var).Waveforms(wf).Start_Freq(1)));
      end
      if any(strcmpi('f1',radar_worksheet_headers))
        row_str = cat(2,row_str, sprintf('\t%.16g',settings(set_idx).(config_var).Waveforms(wf).Stop_Freq(1)));
      end
      if any(strcmpi('ref_fn',radar_worksheet_headers))
        row_str = cat(2,row_str, sprintf('\t'));
      end
      if any(strcmpi('tukey',radar_worksheet_headers))
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
      row_str = cat(2,row_str, ...
        sprintf('\t[%g', tx_weights(1)));
      row_str = cat(2,row_str, ...
        sprintf(' %g', tx_weights(2:end)));
      row_str = cat(2,row_str, ...
        sprintf(']'));
      % Rx paths
      row_str = cat(2,row_str, ...
        sprintf('\t[%d', rx_paths(1)));
      row_str = cat(2,row_str, ...
        sprintf(' %d', rx_paths(2:end)));
      row_str = cat(2,row_str, ...
        sprintf(']'));
      % ADC Gains
      atten = double(settings(set_idx).(config_var).Waveforms(wf).Attenuator_1(1)) ...
        + double(settings(set_idx).(config_var).Waveforms(wf).Attenuator_2(1));
      if 0
        atten = atten * ones(size(rx_paths));
        row_str = cat(2,row_str, ...
          sprintf('\t10.^((52 - [%g', atten(1)));
        row_str = cat(2,row_str, ...
          sprintf(' %g', atten(2:end)));
        row_str = cat(2,row_str, ...
          sprintf('])/20)'));
      else
        row_str = cat(2,row_str, ...
          sprintf('\t10.^((52 - %g*ones(1,%d))/20)', atten(1), length(rx_paths)));
      end
      % Chan Equal DB
      row_str = cat(2,row_str, sprintf('\t%s',chan_equal_dB));
      % Chan Equal deg
      row_str = cat(2,row_str, sprintf('\t%s',chan_equal_deg));
      % Tsys
      row_str = cat(2,row_str, sprintf('\t%s',chan_equal_Tsys));
      if any(strcmpi('DDC_mode',radar_worksheet_headers))
        row_str = cat(2,row_str, sprintf('\t%d',settings(set_idx).DDC_Ctrl.DDC_sel.Val));
      end
      if any(strcmpi('DDC_freq',radar_worksheet_headers))
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
  end
end

return;
