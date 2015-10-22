function success = analysis_task(param)
% success = analysis_task(param)
%
% Cluster task for analysis.  Does the actual data loading and
% generates analysis results.
%
% param = struct controlling the loading, processing, evaluating of psd and
%   calculating of noise power.
%   .profile = structure containg basepath to analysis outputs
%       .out_path
%   .radar_name = string indentifying radar (i.e. "mcords" or "mcords2")
%   .season_name = string identifying season with format
%    "YYYY_Country_Platform" (i.e. "2011_Greenalnd_TO")
%   .day_seg = string identifying day and segment with format "YYYYMMDD_##"
%    (i.e. "20110401_01")
%   .radar = structure containing fields from radar configuration file
%       .rx_path
%           .chan_equal
%           .td
%       .fs = sampling frequency (Hz)
%       .prf = pulses per second (Hz)
%       .sample_size
%       .Vpp_scale
%       .wfs
%           .Tpd = pulse length (seconds). duration of transmit chirp
%           .t0  =
%           .f0  = start frequency of LFM transmit chirp (Hz)
%           .f1  = stop frequency of LFM transmit chirp (Hz)
%           .ref_fn
%           .tukey
%           .blank
%           .tx_weights = [1 x M] vector (where M is the number of adcs)
%           specifying transmit weights for a particular waveform.
%           .rx_paths = [1 x M] vector (where M is the number of adcs)
%           that associates adcs and antenna elements for each transmit
%           waveform.
%           .adc_gains = [1 x M] vector (where M is the number of adcs)
%           that specifes adc gain settings of each transmit waveform.
%   .load = struct specifying which records to load
%       .records_fn = filename of records file
%       .recs = [1 x 2] array containing current records to be loaded
%       .imgs = cell vector of image to load. The image is 1x2 array,
%       wf/adc pair.
%       NOTE: wfs/adc pairs are not indices into anything, they are the
%       absolute waveform/adc numbers.  The records file will be loaded to
%       decode which index each wf/adc belongs to.
%   .debug_level = debug level (scalar integer)
%
%
%   .analysis = structure controlling analysis processing
%       .file
%           .base_dir = string specifying path to raw data stored on server
%           (i.e. '/cresis/data4/MCoRDS/2011_Greenland_TO/20110401')
%           .adc_folder_name = string adc folder name (i.e.
%           '/chan%d/seg##/')
%           .file_prefix
%       .analysis_type = string that specifies type of analysis output
%       ('np' or 'psd')
%       .gps
%           .en = binary value that enables analysis to use gps_time
%           when available.  If no gps is available (which is often the
%           case with runway data or 50 ohm data) this value must be set to
%           0.
%       .out_path
%       .records = [1 x 2] array specifies record start_idx and stop_idx of
%       data being loaded inside task.
%       .ft_bins = [K x 2] array (where K is total number of image pairs
%       used in analysis)  that specifies start and stop bins of fast
%       time window isolating samples containing interference and
%       thermal noise.  To estimate interference spectrum of real data,
%       ft_bins should be chosen to isolate samples of post bed echo data.
%       .block_size = [1 x 2] array that specifies the processing block
%       size in the slow-time dimension.  Since the 'np' type analysis
%       loads data without coherent integrations, this number should be
%       made small (i.e. block_size < 1000).
%       .psd
%           .coh_ave = number coherent integrations used to load data
%       .base_dir = string specifying CSARP_analysis/day_seg directory for
%       analysis outputs (i.e.
%       '/cresis/scratch2/mdce/mcords/2011_Greenland_TO/CSARP_analysis/ ...
%       20110401_01')
%
%
%
%
% analysis fields used by load_mcords_data.m (see that function for
% details) include:
%   .ft_wind
%   .ft_wind_time
%   .trim_vals
%   .pulse_rfi.en
%   .pulse_rfi.inc_ave
%   .pulse_rfi.thresh_scale
%   .radar
%   .Vpp_scale = radar Vpp for full scale quantization
%   .rxs = struct array of receiver equalization coefficients
%       .wfs = struct array for each waveform
%           .chan_equal = scalar complex double (data DIVIDED by this)
%           .td = time delay correction (applied during pulse compression)
%
% Author: Theresa Stumpf
%
% See also analysis.m
global g_data;

physical_constants;

wf          = param.load.imgs{1}(1,1);
adc         = param.load.imgs{1}(1,2);
recs        = param.load.recs;

ct_filename_param               = param.gRadar;
ct_filename_param.radar_name    = param.radar_name;
ct_filename_param.season_name   = param.season_name;
ct_filename_param.day_seg       = param.day_seg;
param.load.records_fn = ct_filename_support(ct_filename_param, ...
  param.load.records_fn,'records');

% =====================================================================
% Determine which records to load with load_mcords_data
%
% Load records on either side of the current block, note if at the
% beginning or end of the segment.  Load with minimal presumming.

load_param.load.recs(1) = param.load.recs(1);
load_param.load.recs(2) = param.load.recs(2);
[records, old_param_records ]= ...
  read_records_aux_files(param.load.records_fn,load_param.load.recs);

% =====================================================================
% Collect waveform information into one structure
%  (used by load_mcords_data)
param.analysis.ft_dec = 0;
param.analysis.ft_wind_time = 0;
param.analysis.ft_wind = '';
param.analysis.pulse_comp = 0;
param.analysis.combine_rx = 0;
param.analysis.trim_vals  = [1 1];
param.analysis.pulse_rfi.en = 1;
param.analysis.pulse_rfi.inc_ave = 101;
param.analysis.pulse_rfi.thresh_scale = 10^(13/10);

[wfs,rec_data_size] = load_mcords_wfs(records.wfs, param, ...
  1:max(old_param_records.file.adcs), param.analysis);
load_param.wfs = wfs;
load_param.load.rec_data_size = rec_data_size;

% =====================================================================
% Collect record file information required for using load_mcords_data
%  - Performs mapping between param.rxs and the records file contents
%  - Translates filenames from relative to absolute
%  - Makes filenames have the correct filesep

% Create a list of unique adcs required by the imgs list
param.load.adcs = [];

for idx = 1:length(param.load.imgs)
  param.load.adcs = unique(cat(2, ...
    abs(param.load.imgs{idx}(:,2)).', param.load.adcs));
end

% recs = param.load.recs - param.load.recs(1) + 1;
for idx = 1:length(param.load.adcs)
  adc_idx = find(old_param_records.file.adcs == param.load.adcs(idx));
  if isempty(adc_idx)
    error('ADC %d not present in records file\n', param.load.adcs(idx));
  end
  
  % Just get the file-information for the records we need
  load_param.load.file_idx{idx} = records.file_idx{adc_idx};
  load_param.load.offset{idx} = records.offset{adc_idx};
  file_idxs = unique(load_param.load.file_idx{idx});
  % Recognize if first record is really from previous file and it is a
  % valid record (i.e. offset does not equal -2^31)
  if sign(load_param.load.offset{idx}(1)) < 0 && load_param.load.offset{idx}(1) ~= -2^31
    file_idxs = [file_idxs(1)-1 file_idxs];
  end
  % Just copy the filenames we need
  load_param.load.filenames{idx}(file_idxs) = records.filenames{adc_idx}(file_idxs);
  
  base_dir = ct_filename_data(ct_filename_param,param.analysis.file.base_dir);
  adc_folder_name = param.analysis.file.adc_folder_name;
  if isfield(param.analysis.file,'file_prefix')
    file_prefix = param.analysis.file.file_prefix;
  else
    file_prefix = '';
  end
  
  % Create sub-folder name for the particular adc
  adc_idx_insert_idxs = strfind(adc_folder_name,'%d');
  board_idx_insert_idxs = strfind(adc_folder_name,'%b');
  all_idx_insert_idxs = sort([adc_idx_insert_idxs board_idx_insert_idxs]);
  mat_cmd = 'adc_folder_name = sprintf(adc_folder_name';
  for all_idx_insert_idx = all_idx_insert_idxs
    if any(all_idx_insert_idx == adc_idx_insert_idxs)
      % Insert adc number
      mat_cmd = [mat_cmd sprintf(', %d',param.load.adcs(idx))];
    else
      % Insert board number
      mat_cmd = [mat_cmd sprintf(', %d',floor((param.load.adcs(idx)-1)/4) )];
      adc_folder_name(all_idx_insert_idx+1) = 'd';
    end
  end
  mat_cmd = [mat_cmd ');'];
  eval(mat_cmd);
  
  filepath = fullfile(base_dir, adc_folder_name);
  
  % Convert relative file paths into absolute file paths if required,
  % also corrects filesep (\ and /)
  for file_idx = 1:length(load_param.load.filenames{idx})
    load_param.load.filenames{idx}{file_idx} ...
      = fullfile(filepath,load_param.load.filenames{idx}{file_idx});
  end
end

% =====================================================================
% Setup control parameters for load_mcords_data
load_param.proc.combine_rx        = param.analysis.combine_rx;
load_param.proc.pulse_comp        = param.analysis.pulse_comp;
load_param.proc.ft_dec            = param.analysis.ft_dec;
load_param.proc.ft_wind           = param.analysis.ft_wind;
load_param.proc.ft_wind_time      = param.ft_wind_time;
load_param.proc.pulse_rfi         = param.analysis.pulse_rfi;
load_param.proc.trim_vals         = param.analysis.trim_vals;
load_param.radar                  = param.radar;
load_param.load.imgs              = param.load.imgs;
load_param.load.adcs              = param.load.adcs;

if strcmpi(param.analysis.analysis_type,'psd')
  load_param.proc.presums         = param.analysis.psd.coh_ave;
elseif strcmpi(param.analysis.analysis_type,'np')
  load_param.proc.presums         = 1;
end

% =========================================================================
% Determine identifier string to tag blocks in psd and/or np filenames.
% When GPS is available, assigns GPS time (in HHssmsms as identifier
start_stamp   = sprintf('%s_%010.0f',param.day_seg,load_param.load.recs(1));

% =====================================================================
% Load and process the image pair specified outside loop
%
% For each waveform, adc combination:
% 1. Load receiver data separately (minimal presumming)
% 2.
% =====================================================================
% Load data into g_data using load_mcords_data
if strcmpi(param.radar_name,'mcords')
  load_mcords_data(load_param);
elseif strcmpi(param.radar_name,'mcords2')
  load_mcords2_data(load_param);
end
g_data = g_data{1};

bins = param.analysis.ft_bins;
start_bin = bins(1);
if bins(2) >= size(g_data,1);
  stop_bin = size(g_data,1) - 1;
else
  stop_bin = bins(2);
end
bins = [start_bin stop_bin];
g_data = g_data(bins(1):bins(2),:,:);


% PSD Analysis
% -------------------------------------------------------------------------
if strcmpi(param.analysis.analysis_type,'psd')
  
  psd_out_dir = fullfile(param.analysis.base_dir, ...
    sprintf('psd_analysis/recs_%010d_%010d/bins_%04d_%04d/coh_ave_%05d/wf_%02d/', ...
    recs(1),recs(2),bins(1),bins(2),param.analysis.psd.coh_ave,wf));
  
  if ~exist(psd_out_dir,'dir')
    mkdir(psd_out_dir);
  end
  
  % Evaluate power spectral density and frequency axis
  Nt = size(g_data,1);
  fs = wfs(wf).fs;
  df = 1/(Nt*wfs(wf).dt);
  
  f0 = wfs(wf).f0;
  f1 = wfs(wf).f1;
  fc = (f0 + f1)/2;
  freq = fs*floor(fc/fs) + (0:df:(Nt-1)*df).';
  
  dt = wfs(wf).dt;
  Tpsd = Nt*dt;
  
  g_data = fft(g_data(:,:),Nt);         % Frequency spectrum of raw voltages
  g_data = abs(g_data).^2;
  g_data = (1/Nt).*mean(g_data,2);      % Periodogram wrt 1 ohm (Watts/bin)
  
  Data = lp(g_data) - lp(fs);           % Estimated one sided PSD wrt 1 ohm
  % (dBW/Hz)
  
  % Build param_analyze
  param_radar                       = param.radar;
  param_radar.wfs                   = wfs;
  param_radar.radar_name            = param.radar_name;
  param_analysis                    = param.analysis;
  param_analysis.img                = [wf adc];
  param_analysis.recs               = load_param.load.recs;
  param_analysis.psd.Nt             = Nt;
  param_analysis.psd.f0             = f0;
  param_analysis.psd.f1             = f1;
  param_analysis.psd.fc             = fc;
  param_analysis.psd.df             = df;
  param_analysis.psd.Tpsd           = Tpsd;
  param_analysis.season_name        = param.season_name;
  param_analysis.day_seg            = param.day_seg;
  
  
  psd_fn = fullfile(psd_out_dir,sprintf('adc_%02d_element_%02d_%s.mat', ...
    adc, wfs(wf).rx_paths(adc),start_stamp));
  fprintf('  Saving PSD to file %s\n', psd_fn);
  save(psd_fn,'Data','freq','param_analysis','param_radar')
  
  
  
  % Noise Power Analysis
  % NOTE:  Noise power calculation is not correct!!  Measured noise powers
  % are less than the expected power of the thermal noise (~97 dBm)
  % -----------------------------------------------------------------------
elseif strcmpi(param.analysis.analysis_type,'np')
  
  tmp_wf_dir = fullfile(param.analysis.base_dir,...
    sprintf('np_analysis/recs_%010d_%010d/tmp/wf_%02d/',param.analysis.records(1),param.analysis.records(2),wf));
  if ~exist(tmp_wf_dir,'dir')
    mkdir(tmp_wf_dir);
  end
  
  %   When associated gps is availabe, create a gps time variable to save
  %   with noise power
  if param.analysis.gps.en
    time_fn = fullfile(tmp_wf_dir,sprintf('np_gps_time_%010.0f.mat',records.gps_time(1)));
    gps_time = records.gps_time;
    save(time_fn,'gps_time');
  end
  
  % Calculate expected noise power
  BW = wfs(wf).f1 - wfs(wf).f0;
  %   nf = 10^(2/10);
  
  nf = 10^(1.91/10);
  exp_np_dBm = 10*log10(BoltzmannConst*290*BW*nf) - 10*log10(10^-3);
  
  % Evaluate noise power of each range line in data array
  presums = wfs(wf).presums;
  meas_power = (abs(g_data(:,:)).^2);
  ave_power = mean(meas_power,1);
  meas_np_dBm = lp(ave_power) - lp(10^-3) - lp(50);
  
  % Save result
  np_fn = fullfile(tmp_wf_dir, ...
    sprintf('adc_%02d_element_%02d_%s.mat', ...
    adc, wfs(wf).rx_paths(adc),start_stamp));
  fprintf('  Saving tmp noise power to file %s\n', np_fn);
  save(np_fn,'meas_np_dBm','exp_np_dBm','presums');
end

success = true;

return




