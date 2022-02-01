function [ param, hdr, frames, records, exec_good ] = flightline_extract(param)

hdr = [];
frames = [];
records = [];
% gps = [];
exec_good = 0;

c = physical_constants('c');

%% PARAMS
fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.sim.day_seg, datestr(now));
fprintf('=====================================================================\n');

params = read_param_xls(ct_filename_param(param,sprintf('%s_param_%s.xls',param.sim.radar_name,param.sim.season_name)));
found_day_seg = 0;
for param_idx = 1:length(params)
  if strcmpi( params(param_idx).day_seg , param.sim.day_seg)
    %     param = merge_structs( param, params(param_idx) );
    param_extracted = params(param_idx);
    found_day_seg = 1;
    break;
  end
end

if ~found_day_seg; fprintf('day_seg not found\n'); return; end;
clear params found_day_seg;

%% Frames, Records, GPS, Layer
fprintf('Loading frames (%s)\n',  datestr(now));
frames = frames_load( param_extracted );

fprintf('Loading records (%s)\n', datestr(now));
records = records_load( param_extracted );

% gps_fn = ct_filename_support(param_extracted,'','gps',1);
% fprintf('Loading gps %s (%s)\n', gps_fn, datestr(now));
% gps    = load( gps_fn );

fprintf('Loading layer (%s)\n', datestr(now));
idx = 1;
layer_params(idx).name = 'surface';
if 0
  layer_params(idx).source = 'echogram';
  layer_params(idx).echogram_source = 'qlook';
else
  layer_params(idx).source = 'layerdata';
end
[layers, new_layer_params] = opsLoadLayers(param.sim,layer_params);


if isfield(param.sim,'start_gps_time') && isfield(param.sim,'stop_gps_time')
  start_idx = find(records.gps_time>=param.sim.start_gps_time,1,'first');
  stop_idx  = find(records.gps_time>=param.sim.stop_gps_time,1,'first');
elseif isfield(param.sim,'frame_idx') % select frame
  start_idx   = frames.frame_idxs(param.sim.frame_idx);
  stop_idx    = frames.frame_idxs(param.sim.frame_idx+1);
elseif 0  % just some default values
  start_idx = 28501;
  stop_idx  = 29200;
elseif 1
  start_idx = 1;
  stop_idx  = 747;
end

% Frames (THESE ARE TEMPORARY DEFAULT VALUES)
frames.frame_idxs = 1;
frames.nyquist_zone = NaN; %1
frames.proc_mode = 0;
frames.quality = 1;

% Records
try % truncates records' structure (instead of read_records_aux_files)
  rec_len = length(records.gps_time);
  records = struct_truncate(records,rec_len,start_idx,stop_idx,0);
catch % fail-safe
  records.gps_time = records.gps_time(start_idx:stop_idx);
  records.lat = records.lat(start_idx:stop_idx);
  records.lon = records.lon(start_idx:stop_idx);
  records.elev = records.elev(start_idx:stop_idx);
  records.roll = records.roll(start_idx:stop_idx);
  records.pitch = records.pitch(start_idx:stop_idx);
  records.heading = records.heading(start_idx:stop_idx);
  %   records.surface = records.surface(start_idx:stop_idx);
  % records.offset = records.offset(start_idx:stop_idx);
  % records.raw.epri = records.raw.epri
  % records.raw.seconds = records.raw.seconds
  % records.raw.fraction = records.raw.fraction
  % records.settings.nyquist_zone = records.settings.nyquist_zone
  % records.settings.nyquist_zone_hw = records.settings.nyquist_zone_hw
  % records.settings.records_mask = records.settings.records_mask
end

% clear frames_fn records_fn gps_fn rec_len %start_idx stop_idx

% raw: load_data
param_extracted.load_data.recs                            = [start_idx stop_idx];
param_extracted.load_data.imgs                            = param.sim.imgs;
param_extracted.load_data.pulse_comp                      = false;
param_extracted.load_data.raw_data                        = true;
param_extracted.load_data.ft_wind                         = @hanning;
param_extracted.load_data.combine_rx                      = param.sim.rx_combine;
param_extracted.radar.wfs(param.sim.wfs).coh_noise_method = '';           %<#######
param_extracted.radar.wfs(param.sim.wfs).rx_paths         = param.sim.rxpath;
[param_extracted.hdr,loaded_data] = load_data(param_extracted);

% To use in data_load.m for file
hdr.bad_rec               = param_extracted.hdr.bad_rec;
hdr.nyquist_zone_hw       = param_extracted.hdr.nyquist_zone_hw;
hdr.nyquist_zone_signal   = param_extracted.hdr.nyquist_zone_signal;
hdr.DDC_dec               = param_extracted.hdr.DDC_dec;
hdr.DDC_freq              = param_extracted.hdr.DDC_freq;
hdr.Nt                    = param_extracted.hdr.Nt ;
hdr.t0_raw                = param_extracted.hdr.t0_raw;
hdr.t_ref                 = param_extracted.hdr.t_ref;

param.load_data.recs = [1 stop_idx-start_idx+1];
param.load_data.imgs = param_extracted.load_data.imgs;

% param_extracted.hdr = [];
% loaded_data=[];
%%
% Waveform
param.signal = [];

% Trajectory
param.gps_source = records.gps_source; % can be used for leverarm

param.gps = [];
param.gps.gps_time = records.gps_time;
param.gps.lat = param_extracted.hdr.records{1, 1}.lat;
param.gps.lon = param_extracted.hdr.records{1, 1}.lon;
param.gps.elev = param_extracted.hdr.records{1, 1}.elev;
param.gps.roll = param_extracted.hdr.records{1, 1}.roll;
param.gps.pitch = param_extracted.hdr.records{1, 1}.pitch;
param.gps.heading = param_extracted.hdr.records{1, 1}.heading;
[param.gps.x, param.gps.y, param.gps.z] = geodetic2ecef(wgs84Ellipsoid,param.gps.lat,param.gps.lon,param.gps.elev);

% Target
param.target = [];

% layer points in specified indices (start,stop)


if ~isfield(param.target,'type')
  param.target.type = 'point'; % 'surface'
end

switch param.target.type
  case 'point'
    mid_idx = ceil(length(records.gps_time)/2);
    param.target.lat  = param_extracted.hdr.records{1, 1}.lat(mid_idx);
    param.target.lon  = param_extracted.hdr.records{1, 1}.lon(mid_idx);
    if isnan(records.elev(mid_idx));
      param.target.elev = 0;
    elseif 0 % outdated
      param.target.elev = param_extracted.hdr.records{1, 1}.elev(mid_idx) - c*records.surface(mid_idx)/2;
    elseif 0 % interpolate target(lat lon elev from layers) for gps_time from records
      if 0
        figure; hold on;
        plot(records.gps_time,records.gps_time,'x');
        plot(layers.gps_time,layers.gps_time,'o');
      end
      if 0 % !!!!!!!! plot lat lon for (records, param_extracter.hdr.records, layers)
        figure; hold on;
        plot(records.gps_time,records.gps_time,'x');
        plot(layers.gps_time,layers.gps_time,'o');
      end
      param.target.elev = param_extracted.hdr.records{1, 1}.elev(mid_idx) - c*param_extracted.hdr.records{1, 1}.twtt(mid_idx)/2;
    elseif 1 % temporary override
      param.target.elev = 0;
    end
    
  case 'layer'
    param.target.lat  = param_extracted.hdr.records{1, 1}.lat;
    param.target.lon  = param_extracted.hdr.records{1, 1}.lon;
    param.target.elev = param_extracted.hdr.records{1, 1}.elev - 664;%3e8*records.surface/2;
end
[param.target.x, param.target.y, param.target.z] = geodetic2ecef(wgs84Ellipsoid,param.target.lat,param.target.lon,param.target.elev);

%%

Ntx = 1;
Nrx = 1;
wfs = struct();

% To load wfs
param_extracted.load.imgs = param.sim.imgs;
load_wfs = data_load_wfs(param_extracted,records);

for idx = 1:Ntx
  param.signal.tx(Ntx).gain           = 1; %[]; % Nt x N_elev x N_azi
  param.signal.tx(Ntx).freq           = 1; %[]; % Nt x 1
  param.signal.tx(Ntx).elev_angle     = 1; %[]; % N_elev x 1
  param.signal.tx(Ntx).azimuth_angle  = 1; %[]; % N_azi x 1
end

for idx = 1:Nrx
  param.signal.rx(Nrx).gain           = 1; %[]; % Nt x N_elev x N_azi
  param.signal.rx(Nrx).freq           = 1; %[]; % Nt x 1
  param.signal.rx(Nrx).elev_angle     = 1; %[]; % N_elev x 1
  param.signal.rx(Nrx).azimuth_angle  = 1; %[]; % N_azi x 1
end

% Populate param.radar
param.radar.fs            = param_extracted.radar.fs;
param.radar.prf           = param_extracted.radar.prf;
param.radar.adc_bits      = param_extracted.radar.adc_bits;
param.radar.Vpp_scale     = param_extracted.radar.Vpp_scale;
param.radar.lever_arm_fh  = param_extracted.radar.lever_arm_fh;

% dt = 1/param.radar.fs;
% pri = 1/param.radar.prf;

% for each waveform
for wf_idx = 1:length(param.sim.wfs)
  wf = param.sim.wfs(wf_idx);
  
  %   wfs(wf).Nt    = load_wfs(wf).Nt;
  wfs(wf).Nt                = [];
  wfs(wf).tx_weights        = load_wfs(wf).tx_weights; % Ntx x 1
  wfs(wf).tukey             = 0.1;%load_wfs(wf).tukey;
  wfs(wf).f0                = load_wfs(wf).f0;
  wfs(wf).f1                = load_wfs(wf).f1;
  wfs(wf).BW_window         = load_wfs(wf).BW_window(1:2);
  wfs(wf).fc                = load_wfs(wf).fc;
  wfs(wf).chirp_rate        = load_wfs(wf).chirp_rate;
  wfs(wf).Tpd               = load_wfs(wf).Tpd;
  wfs(wf).tx_weights        = load_wfs(wf).tx_weights;
  wfs(wf).adc_gains_dB      = load_wfs(wf).adc_gains_dB;
  hdr.t0_raw{1} = hdr.t0_raw{1} - hdr.t0_raw{1}(1);
  wfs(wf).t0_raw            = hdr.t0_raw{1}(1); %load_wfs(wf).t0_raw;  %################### hdr.t0_raw{1}(1) ???
  wfs(wf).t_ref             = hdr.t_ref{1}(1); % load_wfs(wf).t_ref; %################### hdr.t_ref ???
  wfs(wf).Tadc_adjust       = load_wfs(wf).Tadc_adjust;
  wfs(wf).DDC_NCO_delay     = load_wfs(wf).DDC_NCO_delay;
  wfs(wf).prepulse_H        = load_wfs(wf).prepulse_H;
  wfs(wf).coh_noise_method  = load_wfs(wf).coh_noise_method;
  
  [output_dir,radar_type,radar_name] = ct_output_dir(param.sim.radar_name);
  
  if strcmpi(radar_type,'deramp')
    wfs(wf).deramp = 1; % Default 0
    
    %     if wfs(wf).Nt == 0 && ~load_wfs(wf).Nt_raw == 0
    %       wfs(wf).Nt = load_wfs(wf).Nt_raw;
    %       fprintf('deramp: Custom Nt from Nt_raw\n');
    %     elseif load_wfs(wf).Nt_raw == 0;
    %       wfs(wf).Nt = size(data{wf,1},1);
    %       fprintf('deramp: Custom Nt from load_data\n');
    %     else
    %       wfs(wf).Nt = 31400; % temporary fix for time axis
    %       fprintf('deramp: Custom Nt = 31400\n');
    %     end
    
    if isfield(hdr,'Nt') && ~isempty(hdr.Nt{1})
      fprintf('deramp: Custom Nt from param_extracted.hdr.Nt{1}(1)');
      wfs(wf).Nt = hdr.Nt{1}(1);
    end
    if isempty(wfs(wf).Nt) || wfs(wf).Nt == 0
      fprintf(' ### Missing Nt in hdr from load_data\n');
      fprintf('deramp: Custom Nt from size(data{wf,1},1)');
      if ~isempty(loaded_data)
        wfs(wf).Nt = size(loaded_data{wf,1},1);
      end
    end
    if isempty(wfs(wf).Nt) || wfs(wf).Nt == 0
      fprintf(' ### Missed it in data from load_data\n');
      fprintf('deramp: Custom Nt from load_wfs(wf).Nt');
      wfs(wf).Nt = load_wfs(wf).Nt;
    end
    if isempty(wfs(wf).Nt) || wfs(wf).Nt == 0
      fprintf(' ### Missing Nt in wfs from data_load_wfs\n');
      fprintf('deramp: Custom Nt from load_wfs(wf).Nt_raw');
      wfs(wf).Nt = load_wfs(wf).Nt_raw;
    end
    if isempty(wfs(wf).Nt) || wfs(wf).Nt == 0
      fprintf(' ### Missing Nt_raw in wfs from data_load_wfs\n');
      fprintf('deramp: Custom Nt = 31400');
      wfs(wf).Nt = 31400; % temporary fix for time axis
    end
    fprintf('\n');
    
    if isfield(load_wfs(wf),'time')
      wfs(wf).time   = load_wfs(wf).time; % Nt x 1
    else
      wfs(wf).time   = wfs(wf).t0_raw + 1/param_extracted.radar.fs*(0:wfs(wf).Nt-1).'; % Nt x 1
      fprintf('deramp: Custom time\n');
    end
    
    %     wfs(wf).signal = [];
    %
    %     wfs(wf).ref_signal = []; % Nt x Ntx
    %     wfs(wf).ref_time = []; % Nt x 1 % Default is wfs(wf).time
    %     wfs(wf).ref_freq_signal = []; % Nt x Ntx
    %     wfs(wf).freq = []; %( [-wfs(wf).Nt+1:wfs(wf).Nt-1]*param.radar.fs/2/wfs(wf).Nt ).';
    
  elseif strcmpi(radar_type,'pulsed')
    wfs(wf).deramp = 0; % Default 0
    
    wfs(wf).time   = load_wfs(wf).time; % Nt x 1
    wfs(wf).signal = tukeywin_cont(wfs(wf).time/wfs(wf).Tpd-0.5,wfs(wf).tukey) * 0.5 .* exp(1i*pi*wfs(wf).chirp_rate*wfs(wf).time.^2);%  * repmat(wfs(wf).tx_weights,wfs(wf).Nt,1); % Nt x Ntx
    
    wfs(wf).ref_signal = conj(flip(wfs(wf).signal,1)); % Nt x Ntx
    wfs(wf).ref_time = wfs(wf).time; % Nt x 1 % Default is wfs(wf).time
    wfs(wf).ref_freq_signal = fft(wfs(wf).ref_signal/norm(wfs(wf).ref_signal),wfs(wf).Nt); % Nt x Ntx
    wfs(wf).freq = [];%( [-wfs(wf).Nt+1:wfs(wf).Nt-1]*param.radar.fs/2/wfs(wf).Nt ).';
    
  end
  
  wfs(wf).Nt_raw = wfs(wf).Nt;
  
  %   wfs(wf).IF_filter_idx = []; % Nrx x Nx
  %   wfs(wf).NCO = []; % Nrx x Nx
  %   wfs(wf).DDC_filter_idx = []; % Nrx x Nx
  %
  %   wfs(wf).dec = [];
  %   wfs(wf).Vpp_scale = param_extracted.radar.Vpp_scale;
  %   wfs(wf).adc_bits = param_extracted.radar.adc_bits;
  %   wfs(wf).type = [];
  
end

param.signal.wfs = wfs;
param.radar.wfs = wfs;

param.sim.Ntx = Ntx;
param.sim.Nrx = Nrx;

% param.radar.wfs(param.sim.wfs).rx_paths = param.sim.rxpath;

% FINISH param.radar before this
%% Use param_extracted to populate the new param

param.season_name         = sprintf('%ssim',param_extracted.season_name);
param.day_seg             = param_extracted.day_seg;
param.radar_name          = param_extracted.radar_name;
param.user_name           = param_extracted.user_name;
param.param_file_version  = param_extracted.param_file_version;
param.sw_version          = param_extracted.sw_version;
param.old_fn              = param_extracted.fn;
% param.load_data           = param_extracted.load_data;

param.cmd       = param_extracted.cmd;
param.records   = param_extracted.records;
param.qlook     = param_extracted.qlook;
param.sar       = param_extracted.sar;
param.array     = param_extracted.array;
param.post      = param_extracted.post;

if isfield(param_extracted,'analysis_noise')
  param.analysis_noise  = param_extracted.analysis_noise;
end
if isfield(param_extracted,'analysis_spec')
  param.analysis_spec  = param_extracted.analysis_spec;
end

% param.radar was done in prev section

exec_good = 1;


end