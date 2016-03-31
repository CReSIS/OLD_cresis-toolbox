function write_cresis_xml(param)
% write_cresis_xml(param)
%
% Creates NI digital system XML files (configuration files).
%
% param
%  .radar_name = 'mcords3', 'mcords4' or 'mcords5'
%  .num_chan = number of receiving channels (used to determine data
%    rate)
%  .aux_dac = DDS Aux DAC settings
%  .TTL_prog_delay = DDS TTL programming delay
%  .fs = sampling frequency (used to determine start/stop record bins)
%  .fs_sync = sync clock used to determine TTL control settings
%  .fs_dds = used to determine start/stop frequency for down-chirp
%  .TTL_mode = 3 length vector [pre_guard post_guard TTL_delay_from_tx]
%    pre_guard: time before transmit starts to go high (normally ~1 us)
%    post_guard: time after transmit ends to stay high (normally ~100 ns)
%    TTL_delay: offset of transmit start from TTL prog delay (a negative
%      value means that the TTL prog delay is after the transmit starts)
%  .version = NI XML file version
%  .max_tx = maximum DDS amplitude (0 to 2^16-1). If given as a scalar,
%    then a constant max is used for all DDS. If it is a vector is should
%    be the same length as param.tx_weights.
%  .max_data_rate = in MB/sec
%  .flight_hours = used to determine data volume
%  .sys_delay = system delay before transmit starts (used to determine
%     start/stop record bins
%  .max_duty_cycle = used to verify settings do not exceed this duty cycle
%    (0 to 1)
%  .create_IQ = create transmit I and Q waveforms if true
%  .tg = time gate fields
%   .staged_recording = waveform timing set to staggered overlapping stages
%     to cover full time gate if this is true. Low-gain, med-gain,
%     high-gain should be given in chronological order in your param
%     settings since that is what is assumed for staged recording. If false,
%     each waveform will cover the full time gate.
%   .altitude_guard = surface timing will be +/- this many meters. May
%     be overridden by waveform specific fields.
%   .Haltitude = nominal meters above ground level. May
%     be overridden by waveform specific fields.
%   .Hice_thick = meters of ice thickness. May
%     be overridden by waveform specific fields.
%  .fn = output filename
%  .prf = pulse repetition frequency (Hz)
%  .presums = row vector of presums for each waveform
%  .wfs = struct array of waveform definitions
%   .atten = receiver attenuator settings (same for all channels)
%   .delay = 1x8 row vector of optimal transmit delays (ns)
%   .f0 = start frequency (Hz)
%   .f1 = stop frequency (Hz)
%   .tx_mask = 1x8 row vector of logicals (1 disables transmit)
%   .tg = waveform specific time gate fields
%    .altitude_guard = surface timing will be +/- this many meters
%    .Haltitude = nominal meters above ground level
%    .Hice_thick = meters of ice thickness
%   .tx_weights = 1x8 row vector of DDS amplitudes (will be scaled to max_tx)
%   .tukey = Tukey weighting to use on transmit pulse 
%   .Tpd = pulse duration in seconds
%   .phase = 1x8 row vector of optimal transmit phase weights (deg)
%   .f0 = start frequency (Hz) if not specified in .wfs field
%   .f1 = stop frequency (Hz) if not specified in .wfs field
%
% param.wfs.delay = [1 0 0 0 0 0 0 0]: causes channel 1 to transmit 1 ns later
%   than all the other channels.
% param.phase = [30 0 0 0 0 0 0 0]: causes channel 1 to lead in phase by 30
%   deg
% param.tx_weights = [2 1 1 1 1 1 1 1]: causes channel 1 to be twice the
%   voltage of the other channels (i.e. 4x the power)
%
% Examples: see create_settings_SEASON_NAME.m
%
% Author: John Paden
%
% See also read_cresis_xml.m

physical_constants;

% Define waveforms
Tpd_base_duration = 1e-6;

xml_version = param.xml_version;
cresis_xml_mapping;

if isfield(param,'TTL_prog_delay') && ~isempty(param.TTL_prog_delay)
  TTL_prog_delay = param.TTL_prog_delay;
else
  error('param.TTL_prog_delay field missing');
end

if isfield(param,'fs_dds') && ~isempty(param.fs_dds)
  fs_dds = param.fs_dds;
else
  error('param.fs_dds field missing');
end

if isfield(param,'fs') && ~isempty(param.fs)
  fs = param.fs;
else
  error('param.fs field missing');
end

if isfield(param,'fs_sync') && ~isempty(param.fs_sync)
  fs_sync = param.fs_sync;
else
  error('param.fs_sync field missing');
end

if isfield(param,'TTL_mode') && ~isempty(param.TTL_mode)
  % Duration before start of pulse, duration after pulse, offset of pulse
  % from TTL_prog_delay * fs/2
  % For example 2014 Ant DC8: [1.3e-6 0.3e-6 -880e-9]
  TTL_mode = param.TTL_mode;
else
  TTL_mode = {};
end

if isfield(param,'Tpd_base_duration') && ~isempty(param.Tpd_base_duration)
  % Manually set the base pulse duration
  Tpd_base_duration = param.Tpd_base_duration;
else
  Tpd_base_duration = 1e-6;
end

settings_enc(1).('Version') = reshape(char(param.version),[1 length(param.version)]);
if isfield(param,'DDC_freq')
  %% DDC Settings
  settings_enc(1).DDCZ20Ctrl.NCOZ20freq = reshape(uint32(param.DDC_freq/1e6),[1 1]);
  settings_enc(1).DDCZ20Ctrl.DDCZ20sel.Val = reshape(uint32(param.DDC_select),[1 1]);
  settings_enc(1).DDCZ20Ctrl.DDCZ20sel.Choice = {'Non DDC','DDC4 En','DDC8 En'};
  settings_enc(1).DDCZ20Ctrl.samplingZ20freq = reshape(double(fs/1e6),[1 1]);
end
settings_enc(1).(config_var_enc)(1).('SYNCZ20CLK') = reshape(double(fs_sync),[1 1]);
settings_enc(1).(config_var_enc)(1).('PRF') = reshape(double(param.prf),[1 1]);
settings_enc(1).(config_var_enc)(1).('Delay') = reshape(uint16(TTL_prog_delay),[1 1]);
settings_enc(1).(config_var_enc)(1).('BaseZ20Len') = reshape(double(Tpd_base_duration),[1 1]);
settings_enc(1).(config_var_enc)(1).('RAMZ20Taper') = reshape(double(param.tukey),[1 1]);
settings_enc(1).(config_var_enc)(1).('AUXZ20DACZ20Z28HEXZ29') = reshape(uint8(param.aux_dac),[1 8]);
if any(param.tx_weights > param.max_tx)
  error('Tx weights too high');
end
settings_enc(1).(config_var_enc)(1).(ram_var_enc) = reshape(uint16(param.tx_weights),[1 8]);

if length(param.wfs)*(1+param.create_IQ) > 16
  error('Too many waveforms');
end
if param.create_IQ
  settings_enc(1).(config_var_enc)(1).(wave_var_enc) = reshape(int32(2*length(param.wfs)),[1 1]);
else
  settings_enc(1).(config_var_enc)(1).(wave_var_enc) = reshape(int32(length(param.wfs)),[1 1]);
end

duty_cycle = 0;
for wf = 1:length(param.wfs)
  if isfield(param,'Tpd')
    param.wfs(wf).Tpd = param.Tpd;
  end
  param.wfs(wf).Tpd = floor(param.wfs(wf).Tpd/Tpd_base_duration) * Tpd_base_duration;
  duty_cycle = duty_cycle + param.presums(wf)/sum(param.presums)*param.wfs(wf).Tpd*param.prf;
end
if duty_cycle > param.max_duty_cycle
  error('Duty cycle too high');
end

%% Convert staged_recording logical setting (0 or 1) into array of
% waveform positions in the staged recording. Staged recording is used to
% record overlapping range gates with different gains to provide a single
% long record with time varying dynamic range. The pulse durations and
% height information are used to do this. Usually the first waveform
% is a low gain waveform for the first part of the return and a second high
% gain waveform is for the second part of the return.  The return can be
% broken into any number of waveforms.  All returns in a given stage are
% assumed to have the same pulse duration and only the first waveform in
% each stage is examined to determine the range gates.
% First stage is assumed to desire capture of ice surface (Haltitude)
% Last stage is assumed to desire capture of ice bottom (Haltitude,
% Hice_thick)
%
% Examples:
%
% An empty array or zero causes every waveform to capture the whole range
% gate. An array with just a one in it, causes the default staged recording
% approach which is waveform N is the Nth stage.
%
% param.tg.staged_recording: [1 2 3]
%   This would cause waveform 1 to be the first stage
%   waveform 2 to be the second stage
%   waveform 3 to be the third stage
%
% param.tg.staged_recording: [3 2 1]
%   This would cause waveform 3 to be the first stage
%   waveform 2 to be the second stage
%   waveform 1 to be the third stage
%
% param.tg.staged_recording: [1 2 3 3 3]
%   This would cause waveform 1 to be the first stage
%   waveform 2 to be the second stage
%   waveform 3, 4, and 5 to be the third stage
if ~isempty(param.tg.staged_recording)
  if param.tg.staged_recording(1) == 0
    param.tg.staged_recording = [];
  elseif length(param.tg.staged_recording) ~= length(param.wfs)
    param.tg.staged_recording = 1:length(param.wfs);
  end
end

%% Create the outputs for each waveform
for wf = 1:length(param.wfs)
  Tpd = param.wfs(wf).Tpd;
  
  if isfield(param,'phase')
    phase_setting = param.phase;
  else
    phase_setting = param.wfs(wf).phase;
  end
  settings_enc(1).(config_var_enc)(1).('Waveforms')(wf).(phase_var_enc) = reshape(double(phase_setting),[1 8]);
  
  %% TTL Settings
  if ~isempty(TTL_mode)
    if ischar(TTL_mode)
      % TTL_mode = 'rx_only'
      % Receive only mode (TTL transmit mode is never asserted)
      TTL_start = 0;
      TTL_duration(1:8) = 0;
    elseif any(strcmpi(param.radar_name,{'mcords3','mcords4','mcords5'}))
      % Duration before start of pulse, duration after pulse, offset of pulse
      % from TTL_prog_delay * fs/2
      % For example 2015 Gr LC130: [2.5e-6 260e-9 -1100e-9]
      TTL_start = round((TTL_prog_delay/fs_sync + TTL_mode(3) - TTL_mode(1))*fs_sync);
      TTL_duration(1:8) = round((TTL_mode(1) + TTL_mode(2) + Tpd) * fs_sync);
    else
      error('Not supported');
    end
    
    settings_enc(1).(config_var_enc)(1).('Waveforms')(wf).(ttl_start_var_enc) = reshape(uint16(TTL_start*ones([1 8],'uint16')),[1 8]);
    settings_enc(1).(config_var_enc)(1).('Waveforms')(wf).(ttl_length_var_enc) = reshape(uint16(TTL_duration),[1 8]);
    
  else
    if strcmpi(param.radar_name,'mcords3')
        % TTL start time (317 for 1e9/9 fs clock and 650 TTL delay)
        TTL_start = TTL_prog_delay + round((317-650) / (1e9/9/2) * fs/2);
        
        % TTL1 length (original durations were for 1e9/9 fs clock)
        x = [1 3 10]*1e-6;
        y = [354 472 845] / (1e9/9/2) * fs/2;
        %plot(x,y);
        p = polyfit(x,y,1);
        TTL_duration = round(polyval(p,Tpd));
        
        % TTL2 length (original durations were for 1e9/9 fs clock)
        x = [1 3 10]*1e-6;
        y = [390 495 888] / (1e9/9/2) * fs/2;
        %plot(x,y);
        p = polyfit(x,y,1);
        TTL_duration(2:8) = round(polyval(p,Tpd));
      
      settings_enc(1).(config_var_enc)(1).('Waveforms')(wf).(ttl_start_var_enc) = reshape(uint16(TTL_start*ones([1 8],'uint16')),[1 8]);
      settings_enc(1).(config_var_enc)(1).('Waveforms')(wf).(ttl_length_var_enc) = reshape(uint16(TTL_duration),[1 8]);
      
    elseif strcmpi(param.radar_name,'mcords4')
      settings_enc(1).(config_var_enc)(1).('Waveforms')(wf).(ttl_start_var_enc) = reshape(uint16(542*ones([1 8],'uint16')),[1 8]);
      
      x = [1 3 10]*1e-6;
      y = [126 251 689];
      %plot(x,y);
      p = polyfit(x,y,1);
      TTL_duration = round(polyval(p,Tpd));
      
      settings_enc(1).(config_var_enc)(1).('Waveforms')(wf).(ttl_length_var_enc) = reshape(uint16(TTL_duration*ones([1 8],'uint16')),[1 8]);
      
    elseif strcmpi(param.radar_name,'mcords5')
      error('Not supported');
    end
    
  end

  atten = param.wfs(wf).atten;
  atten = round(atten*2)/2;
  atten1 = uint8(atten - 31.5);
  atten2 = uint8(atten - double(atten1));
  
  if length(atten1) == 1
    if strcmpi(param.radar_name,'mcords3')
%       settings_enc(1).Atten.('AttenZ20Z30') = reshape(uint8(atten1*ones([1 8],'uint8')),[1 8]);
%       settings_enc(1).Atten.('AttenZ31') = reshape(uint8(atten2*ones([1 8],'uint8')),[1 8]);
    end
    settings_enc(1).(config_var_enc)(1).('Waveforms')(wf).('AttenuatorZ20Z31') = reshape(uint8(atten1*ones([1 8],'uint8')),[1 8]);
    settings_enc(1).(config_var_enc)(1).('Waveforms')(wf).('AttenuatorZ20Z32') = reshape(uint8(atten2*ones([1 8],'uint8')),[1 8]);
  else
    if strcmpi(param.radar_name,'mcords3')
%       settings_enc(1).Atten.('AttenZ20Z30') = reshape(uint8(atten1),[1 8]);
%       settings_enc(1).Atten.('AttenZ31') = reshape(uint8(atten2),[1 8]);
    end
    settings_enc(1).(config_var_enc)(1).('Waveforms')(wf).('AttenuatorZ20Z31') = reshape(uint8(atten1),[1 8]);
    settings_enc(1).(config_var_enc)(1).('Waveforms')(wf).('AttenuatorZ20Z32') = reshape(uint8(atten2),[1 8]);
  end
  
  settings_enc(1).(config_var_enc)(1).('Waveforms')(wf).('LenZ20Mult') ...
    = reshape(uint16(round(Tpd/Tpd_base_duration)),[1 1]);
  settings_enc(1).(config_var_enc)(1).('Waveforms')(wf).('Presums') ...
    = reshape(uint16(param.presums(wf)),[1 1]);
  
  if isfield(param,'tx_mask')
    tx_mask = bin2dec(char('0' + param.tx_mask));
  else
    tx_mask = bin2dec(char('0' + param.wfs(wf).tx_mask));
  end
  settings_enc(1).(config_var_enc)(1).('Waveforms')(wf).('TXZ20Mask') = reshape(uint8(tx_mask),[1 1]);
  
  if isfield(param,'sys_delay')
    Tsystem_delay = param.sys_delay;
  else
    error('You need to set param.sys_delay. Old value was 10.6e-6.');
    %Tsystem_delay = 10.6e-6; % 10.6 us to program DDS
  end
  
  Tguard = [];
  Haltitude = [];
  Hice_thick = [];
  if isfield(param,'wfs') && length(param.wfs) >= wf && isfield(param.wfs(wf),'tg')
    if isfield(param.wfs(wf).tg,'altitude_guard') && ~isempty(param.wfs(wf).tg.altitude_guard)
      Tguard = param.wfs(wf).tg.altitude_guard / (3e8/2);
    else
      Tguard = param.tg.altitude_guard / (3e8/2);
    end
    if isfield(param.wfs(wf).tg,'Haltitude') && ~isempty(param.wfs(wf).tg.Haltitude)
      Haltitude = param.wfs(wf).tg.Haltitude;
    else
      Haltitude = param.tg.Haltitude;
    end
    if isfield(param.wfs(wf).tg,'Hice_thick') && ~isempty(param.wfs(wf).tg.Hice_thick)
      Hice_thick = param.wfs(wf).tg.Hice_thick;
    else
      Hice_thick = param.tg.Hice_thick;
    end
  end
  if isempty(Tguard)
    Tguard = param.tg.altitude_guard / (3e8/2);
  end
  if isempty(Haltitude)
    Haltitude = param.tg.Haltitude;
  end
  if isempty(Hice_thick)
    Hice_thick = param.tg.Hice_thick;
  end
  
  if isfield(param.tg,'look_angle')
    look_angle = param.tg.look_angle(wf);
  else
    look_angle = 0;
  end
  if ~isempty(param.tg.staged_recording)
    current_stage = param.tg.staged_recording(wf);
    next_stage_wf = find(param.tg.staged_recording == current_stage+1,1);
    if current_stage == 1
      % First Stage
      Tstart = Haltitude / (3e8/2);
      Tend = Haltitude / (3e8/2) + param.wfs(next_stage_wf).Tpd;
    elseif isempty(next_stage_wf)
      % Last Stage
      Tstart = Haltitude / (3e8/2) + param.wfs(wf).Tpd;
      Tend = Haltitude / (3e8/2)  + Hice_thick / (3e8/2/sqrt(3.15));
    else
      % In between first and last stage
      Tstart = Haltitude / (3e8/2) + param.wfs(wf).Tpd;
      Tend = Haltitude / (3e8/2) + param.wfs(next_stage_wf).Tpd;
    end
  else
    look_angle_ice = asind(sind(look_angle)/sqrt(er_ice));
    Tstart = Haltitude / (3e8/2);
    Tend = Haltitude/cosd(look_angle) / (3e8/2)  + Hice_thick/cosd(look_angle_ice) / (3e8/2/sqrt(3.15));
  end
  if strcmpi(param.radar_name,'mcords3')
    bin_start = round((Tstart + Tsystem_delay - Tguard)  * fs);
    bin_stop = round((Tend + param.wfs(wf).Tpd + Tguard + Tsystem_delay)  * fs);
  elseif any(strcmpi(param.radar_name,{'mcords4','mcords5'}))
    bin_start = round((Tstart + Tsystem_delay - Tguard)  * fs/8);
    bin_stop = round((Tend + param.wfs(wf).Tpd + Tguard + Tsystem_delay)  * fs/8);
  else
    error('Unsupported radar')
  end
  
  settings_enc(1).(config_var_enc)(1).('Waveforms')(wf).('RecordZ20Stop') = reshape(uint16(bin_stop),[1 1]);
  if uint16(bin_start) < 32
      % Force record start to be at least 32
      settings_enc(1).(config_var_enc)(1).('Waveforms')(wf).('RecordZ20Start') = reshape(uint16(32),[1 1]);
  else
      settings_enc(1).(config_var_enc)(1).('Waveforms')(wf).('RecordZ20Start') = reshape(uint16(bin_start),[1 1]);
  end
  
  if isfield(param,'f0')
    f0 = param.f0;
  else
    f0 = param.wfs(wf).f0;
  end
  if isfield(param,'f1')
    f1 = param.f1;
  else
    f1 = param.wfs(wf).f1;
  end
  if f0 > f1
    f0 = fs_dds - f0;
    f1 = fs_dds - f1;
  end
  settings_enc(1).(config_var_enc)(1).('Waveforms')(wf).('StartZ20Freq') ...
    = reshape(double(f0*ones([1  8])),[1 8]);
  settings_enc(1).(config_var_enc)(1).('Waveforms')(wf).('StopZ20Freq') ...
    = reshape(double(f1*ones([1  8])),[1 8]);
  
  if isfield(param,'delay')
    delay = param.delay;
  else
    delay = param.wfs(wf).delay;
  end
  settings_enc(1).(config_var_enc)(1).('Waveforms')(wf).('Delay') ...
    = reshape(double(delay),[1 8]);
end

if param.create_IQ
  settings_enc(1).(config_var_enc)(1).('Waveforms') ...
    = settings_enc(1).(config_var_enc)(1).('Waveforms')(reshape([1:length(param.wfs); 1:length(param.wfs)],[1 2*length(param.wfs)]));
  for wf = 2:2:length(settings_enc(1).(config_var_enc)(1).('Waveforms'))
    settings_enc(1).(config_var_enc)(1).('Waveforms')(wf).(phase_var_enc) ...
      = angle(exp(j*settings_enc(1).(config_var_enc)(1).('Waveforms')(wf).(phase_var_enc)/180*pi + j*pi/2))*180/pi;
  end
end

if strcmpi(param.radar_name,'mcords4')
  settings_enc(1).('FPGAZ20Configuration')(1).('HIZ20Presums') = reshape(uint8([50 ]),[1 1]);
  settings_enc(1).('FPGAZ20Configuration')(1).('HIZ20UpdateZ23') = reshape(uint8([0 ]),[1 1]);
  settings_enc(1).('FPGAZ20Configuration')(1).('HIZ20Z20WaveformZ23') = reshape(uint8([0 ]),[1 1]);
  settings_enc(1).('FPGAZ20Configuration')(1).('HIZ20cntZ20rst') = reshape(logical([0 ]),[1 1]);
end

out_xml_fn = param.fn;

%% Write output
fprintf('Writing %s\n', out_xml_fn);

% Check data rate
num_sam = 0;
presums = 0;
for wf = 1:length(settings_enc.(config_var_enc).Waveforms)
  num_sam = num_sam + double(settings_enc.(config_var_enc).Waveforms(wf).RecordZ20Stop - settings_enc.(config_var_enc).Waveforms(wf).RecordZ20Start);
  presums = presums + double(settings_enc.(config_var_enc).Waveforms(wf).Presums);
end
EPRF = double(settings_enc.(config_var_enc).PRF) / presums;

sample_size = 2;

if strcmpi(param.radar_name,'mcords3')
  data_rate = EPRF * num_sam * sample_size * param.num_chan;
elseif any(strcmpi(param.radar_name,{'mcords4','mcords5'}))
  data_rate = EPRF * num_sam * 8 * sample_size * param.num_chan / 2^param.DDC_select;
else
  error('Unsupported radar')
end
fprintf('  EPRF: %.1f Hz\n', EPRF);
lambda_fc = 3e8/ (0.5*(settings_enc.(config_var_enc).Waveforms(1).StopZ20Freq(1) + settings_enc.(config_var_enc).Waveforms(1).StartZ20Freq(1)));
fprintf('  fc Nyquist sampling rate velocity: %.1f m/s (%.1f knots)\n', lambda_fc/4 * EPRF, lambda_fc/4 * EPRF / 0.5515)
fprintf('  Data rate: %.1f MB/sec (%.1f hours: %.1f TB)\n', data_rate / 2^20, param.flight_hours, param.flight_hours*3600*data_rate/2^40);

Tguard = 1e-6;
if 1/settings_enc.(config_var_enc).PRF - Tguard - double(settings_enc.(config_var_enc).Waveforms(wf).RecordZ20Stop)/fs <= 0
  error('Data recording window, %.2f us, is too long for this PRI, %.2f us, including %.2f us time guard', ...
    double(settings_enc.(config_var_enc).Waveforms(wf).RecordZ20Stop)/fs*1e6, 1/settings_enc.(config_var_enc).PRF*1e6, Tguard*1e6);
end

if data_rate / 2^20 > param.max_data_rate
  error('Data rate %f MB/sec is too high', data_rate/2^20);
end

if xml_version >= 2.0
  settings_enc.sys.DDCZ20Ctrl = settings_enc.DDCZ20Ctrl;
  settings_enc.sys.DDSZ5FSetup = settings_enc.DDSZ5FSetup;
  settings_enc.sys.XMLZ20FileZ20Path = {struct('type','Path','values',[])};
  settings_enc.sys.XMLZ20FileZ20Path{1}.values = {out_xml_fn};
  settings_enc.sys.xmlversion = {struct('type','String','values',[])};
  settings_enc.sys.xmlversion{1}.values = {'2.0'};
  settings_enc = rmfield(settings_enc,'DDCZ20Ctrl');
  settings_enc = rmfield(settings_enc,'DDSZ5FSetup');
end

fid = fopen(out_xml_fn,'w');
fprintf(fid,'<?xml version=''1.0'' standalone=''yes'' ?>\n');
fprintf(fid,'<LVData xmlns="http://www.ni.com/LVData">\n');
write_ni_xml_object(settings_enc,fid,true,struct('array_list','Waveforms','enum_list','DDCZ20sel'));
fprintf(fid,'</LVData>');
fclose(fid);

return;

