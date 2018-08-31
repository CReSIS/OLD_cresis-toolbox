function configs = write_radar_xml(param)
% configs = write_radar_xml(param)
%
% Creates Arena and NI digital system XML files (configuration files).
%
% param
%  .num_chan = number of receiving channels (used to determine data
%    rate)
%  .TTL_prog_delay = DDS TTL programming delay
%  .fs: ADC sampling frequency (used to determine start/stop record bins)
%  .fs_TTL: TTL clock used to determine TTL control settings
%  .fs_dac: DAC sampling frequency
%  .TTL_mode = 3 length vector [pre_guard post_guard TTL_delay_from_tx]
%    pre_guard: time before transmit starts to go high (normally ~3 us)
%    post_guard: time after transmit ends to stay high (normally ~350 ns)
%    TTL_delay: offset of transmit start from TTL prog delay (a negative
%      value means that the TTL prog delay is after the transmit starts)
%  .max_data_rate = in MB/sec
%  .flight_hours = used to determine data volume
%  .sys_delay = system delay before transmit starts (used to determine
%     start/stop record bins
%  .max_duty_cycle = used to verify settings do not exceed this duty cycle
%    (0 to 1)
%
%  .version = NI XML file version
%  .max_DDS_amp = maximum DDS amplitude (0 to 2^16-1). If given as a scalar,
%    then a constant max is used for all DDS. If it is a vector is should
%    be the same length as param.tx_weights. Alternate field name is max_tx.
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
%   .rg_start_offset = adjust the range gate start by this amount in meters.
%     Positive values make the start later (i.e. shorter range gate).
%   .rg_stop_offset: adjust the range gate stop by this amount in meters.
%     Positive values make the stop later (i.e. lengthen the range gate).
%   .start_ref: cell array containing 'surface' or 'bottom' for each
%     waveform. The start time is relative to the specified reference. The
%     default is the surface.
%   .stop_ref: cell array containing 'surface' or 'bottom' for each
%     waveform. The stop time is relative to the specified reference. The
%     default is the bottom.
%  .fn = output filename
%  .prf = pulse repetition frequency (Hz)
%  .presums = row vector of presums for each waveform
%  .wfs = struct array of waveform definitions
%   .atten = receiver attenuator settings (same for all channels)
%   .delay = 1 by N_tx row vector of optimal transmit delays (ns)
%   .f0 = start frequency (Hz)
%   .f1 = stop frequency (Hz)
%   .tx_enable = 1 by N_tx row vector of logicals (1 disables transmit)
%   .tg = waveform specific time gate fields
%    .altitude_guard = surface timing will be +/- this many meters
%    .Haltitude = nominal meters above ground level
%    .Hice_thick = meters of ice thickness
%   .tx_weights = 1x8 row vector of DDS amplitudes (will be scaled to max_tx)
%   .tukey = Tukey weighting to use on transmit pulse
%   .Tpd = pulse duration in seconds
%   .phase = 1 by N_tx row vector of optimal transmit phase weights (deg)
%   .f0 = start frequency (Hz) if not specified in .wfs field
%   .f1 = stop frequency (Hz) if not specified in .wfs field
%
% param.wfs(wf).delay = [1 0 0 0 0 0 0 0]: causes channel 1 to transmit 1
%   ns later than all the other channels for waveform wf.
%
% param.phase = [30 0 0 0 0 0 0 0]: causes channel 1 to lead in phase by 30
%   deg for all waveforms
% param.tx_weights = [2 1 1 1 1 1 1 1]: causes channel 1 to be twice the
%   voltage of the other channels (i.e. 4x the power) for all waveforms
%
% Examples: see create_settings_SEASON_NAME_RADAR_NAME.m
%
% Author: John Paden
%
% See also: read_arena_xml, read_cresis_xml, write_arena_xml, write_ni_xml,
% write_radar_xml

physical_constants;

for wf = 1:length(param.wfs)
  if ~isfield(param.wfs(wf),'f0') || isempty(param.wfs(wf).f0)
    param.wfs(wf).f0 = param.f0;
  end
  if ~isfield(param.wfs(wf),'f1') || isempty(param.wfs(wf).f1)
    param.wfs(wf).f1 = param.f1;
  end
  if ~isfield(param.wfs(wf),'DDC_freq') || isempty(param.wfs(wf).DDC_freq)
    param.wfs(wf).DDC_freq = param.DDC_freq;
  end
  if ~isfield(param.wfs(wf),'tukey') || isempty(param.wfs(wf).tukey)
    param.wfs(wf).tukey = param.tukey;
  end
  if ~isfield(param.wfs(wf),'Tpd') || isempty(param.wfs(wf).Tpd)
    param.wfs(wf).Tpd = param.Tpd;
  end
  if ~isfield(param.wfs(wf),'zeropimods') || isempty(param.wfs(wf).zeropimods)
    param.wfs(wf).zeropimods = param.zeropimods;
  end
  if ~isfield(param.wfs,'name')
    param.wfs(wf).name = '';
  end
  if ~isfield(param.wfs(wf),'phase') || isempty(param.wfs(wf).phase)
    param.wfs(wf).phase = param.phase;
  end
  if ~isfield(param.wfs(wf),'delay') || isempty(param.wfs(wf).delay)
    param.wfs(wf).delay = param.delay;
  end
  if ~isfield(param.wfs(wf),'tx_enable') || isempty(param.wfs(wf).tx_enable)
    param.wfs(wf).tx_enable = param.tx_enable;
  end
  if ~isfield(param.wfs(wf),'tx_weights') || isempty(param.wfs(wf).tx_weights)
    param.wfs(wf).tx_weights = param.tx_weights;
  end
  if any(param.wfs(wf).tx_weights > param.max_tx)
    error('Tx weights param.wfs(%d).tx_weights too high for waveform %d', wf, wf);
  end
end

duty_cycle = 0;
for wf = 1:length(param.wfs)
  duty_cycle = duty_cycle + param.wfs(wf).presums/sum(cell2mat({param.wfs.presums}))*param.wfs(wf).Tpd*param.prf;
end
if duty_cycle > param.max_duty_cycle
  error('Duty cycle too high');
end

if isfield(param,'prf_multiple') && ~isempty(param.prf_multiple)
  if any(mod(param.prf_multiple/param.prf,1))
    error('param.prf (%g) must be a factor of default.prf_multiple (%s) for coherent noise cancelling to work.', param.prf, mat2str_generic(param.prf_multiple));
  end
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
  if all(param.tg.staged_recording == 0)
    param.tg.staged_recording = zeros(size(param.wfs));
  elseif length(param.tg.staged_recording) ~= length(param.wfs)
    param.tg.staged_recording = 1:length(param.wfs);
  end
end

%% Create the outputs for each waveform
for wf = 1:length(param.wfs)
  Tpd = param.wfs(wf).Tpd;

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
  
  if isfield(param.tg,'look_angle_deg')
    look_angle_deg = param.tg.look_angle_deg(wf);
  else
    look_angle_deg = 0;
  end
  
  if isfield(param.tg,'start_ref')
    start_ref = param.tg.start_ref{wf};
  else
    start_ref = 'surface';
  end
  if isfield(param.tg,'stop_ref')
    stop_ref = param.tg.stop_ref{wf};
  else
    stop_ref = '';
  end
  
  if isfield(param.tg,'rg_start_offset')
    Trg_start_offset = param.tg.rg_start_offset(wf) / (3e8/2);
  else
    Trg_start_offset = 0;
  end
  if isfield(param.tg,'rg_stop_offset')
    Trg_stop_offset = param.tg.rg_stop_offset(wf) / (3e8/2);
  else
    Trg_stop_offset = 0;
  end
  if isfield(param.tg,'Hice_thick_min')
    Hice_thick_min = param.tg.Hice_thick_min;
  else
    Hice_thick_min = 0;
  end
  
  look_angle_ice_deg = asind(sind(look_angle_deg)/sqrt(er_ice));
  current_stage = param.tg.staged_recording(wf);
  next_stage_wf = find(param.tg.staged_recording == current_stage+1,1);
  if strcmpi(start_ref,'surface')
    Tstart_ref = Haltitude / (3e8/2);
  elseif strcmpi(start_ref,'bottom')
    Tstart_ref = Haltitude / (3e8/2) + Hice_thick_min / (3e8/2/sqrt(er_ice));
  else
    error('Invalid start ref %s', start_ref);
  end
  if strcmpi(stop_ref,'surface')
    Tstop_ref = Haltitude/cosd(look_angle_deg) / (3e8/2);
  elseif strcmpi(stop_ref,'bottom')
    Tstop_ref = Haltitude/cosd(look_angle_deg) / (3e8/2) + Hice_thick/cosd(look_angle_ice_deg) / (3e8/2/sqrt(er_ice));
  elseif isempty(stop_ref)
    if current_stage > 0 && ~isempty(next_stage_wf)
      % First and in between stages default references surface
      Tstop_ref = Haltitude/cosd(look_angle_deg) / (3e8/2);
    else
      % Last stage default reference is bottom
      Tstop_ref = Haltitude/cosd(look_angle_deg) / (3e8/2) + Hice_thick/cosd(look_angle_ice_deg) / (3e8/2/sqrt(er_ice));
    end
  else
    error('Invalid start ref %s', start_ref);
  end
  
  if current_stage == 0
    % Record the whole range gate
    Tstart = Tstart_ref;
    Tend = Tstop_ref + param.wfs(wf).Tpd;
  elseif current_stage == 1
    % First Stage
    Tstart = Tstart_ref;
    if ~isempty(next_stage_wf)
      % Cover the pulse duration of the next waveform after the surface
      % return + any off nadir scattering that is to be caught.
      Tend = max(Tstart_ref + param.wfs(next_stage_wf).Tpd + param.wfs(wf).Tpd, Tstop_ref);
    else
      % The first stage is also the last stage
      Tend = Tstop_ref + param.wfs(wf).Tpd;
    end
  elseif isempty(next_stage_wf)
    % Last Stage
    Tstart = Tstart_ref + param.wfs(wf).Tpd;
    Tend = Tstop_ref + param.wfs(wf).Tpd;
  else
    % In between first and last stage
    Tstart = Tstart_ref + param.wfs(wf).Tpd;
    Tend = max(Tstart_ref + param.wfs(next_stage_wf).Tpd + param.wfs(wf).Tpd, Tstop_ref);
  end
  
  param.wfs(wf).Tstart = Tstart-Tguard;
  param.wfs(wf).Tend = Tend+Tguard;
end

if param.create_IQ
  num_wfs = numel(param.wfs);
  for wf = 1:num_wfs
    param.wfs(num_wfs+wf) = param.wfs(wf);
    param.wfs(num_wfs+wf).phase = param.wfs(num_wfs+wf).phase + 90;
  end
end

%% Check data rate (move to write_arena_xml)
% =========================================================================
param.num_chan = 2;
param.sample_size = 4;
param.decimation = 2;
PRI_guard = 1e-6;

num_sam = 0;
presums = 0;
for wf = 1:numel(param.wfs)
  num_sam = num_sam + (param.wfs(wf).Tend-param.wfs(wf).Tstart)*param.fs;
  presums = presums + param.wfs(wf).presums;
end
EPRF = double(param.prf) / presums;

data_rate = EPRF * num_sam * param.sample_size * param.num_chan / param.decimation;

fprintf('  EPRF: %.1f Hz\n', EPRF);
fc = (param.wfs(1).f0+param.wfs(wf).f1)/2;
lambda_fc = 3e8/fc;
fprintf('  fc Nyquist sampling rate velocity: %.1f m/s (%.1f knots)\n', lambda_fc/4 * EPRF, lambda_fc/4 * EPRF / 0.5515)
fprintf('  Data rate: %.1f MB/sec (%.1f hours: %.1f TB)\n', data_rate / 2^20, param.flight_hours, param.flight_hours*3600*data_rate/2^40);

for wf = 1:numel(param.wfs)
  if 1/param.prf - PRI_guard - param.wfs(wf).Tend <= 0
    error('Data recording window, %.2f us, is too long for this PRI, %.2f us, including %.2f us time guard', ...
      param.wfs(wf).Tend*1e6, 1/param.prf*1e6, PRI_guard*1e6);
  end
end

if data_rate / 2^20 > param.max_data_rate
  error('Data rate %f MB/sec is too high', data_rate/2^20);
end

%% Write output
if isfield(param,'ni') && ~isempty(param.ni)
  fprintf('Writing %s\n', param.fn);
  out_xml_fn_dir = fileparts(out_xml_fn);
  if ~exist(out_xml_fn_dir,'dir')
    mkdir(out_xml_fn_dir);
  end
  fid = fopen(out_xml_fn,'w');
  fprintf(fid,'<?xml version=''1.0'' standalone=''yes'' ?>\n');
  fprintf(fid,'<LVData xmlns="http://www.ni.com/LVData">\n');
  write_ni_xml_object(settings_enc,fid,true,struct('array_list','Waveforms','enum_list','DDCZ20sel'));
  fprintf(fid,'</LVData>');
  fclose(fid);
end

%% Write RSS Arena XML config file
if isfield(param,'arena') && ~isempty(param.arena)
  % Create XML document
  doc = write_arena_xml([],param);

  out_str = xmlwrite(doc);
  out_str = ['<!DOCTYPE systemXML>' out_str(find(out_str==10,1):end)];
  arena_fn_dir = fileparts(param.arena.fn);
  if ~exist(arena_fn_dir,'dir')
    mkdir(arena_fn_dir);
  end
  fprintf('  Writing Arena XML: %s\n', param.arena.fn);
  fid = fopen(param.arena.fn,'w');
  fwrite(fid,out_str,'char');
  fclose(fid);
end

return;

