function out = fmcw3_NCO(in,param)
% out = fmcw3_NCO(in,param)
% 
% Function that converts between NCO counts and IF frequency or range
%
% param = struct controlling conversion
%  .NCO_step_size (default 20)
%  .fs (default 125e6)
%  .chirp_rate (default 6e9/240e-6)
%  .units (default 'ft'), options are 'm', 'Hz', 'NCO_file', 'NCO'
%
% Examples:
% NCO = fmcw3_NCO(1500)
% NCO = fmcw3_NCO(500,struct('units','m'))
% NCO = fmcw3_NCO([0 1 2 3]*125e6/2,struct('units','Hz'))
% f_IF = fmcw3_NCO(17900,struct('units','NCO_file'))
% f_IF = fmcw3_NCO(895,struct('units','NCO'))

if ~exist('param','var')
  param.NCO_step_size = 20;
  param.fs = 125e6;
  param.chirp_rate = 6e9/240e-6;
  param.units = 'ft';
end

if ~isfield(param,'NCO_step_size') || isempty(param.NCO_step_size)
  param.NCO_step_size = 20;
end
if ~isfield(param,'fs') || isempty(param.fs)
  param.fs = 125e6;
end
if ~isfield(param,'chirp_rate')
  param.chirp_rate = 6e9/240e-6;
end
if ~isfield(param,'units')
  param.units = 'ft';
end

c = 2.997924580003452e+08;
if strcmpi(param.units,'ft')
  td = in*12*2.54/100 / (c/2);
  f_IF = td * param.chirp_rate;
  out = round(((f_IF+param.fs/2) / param.fs / 2) * 2^15/param.NCO_step_size);
elseif strcmpi(param.units,'m')
  td = in / (c/2);
  f_IF = td * param.chirp_rate;
  out = round(((f_IF+param.fs/2) / param.fs / 2) * 2^15/param.NCO_step_size);
elseif strcmpi(param.units,'Hz')
  out = round(((in+param.fs/2) / param.fs / 2) * 2^15/param.NCO_step_size);
elseif strcmpi(param.units,'NCO')
  out = in*param.NCO_step_size*2*param.fs/2^15 - param.fs/2;
elseif strcmpi(param.units,'NCO_file')
  out = in*2*param.fs/2^15 - param.fs/2;
end

return

NCO = fmcw3_NCO(1500)
NCO = fmcw3_NCO(500,struct('units','m'))
NCO = fmcw3_NCO(0e6,struct('units','Hz'))
NCO = fmcw3_NCO(62.5e6,struct('units','Hz'))
NCO = fmcw3_NCO(125e6,struct('units','Hz'))
NCO = fmcw3_NCO(187.5e6,struct('units','Hz'))
f_IF = fmcw3_NCO(17900,struct('units','NCO_file'))
f_IF = fmcw3_NCO(895,struct('units','NCO'))

% ========================================================
% Notes for tracking aliased peaks in the DDC
% ========================================================

[max_val max_freq] = max(lp(g_data(1:4300,:)));
signal_freqs = freq_raw(max_freq);
%signal_freqs = [0 92.6]*1e6;
NCO_freq = 90.5e6;
fs = 125e6;
image_freqs = fs - signal_freqs;
DDC_signals = image_freqs - NCO_freq;
DDC_filter = fs/4;
DDC_alias = floor((DDC_signals--DDC_filter)/(2*DDC_filter));
DDC_signals_aliased = DDC_signals - DDC_alias * 2*DDC_filter;
aliased_freq = NCO_freq+DDC_signals_aliased;

figure(1); clf;
imagesc([],freq_raw/1e6,lp(g_data(1:num_sam(1),:)))
hold on;
plot(signal_freqs/1e6,'b')
plot(aliased_freq/1e6,'k')
hold off;
xlabel('Range lines');
ylabel('IF freq (MHz)');


% ========================================================
% Notes for tracker
% ========================================================
fs = 125e6;
chirp_rate = 6e9/240e-6;

NCO_step_size = 20; % GUI setting

NCO_freq_change_threshold = 20; % NEW GUI setting
power_threshold = 20; % NEW GUI setting
NCO_tracking = true; % NEW GUI setting

% Determine which range bins to search in (account for the fact that
% the digital hardware is oversampled by 2x)
oversample = 2;
Nt = hdr.stop_idx - hdr.start_idx;
search_bins = 1 + round(oversample*Nt*[0.51 0.74]); % NEW GUI setting

% Convert data to power
data = abs(data).^2;
% Incoherently average (multilook) data
data = fir_dec(data,5);

% Search for max
[max_val max_idx] = data(search_bins,rline);

% Convert max_idx to IF frequency (it is important that max_idx is 
% referenced to raw data DC)
new_NCO_freq = max_idx / Nt * fs;
% Convert the new IF frequency into DDC NCO steps
new_NCO_freq = round(((new_NCO_freq+param.fs/2) / param.fs / 2) * (2^15/2)/NCO_step_size);

if NCO_tracking == true ...
    && abs(new_NCO_freq - NCO_freq) > NCO_freq_change_threshold ...
    && max_val > power_threshold
  NCO_freq = new_NCO_freq;
end

% Determine which range bins the DDC is selecting

% NCO_freq(rline) = NCO_freq_step / 2^15 * param.wfs(wf).fs_raw * 2 - 62.5e6;
NCO_freq*NCO_step_size / 2^15 * param.fs/2

new_NCO_freq = max_idx / Nt * fs;
% Convert the new IF frequency into DDC NCO steps
new_NCO_freq = round(((new_NCO_freq+param.fs/2) / param.fs / 2) * 2^15/NCO_step_size);

