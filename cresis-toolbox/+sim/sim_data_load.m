%OM
try hm; end;
[c, WGS84] = physical_constants('c', 'WGS84');

% param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/snow/2012_Greenland_P3sim/20120330/param.mat';
% param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/snow/2013_Greenland_P3sim/20130327/param.mat';
% param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/snow/2018_Antarctica_DC8sim/20181010/param.mat';
param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/rds/2014_Greenland_P3sim/20140325/param.mat';

load(param_fn);

%% FullSim load data

param.load_data.pulse_comp            = 1;
param.load_data.raw_data              = 0;
[hdr,data] = load_data(param);

%% FullSim hdr checks

img         = 1; % support only one for now
wf_adc_idx  = 1;
wf          = param.load_data.imgs{img}(wf_adc_idx,1);
adc         = param.load_data.imgs{img}(wf_adc_idx,2);

% Expected values for traj, range, twtt, time axis
traj = hdr.records{img};
[traj.x, traj.y, traj.z] = geodeticD2ecef(traj.lat,traj.lon,traj.elev, WGS84.ellipsoid);

target  = hdr.param_load_data.target;
range   = sqrt( (traj.x - target.x).^2 + (traj.y - target.y).^2 + (traj.z - target.z).^2 );
twtt    = 2*range/c;
time    = hdr.time{img};
freq    = hdr.freq{img};
freqs   = fftshift(freq);

test_range  = param.sim.range - range;
test_twtt   = param.sim.twtt - twtt;
test_time   = param.radar.wfs(wf).time - hdr.time{img};

if any(test_range); warning('!!! Range mismatch\n'); end
if any(test_twtt); warning('!!! twtt mismatch\n'); end
if any(test_time); warning('!!! time mismatch\n'); end

for compressing_this = 1 % continue; % FullSim hdr checks
  
  fig_title = sprintf('%s_%s',mfilename, 'FullSim hdr checks');
  fig_h = figure('Name',fig_title);
  subplot(131);
  plot(param.traj.lon, param.traj.lat, 'x'); hold on; grid on;
  plot(traj.lon, traj.lat, 'o');
  plot(param.gps.lon, param.gps.lat, 's');
  xlabel('Longitude, Degrees');
  ylabel('Latitude, Degrees');
  legend('Simulated', 'Loaded', 'Reference');
  title('Trajectories');
  subplot(132);
  plot(param.sim.range, 'x'); hold on; grid on;
  plot(range, 'o');
  xlabel('Along-track position, rlines');
  ylabel('Range, meters');
  legend('Simulated', 'Calculated');
  title('Range to target');
  subplot(133);
  plot(param.sim.twtt./1e-6, 'x'); hold on; grid on;
  plot(twtt./1e-6, 'o');
  xlabel('Along-track position, rlines');
  ylabel('TWTT, us');
  legend('Simulated', 'Calculated');
  title('TWTT to target');
  try
    sgtitle(fig_title,'FontWeight','Bold','FontSize',14,'Interpreter','None');
  end
  set(findobj(fig_h,'type','axes'),'FontWeight', 'Bold', 'FontSize',14);
  set(fig_h, 'Position', get(0, 'Screensize'));
  print(gcf, '-dpng', fig_title, '-r300');
  
end % for compressing_this % FullSim hdr checks

%% FullSim data checks

data    = data{img};
tmp     = 20*log10(abs(data));
fdata   = fftshift(fft(data));
ftmp    = 20*log10(abs(fdata));

[data_max, data_max_idx]    = max(data,[],1);
[fdata_max, fdata_max_idx]  = max(fdata,[],1);
[tmp_max, tmp_max_idx]      = max(tmp,[],1);
[ftmp_max, ftmp_max_idx]    = max(ftmp,[],1);
[twtt_min, twtt_min_idx]    = min(twtt);

[Nt,Nx] = size(data);
y       = 1:Nx;

% To compare data_max with expected phase
wave_number = 4*pi/ (c/hdr.param_load_data.radar.wfs(wf).fc);
% expected_phase = exp(1i * wave_number * (-range + range(twtt_min_idx)) );
% expected phase is data dependent
expected_phase = exp(1i * wave_number * (-range + range(twtt_min_idx)) ) * data_max(twtt_min_idx);

for compressing_this = 1 % continue; % FullSim Time Domain
  
  call_sign = 'FullSim Time Domain';
  fig_title = sprintf('%s_%s',mfilename, call_sign);
  fig_h = figure('Name',fig_title);
  subplot(131);
  imagesc(y, time/1e-6, tmp-max(tmp(:)) ,[-30,0] );
  cb = colorbar; cb.Label.String = 'Relative Power, dB';
  grid on; hold on; axis tight;
  plot(y, twtt/1e-6, '--', 'Color', 'g');
  plot(twtt_min_idx, twtt(twtt_min_idx)/1e-6, 's','LineWidth',2);
  xlabel('Along-track position, rlines'); ylabel('Fast-time, us');
  title('PhaseHistory Magnitude');
  subplot(132);
  imagesc(y, time/1e-6, angle(data) );
  cb = colorbar; cb.Label.String = 'Radians';
  grid on; hold on; axis tight;
  plot(y, twtt/1e-6, '--', 'Color', 'g');
  plot(twtt_min_idx, twtt(twtt_min_idx)/1e-6, 's','LineWidth',2);
  xlabel('Along-track position, rlines'); ylabel('Fast-time, us');
  title('PhaseHistory Phase');
  subplot(133);
  imagesc(y, time/1e-6, unwrap(angle(data)) );
  cb = colorbar; cb.Label.String = 'Radians';
  grid on; hold on; axis tight;
  plot(y, twtt/1e-6, '--', 'Color', 'g');
  plot(twtt_min_idx, twtt(twtt_min_idx)/1e-6, 's','LineWidth',2);
  xlabel('Along-track position, rlines'); ylabel('Fast-time, us');
  title('PhaseHistory Phase Unwrap');
  if 0
    plot(time/1e-6, angle(data(:,twtt_min_idx)) ); hold on;
    plot(twtt_min/1e-6, angle(data(data_max_idx(twtt_min_idx), twtt_min_idx)) , 's','LineWidth',2);
    xlabel('Fast-time, us');
    ylabel('Phase, Radians');
    legend('rline','expected');
    title('At closest range');
  end
  try
    sgtitle(fig_title,'FontWeight','Bold','FontSize',14,'Interpreter','None');
  end
  set(findobj(fig_h,'type','axes'),'FontWeight', 'Bold', 'FontSize',14);
  set(fig_h, 'Position', get(0, 'Screensize'));
  print(gcf, '-dpng', fig_title, '-r300');
  
end % for compressing_this % FullSim Time Domain

for compressing_this = 1 % continue; % FullSim Freq Domain
  
  call_sign = 'FullSim Freq Domain';
  fig_title = sprintf('%s_%s',mfilename, call_sign);
  fig_h = figure('Name',fig_title);
  subplot(131);
  imagesc(y, freqs/1e6, ftmp-max(ftmp(:)) ,[-30,0] );
  cb = colorbar; cb.Label.String = 'Relative Power, dB';
  grid on; hold on; axis tight;
  xlabel('Along-track position, rlines'); ylabel('Freq, MHz');
  title('fft(PhaseHistory) Magnitude');
  subplot(132);
  imagesc(y, freqs/1e6, angle(fdata) );
  cb = colorbar; cb.Label.String = 'Radians';
  grid on; hold on; axis tight;
  xlabel('Along-track position, rlines'); ylabel('Freq, MHz');
  title('fft(PhaseHistory) Phase');
  subplot(133);
  imagesc(y, freqs/1e6, unwrap(angle(fdata)) );
  cb = colorbar; cb.Label.String = 'Radians';
  grid on; hold on; axis tight;
  xlabel('Along-track position, rlines'); ylabel('Freq, MHz');
  title('fft(PhaseHistory) Phase Unwrap');
  try
    sgtitle(fig_title,'FontWeight','Bold','FontSize',14,'Interpreter','None');
  end
  set(findobj(fig_h,'type','axes'),'FontWeight', 'Bold', 'FontSize',14);
  set(fig_h, 'Position', get(0, 'Screensize'));
  print(gcf, '-dpng', fig_title, '-r300');
  
end % for compressing_this % FullSim Freq Domain

for compressing_this = 1 % continue; % FullSim TD FD windows
  
  call_sign = 'FullSim TD FD windows';
  fig_title = sprintf('%s_%s',mfilename, call_sign);
  fig_h = figure('Name',fig_title);
  subplot(121);
  plot(time/1e-6, tmp(:,twtt_min_idx) ); hold on;
  plot(twtt_min/1e-6, tmp(tmp_max_idx(twtt_min_idx), twtt_min_idx) , 's','LineWidth',2);
  xlabel('Fast-time, us');
  ylabel('Magnitude, dB');
  legend('rline','expected');
  grid on;
  title('At closest range');
  subplot(122);
  plot(freqs/1e6, ftmp(:,twtt_min_idx) ); hold on;
  xlabel('Freq, MHz');
  ylabel('Magnitude, dB');
  grid on;
  title('At closest range');
  try
    sgtitle(fig_title,'FontWeight','Bold','FontSize',14,'Interpreter','None');
  end
  set(findobj(fig_h,'type','axes'),'FontWeight', 'Bold', 'FontSize',14);
  set(fig_h, 'Position', get(0, 'Screensize'));
  print(gcf, '-dpng', fig_title, '-r300');
  
end % for compressing_this % FullSim TD FD windows

for compressing_this = 1 % continue; % FullSim Time Phase checks
  
  call_sign = 'FullSim Time Phase checks';
  fig_title = sprintf('%s_%s',mfilename, call_sign);
  fig_h = figure('Name',fig_title);
  subplot(311);
  plot(y, twtt/1e-6); hold on;
  plot(y, time(data_max_idx)/1e-6)
  xlabel('Along-track position, rlines');
  ylabel('Fast-time, us');
  legend('TWTT', 'fast-time(max(data))');
  grid on; axis tight;
  title('Time checks');
  subplot(312);
  plot(y, unwrap(angle(data_max)),'x'); hold on;
  plot(y, unwrap(angle(expected_phase)),'s');
  xlabel('Along-track position, rlines');
  ylabel('Radians');
  legend('UnwrapPhase(max(data))', 'Expected');
  grid on; axis tight;
  title('Phase Checks');
  subplot(313);
  plot(y, angle(data_max),'x'); hold on;
  plot(y, angle(expected_phase),'s');
  xlabel('Along-track position, rlines');
  ylabel('Radians');
  legend('Phase(max(data))', 'Expected');
  grid on; axis tight;
  title('Phase Checks');
  try
    sgtitle(fig_title,'FontWeight','Bold','FontSize',14,'Interpreter','None');
  end
  set(findobj(fig_h,'type','axes'),'FontWeight', 'Bold', 'FontSize',14);
  set(fig_h, 'Position', get(0, 'Screensize'));
  print(gcf, '-dpng', fig_title, '-r300');
  
end % for compressing_this % FullSim Time Phase checks
