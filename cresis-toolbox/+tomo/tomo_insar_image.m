if 0
  fn = '/cresis/snfs1/dataproducts/ct_data/rds/2011_Greenland_P3/CSARP_insar/rds_thule_2011_2014_insar.mat';
  fn = '/cresis/snfs1/dataproducts/ct_data/rds/2011_Greenland_P3/CSARP_insar/rds_thule_2011_2014_sameday_insar.mat';
  set1 = 1:15;
  set2 = 16:30;
  title_str = '2011 to 2014';
  time_year = 3;
  
elseif 0
  fn = '/cresis/snfs1/dataproducts/ct_data/rds/2012_Greenland_P3/CSARP_insar/rds_thule_2012_2014_wf1_insar.mat';
  set1 = 1:15;
  set2 = 16:30;
  title_str = '2012 to 2014';
  time_year = 2;
elseif 0
  fn = '/cresis/snfs1/dataproducts/ct_data/rds/2013_Greenland_P3/CSARP_insar/rds_thule_2013_2014_wf2_insar.mat';
  set1 = 1:15;
  set2 = 16:22;
  title_str = '2013 to 2014';
  time_year = 1;
elseif 1
  %   fn = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/rds_thule_combine_wf2_insar.mat';
  fn = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/rds_thule_combine_wf2_sameday_insar.mat';
  %   fn = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/rds_thule_combine_wf2_month_insar.mat';
  %   fn = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/rds_thule_combine_wf2_2week_insar.mat';
  set1 = 1:15;
  set2 = 16:30; % none
  %   set2 = 31:45; % Tiny
  %   set2 = 46:60; % Tiny
  %   set2 = 61:75; % Tiny
  %   set2 = 76:90; % Tiny
  %   set2 = 91:105; % Biggest (snow fall?)
  title_str = '2014 to 2014 (same day)';
  time_year = 1;
elseif 0
  %   fn = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/rds_thule_combine_wf2_insar.mat';
  %   fn = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/rds_thule_combine_wf2_sameday_insar.mat';
  %     fn = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/rds_thule_combine_wf2_month_insar.mat';
  fn = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/rds_thule_combine_wf2_2week_insar.mat';
  set1 = 1:15;
  set2 = 16:30; % none
  %   set2 = 31:45; % Tiny
  %   set2 = 46:60; % Tiny
  %   set2 = 61:75; % Tiny
  %   set2 = 76:90; % Tiny
  %   set2 = 91:105; % Biggest (snow fall?)
  title_str = '2014 to 2014 (1 month)';
  title_str = '2014 to 2014 (2 week)';
  time_year = 1;
end
rbins = [240:420]; % Good range 
% rbins = [160:440]; % Test range

rbins_baseline = [320:350]; % Bins to use to estimate phase error due to GPS errors
rbins_baseline = rbins_baseline - rbins(1);

tmp = load(fn);

%% Process data
insar1 = mean(tmp.data(rbins,:,set1),3);
insar2 = mean(tmp.data(rbins,:,set2),3);

insar_data = insar2 .* conj(insar1) ./ (abs(insar1).*abs(insar2));

% Surface flattening
% Remove surface variations by moving all bins down so that the surface
% aligns for every range line and the surface bin is the max surface bin
max_surface_bin = ceil(max(tmp.ref.surface_bin));
dbin = tmp.ref.surface_bin - max_surface_bin;
Nt = size(insar_data,1);
freq_bin = 1i*2*pi/Nt * ifftshift( -floor(Nt/2) : floor((Nt-1)/2) ).';
insar_data = ifft(fft(insar_data) .* exp(bsxfun(@times,dbin,freq_bin)));

% Remove bad circular deconvolution bins from surface flattening
[max_dbin,max_dbin_rline] = max(abs(dbin));
max_dbin = ceil(max_dbin);
insar_data = insar_data(1+max_dbin:end,:);
rbins = rbins(1+max_dbin:end);

% Baseline phase error estimation (caused by GPS position errors)
% baseline_correction = fir_dec(mean(insar_data(1:30,:),1),ones(1,51)/51,1);
baseline_correction = fir_dec(mean(insar_data(80:110,:),1),ones(1,51)/51,1);
insar_data = bsxfun(@times,insar_data,exp(-1i*angle(baseline_correction)));

insar_data_filt = fir_dec(fir_dec(insar_data,ones(1,31)/31,1).',ones(1,3)/3).';

meter_units = true;

range = (tmp.ref.wfs(tmp.ref.wf).time(rbins) - tmp.ref.wfs(tmp.ref.wf).time(max_surface_bin)) * c/2 / sqrt(er_ice);

along_track = geodetic_to_along_track(tmp.ref.lat,tmp.ref.lon);

figure(3001); clf;
if meter_units
  imagesc(along_track/1e3,range,hsv_plot(insar_data_filt,-2));
else
  imagesc(hsv_plot(insar_data_filt,-2));
end
colormap(hsv(256))
h_colorbar = colorbar;
caxis([-pi pi])
set(get(h_colorbar,'ylabel'),'string','Angle (radians)');
grid on;
if meter_units
  xlabel('Along-track (km)');
  ylabel('Depth, \epsilon_{r,ice}=3.15 (m)');
else
  xlabel('Range line');
  ylabel('Range bin');
end
title(title_str);
set(3001,'Position',[56   470   867   426]);

figure(3000); clf;
if meter_units
  imagesc(along_track/1e3,range,abs(insar_data_filt));
else
  imagesc(abs(insar_data_filt));
end
h_colorbar = colorbar;
set(get(h_colorbar,'ylabel'),'string','Coherence');
grid on;
if meter_units
  xlabel('Along-track (km)');
  ylabel('Depth, \epsilon_{r,ice}=3.15 (m)');
else
  xlabel('Range line');
  ylabel('Range bin');
end
title(title_str);
set(3000,'Position',[56   470   867   426]);
caxis([0 1]);

[fn_dir,fn_name] = fileparts(fn);
fn_insar_phase = fullfile(fn_dir,[fn_name '_phase.fig']);
fprintf('Saving %s (%s)\n', fn_insar_phase, datestr(now));
saveas(3001,fn_insar_phase);
fn_insar_phase = fullfile(fn_dir,[fn_name '_phase.png']);
fprintf('Saving %s (%s)\n', fn_insar_phase, datestr(now));
saveas(3001,fn_insar_phase);

[fn_dir,fn_name] = fileparts(fn);
fn_insar_coherence = fullfile(fn_dir,[fn_name '_coherence.fig']);
fprintf('Saving %s (%s)\n', fn_insar_coherence, datestr(now));
saveas(3000,fn_insar_coherence);
fn_insar_coherence = fullfile(fn_dir,[fn_name '_coherence.png']);
fprintf('Saving %s (%s)\n', fn_insar_coherence, datestr(now));
saveas(3000,fn_insar_coherence);

link_figures([3000 3001]);

%%

figure(3003);
% clf;
Nt = size(insar_data_filt,1);
Mt = 200;
dd=abs(fftshift(fft(insar_data_filt,Mt*Nt),1)).^2;
dd=fir_dec(dd,ones(1,501)/501,20);
along_track_dec=fir_dec(along_track,ones(1,501)/501,20);

% imagesc(lp(dd));

% time_years = (tmp.pass(set1(1)).gps_time(1) - tmp.pass(set2(1)).gps_time(1))/86400/365.25;

[dd_val,dd_idx] = max(dd);
max_freq = 1 + dd_idx - Nt*Mt/2;

lambda = 3e8/195e6;
dt = tmp.ref.wfs(tmp.ref.wf).time(2)-tmp.ref.wfs(tmp.ref.wf).time(1);
dr = dt*3e8/2;
if meter_units
  xlabel('Along-track (km)');
  ylabel('Vertical strain rate (m\cdotm^{-1}\cdota^{-1})');
  plot(along_track_dec/1e3, 2*pi*max_freq/(Nt*Mt) * lambda/(4*pi) / dr / time_year)
else
  xlabel('Range line');
  ylabel('Radians/range bin');
  plot(2*pi*max_freq/(Nt*Mt) * lambda/(4*pi) / dr / time_year)
end
grid on;
hold on;
xlim(along_track_dec([1 end])/1e3);

mean(2*pi*max_freq/(Nt*Mt) * lambda/(4*pi) / dr / time_year)

mean(2*pi*max_freq(end-90:end-10)/(Nt*Mt) * lambda/(4*pi) / dr / time_year)

strain_rate = 2*pi*max_freq/(Nt*Mt) * lambda/(4*pi) / dr / time_year;
