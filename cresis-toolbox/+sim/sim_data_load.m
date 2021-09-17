
physical_constants;

% param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/snow/2012_Greenland_P3sim/20120330/param.mat';

%
param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/snow/2013_Greenland_P3sim/20130327/param.mat';
% param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/snow/2013_Greenland_P3sim/20130327/param.mat';
% % %
% param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/snow/2018_Antarctica_DC8sim/20181010/param.mat';

% param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/rds/2014_Greenland_P3sim/20140325/param.mat';

load(param_fn);
if ~strcmpi(param.sim.out_fn_param,param_fn); return; fprintf('hmm..something is not good\n'); end;

%%
param.load_data.pulse_comp            = 1;
param.load_data.raw_data              = 0;

[hdr,data] = load_data(param);

[Nt,Nx] = size(data{1});

%
[max_val, max_idx] = max(data{1},[],1);
% [closest_idx, closest_rline] = min(max_idx);


[traj.x, traj.y, traj.z] = geodetic2ecef(wgs84Ellipsoid,hdr.records{1}.lat,hdr.records{1}.lon,hdr.records{1}.elev);
% [traj.x, traj.y, traj.z] = geodetic2ecef(wgs84Ellipsoid,hdr.param_load_data.gps.lat,hdr.param_load_data.gps.lon,hdr.param_load_data.gps.elev);

target = hdr.param_load_data.target;

% range = 660 * ones(Nx,1);
range = sqrt( (traj.x - target.x).^2 + (traj.y - target.y).^2 + (traj.z - target.z).^2 );

td = 2*range/c;
[closest_td, closest_rline] = min(td);

v_aircraft = 125;
prf = hdr.param_load_data.radar.prf;
% fc = hdr.param_load_data.radar.wfs.fc;
fc = hdr.freq{1}(1);

chirp_rate = hdr.param_load_data.radar.wfs.chirp_rate;
pri = 1 / (prf/hdr.param_load_data.radar.wfs.presums );
lambda = c/fc;
wave_number = 4*pi/lambda;
sample_dist = 125 * ( pri );
expect = exp(-1i*wave_number*range)*max_val(closest_rline)/ exp(-1i*wave_number*range(closest_rline));

%%

fig_h = figure(3333); clf(3333);
subplot(241);imagesc(1:Nx,hdr.time{1}/1e-6,lp(data{1}));colorbar;
hold on;
plot(td/1e-6,'--','Color','k');
ylabel('Time,us');
title(param.sim.season_name,'Interpreter', 'none');

subplot(242);plot(hdr.time{1}/1e-6,lp(data{1}(:,ceil(size(data{1},2)/2))));
title(param.sim.day_seg,'Interpreter', 'none');

subplot(243);imagesc(unwrap(angle(data{1})));colorbar;
title('unwrapped phase');

subplot(244);
plot(unwrap(angle(max_val))); hold on;
plot(unwrap(angle(expect)));
xlabel('rlines');
axis tight;
title('unwrapped phase MAX');

subplot(248);
plot(td/1e-6); hold on;
plot(hdr.time{1}(max_idx)/1e-6)
xlabel('rlines');
axis tight;
title('tdelay, us');


subplot(245);imagesc(lp(fft(data{1})));colorbar;
title('fft');
subplot(246);plot(lp(fft(data{1}(:,ceil(size(data{1},2)/2)))));
title('fft');
subplot(247);imagesc(unwrap(angle(fft(data{1}))));colorbar;
title('fft');
%%
% figure(469); clf(469);
% plot(td/1e-6);hold on;
% plot(hdr.time{1}(max_idx)/1e-6)
% xlabel('rlines');
% axis tight;
% title('tdelay, us');
% % keyboard;


figure(2);clf(2);
plot((angle(max_val)),'-o'); hold on;
plot((angle(expect)),'-x');
plot(closest_rline,angle(expect(closest_rline)) ,'o','LineWidth',4,'Color','g');
xlim([closest_rline-80 closest_rline+80]);
grid on
% legend({});

test = param.sim.range - range;
if any(test);  warning('Please check range in simulator and loader'); end
test_td = param.sim.twtt - td;
if any(test_td)  warning('Please check td in simulator and loader'); end

%%

% fig_title = sprintf('load %d %s %s',param.load_data.pulse_comp, param.sim.season_name, param.sim.day_seg);
% print(fig_h, '-dpng', fig_title, '-r300');
% % % % % % saveas(fig_h,fig_title);
