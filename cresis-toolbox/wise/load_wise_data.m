path = 'F:\data\WISE\Helheim_2008\20080507T140219\';
fn.data = 'tone_2p5_2c_1200s_20080507T140219.ch0.cmp.uniform_alt_corr.img.dat';
fn.gps = 'tone_2p5_2c_1200S_20080507T140219.datgps.dat.CSV';
fn.bed = 'bed_20080507T140219.txt';

path = '/cresis/data2/WISE/Kanger_2009/20090409T073404';
fn.data = 'chirp_0.8v_1.8_3.8._1200s_20090409T073404.ch0.cmp.uniform_alt_corr.img.dat';
fn.gps = 'chirp_0.8v_1.8_3.8._1200s_20090409T073404.datgps.dat.CSV';
fn.bed = 'bed_20090409T073404.txt';

param.year = 2009;
param.month = 04;
param.day = 09;

% Load GPS data
raw = dlmread(fullfile(path,fn.gps));
lat = raw(:,2);
lon = raw(:,3);
elev = raw(:,4);

figure(1); clf;
plot(lon,lat);
xlabel('Longitude, E (deg)');
ylabel('Latitude, N (deg)');

figure(2); clf;
plot(elev);
ylabel('Elevation (m)');

% Load surface and bed picks
bed_data = load(fullfile(path,fn.bed));
bed_data(bed_data == -9999) = NaN;

figure(3); clf;
hs = plot(bed_data(:,2),'r');
hold on;
hb = plot(bed_data(:,3),'b');
hold off;
legend([hs hb], 'Surface', 'Bed');
grid on;

% Load radar echogram
fid = fopen(fullfile(path,fn.data),'r');

Nx = fread(fid,1,'long');
Ny = fread(fid,1,'long');

xAxis = fread(fid,Nx,'float32');
yAxis = fread(fid,Ny,'float32');

data = fread(fid,[Ny Nx],'float32');

fclose(fid);

figure(4); clf;
%imagesc(lp(local_detrend(flipud(data),[1 1], [5 10], 1)));
imagesc(lp(flipud(data)));
set(gca,'YDir','normal');
colormap(gray(256));
hold on;
hs = plot(bed_data(:,2),'r');
hb = plot(bed_data(:,3),'b');
hold off;
legend([hs hb], 'Surface', 'Bed');

return;

