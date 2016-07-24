function [hdr,data] = load_fmcw_bin(fn)
%
% Reader for old CReSIS custom binary FMCW products (2009-2010 OIB seasons).
%
% fn: string containing filename to load
%
% hdr: header structure with header fields for each range line
%  .utc_time_sod: UTC time (seconds of day)
%  .lat: latitude (deg, N)
%  .lon: longitude (deg, E)
%  .elev: elevation (m)
%  .cshift: circular shift (applied to data)
%  .dt: range bin size in two way travel time (sec)
% data: linear power data
% 
% Data files contain a continuous stream of records. The first record
% begins at the beginning of the file. The last record ends at the end of
% the file. In other words records do not span across files.
%
% Record Format:
%  1. 28 byte header with 7 fields
%  2. Radar range line data samples (linear power)
%
% Record Format Detailed:
%  Number of bytes in record (int32)
%    Number of samples in record = (Number of bytes in record - 28)/4
%  UTC time*1e3 (int32)
%  Latitude*1e6 (int32)
%  Longitude*1e6 (int32)
%  Platform-Elevation*1e3 (int32)
%  Circular Shifts (int32)
%  Fast time bin size in one way travel time*1e12, dt (float32)
%  Data samples (Number of samples in record * float32)
%
% Example:
%   fn = 'data00.0096.bin';
%   [hdr,data] = load_fmcw_bin(fn);
%   
% Author: John Paden

% Preparation
hdr = [];
HEADER_SIZE_IN_BYTES = 28;
[fid,msg] = fopen(fn,'r','ieee-le');
if fid<0
  error('Could not open %s: %s.', fn, msg);
end

%% Determine the record size
rec_size_bytes = fread(fid,1,'int32');
num_sam = (rec_size_bytes-HEADER_SIZE_IN_BYTES)/4;

% Determine size of file and number of records
fseek(fid,0,1);
file_size = ftell(fid);
num_rec = floor(file_size / num_sam);

%% Read in file
fseek(fid,0,-1);
hdr_data = fread(fid,[6 num_rec],'6*int32',rec_size_bytes-6*4);
hdr.utc_time_sod = hdr_data(2,:)*1e-3;
hdr.lat = hdr_data(3,:)*1e-6;
hdr.lon = hdr_data(4,:)*1e-6;
hdr.elev = hdr_data(5,:)*1e-3;
hdr.cshift = hdr_data(6,:);
clear hdr_data;

fseek(fid,24,-1);
dt_and_data = fread(fid,[1+num_sam num_rec],sprintf('%d*float32',num_sam+1),6*4);
hdr.dt = dt_and_data(1,:)*2e-12;
data = dt_and_data(2:end,:);
clear dt_data;

fclose(fid);

%% Unfo circular shift
for rline = 1:size(data,2)
  data(:,rline) = circshift(data(:,rline),-hdr.cshift(rline));
end

time = hdr.dt(1)*(0:size(data,1)-1);

%% Plot results
figure(1); clf;
imagesc([],time*1e6,db(data));
xlabel('Range line');
ylabel('Relative fast time (us)');
cc = caxis; caxis(cc(2)+[-45 0]);
colormap(1-gray(256));

figure(2); clf;
imagesc([],time*3e8/2,db(data));
xlabel('Range line');
ylabel('Relative range (m)');
cc = caxis; caxis(cc(2)+[-45 0]);
colormap(1-gray(256));

figure(3); clf;
plot(hdr.utc_time_sod);
xlabel('Range line');
ylabel('UTC time (seconds of day)');
grid on;

figure(4); clf;
plot(hdr.lon, hdr.lat);
title('Position')
xlabel('Longitude (deg, E)');
ylabel('Latitude (deg, N)');
grid on;

figure(5); clf;
plot(hdr.elev);
xlabel('Range line');
ylabel('Elevation (m)');
grid on;
