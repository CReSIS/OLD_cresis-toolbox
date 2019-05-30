function mdata = load_L1B_L2(echo_fn,lay_fn)
% mdata = load_L1B_L2(echo_fn,lay_fn)
%
% Loads the Level 1B echogram file and Level 2 layer file specified, aligns
% the layer data to the echogram data, shifts each column so that the
% surface is at the ref_bin, low pass filters the echogram along each row,
% and returns the resulting echogram structure with layers, filtered data,
% and fields to undo the time shift.
%
% echo_fn: Echogram filename (e.g. CSARP_qlook)
% lay_fn: Layerdata filename (e.g. CSARP_layerData)
%
% mdata: Structure containing echogram fields with the following changes:
% .Data: shifted and filtered Nt by Nx echogram image
% .t0: 1 by Nx vector, two way travel time for the first row of each column
% .dt: scalar with time spacing between rows
% .lay: every layer from lay_fn shifted and aligned with .Data field, units
%   also converted to .Data row/bin.
%
% Example:
%
% echo_fn = '/cresis/snfs1/dataproducts/ct_data/snow/2012_Greenland_P3/CSARP_post/CSARP_qlook/20120330_04/Data_20120330_04_100.mat';
% lay_fn = '/cresis/snfs1/dataproducts/ct_data/snow/2012_Greenland_P3/CSARP_layerData/20120330_04/Data_20120330_04_100.mat';
% mdata = load_L1B_L2(echo_fn,lay_fn);
% figure; clf;
% imagesc(10*log10(mdata.Data));
% colormap(1-gray(256));
% hold on;
% plot(mdata.lay.');
%
% Author: John Paden

%% Setup
ref_bin = 50;
max_time = 200e-9;
detrend_poly_order = 4;
B = tukeywin(51,0.2).';
B = B / sum(B);

%% Load Data
mdata = load_L1B(echo_fn);
lay = load(lay_fn);

% dt: time increment from one row to the next
dt = mdata.Time(2)-mdata.Time(1);
% t0: two way travel time to first row of image
t0 = mdata.Time(1);

%% Flatten Data

% Convert surface to rows/bins relative to first row of image
lay_idx = 1;
surf_twtt_bin = (interp1(lay.GPS_time,lay.layerData{lay_idx}.value{2}.data,mdata.GPS_time)-t0)/dt - ref_bin;

% Flatten data with surface
mdata.Data(isnan(mdata.Data)) = 0;
mdata.Data = fft(mdata.Data);
Nt = size(mdata.Data,1);
freq = ifftshift( 1/Nt*(-floor(Nt/2) : floor((Nt-1)/2)) ).';
for rline=1:size(mdata.Data,2)
  mdata.Data(:,rline) = mdata.Data(:,rline) .* exp(1i*2*pi*freq*surf_twtt_bin(rline));
end
mdata.Data = abs(ifft(mdata.Data));

% Filter data
mdata.Data = fir_dec(mdata.Data,B,1);

% Truncate data
Nt_truncate = round(max_time/dt) + 50;
if Nt_truncate > size(mdata.Data,1)
  mdata.Data = mdata.Data(1:end,:);
else
  mdata.Data = mdata.Data(1:Nt_truncate,:);
end

% Convert layers to rows/bins and flatten with smoothed surface
mdata.lay = [];
for lay_idx = 1:length(lay.layerData)
  mdata.lay(lay_idx,:) = (interp1(lay.GPS_time,lay.layerData{lay_idx}.value{2}.data,mdata.GPS_time)-t0)/dt - surf_twtt_bin;
end

% Create remaining output fields
mdata.dt = dt;
mdata.t0 = t0 + surf_twtt_bin*dt;

%% Detrend Data
if detrend_poly_order > 0
  % Find the along-track mean of the power detected data
  mean_data = log(mean(mdata.Data,2));
  
  % Create a mask of good points to use in polynomial fitting of the data
  mask = false(size(mdata.Data,1),1);
  mask(ref_bin:size(mdata.Data,1)) = true;
  mask(~isfinite(mean_data)) = false;
  
  % Fit a polynomial to the trend of the data and remove this trend
  p = polyfit(find(mask),mean_data(mask),detrend_poly_order);
  trend = polyval(p,(1:size(mdata.Data)).');
  trend(~mask) = NaN;
  trend = interp_finite(trend,0);
  mdata.Data = bsxfun(@rdivide,mdata.Data,exp(trend));
end

%% Debug Plots
if 0
  figure(1); clf;
  imagesc(10*log10(mdata.Data));
  hold on;
  plot(mdata.lay.');
  
  figure(2); clf;
  plot(mdata.t0 + ref_bin*dt);
  hold on;
  lay_idx = 1;
  plot(interp1(lay.GPS_time,lay.layerData{lay_idx}.value{2}.data,mdata.GPS_time),'--');
end
