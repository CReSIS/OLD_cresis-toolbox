%% Setup
[B,A] = butter(2,1/50);
figure(1000); clf;
cdata.t0 = [];
cdata.dt = [];
cdata.Data = [];
cdata.GPS_time = [];
cdata.layer = [];
max_time = 200e-9;
B = tukeywin(51,0.2).';
B = B / sum(B);
frms = 231:253;
frms = 100:120;

%% Load Layer Data
lay_all.GPS_time = [];
lay_all.surf_twtt = [];
for frm = frms
  %lay = load(sprintf('/cresis/snfs1/dataproducts/ct_data/snow/2011_Greenland_P3/CSARP_layerData/20110329_01/Data_20110329_01_%03d.mat',frm));
  lay = load(sprintf('/cresis/snfs1/dataproducts/ct_data/snow/2012_Greenland_P3/CSARP_layerData/20120330_04/Data_20120330_04_%03d.mat',frm));
  lay_all.GPS_time(end+(1:length(lay.GPS_time))) = lay.GPS_time;
  lay_all.surf_twtt(end+(1:length(lay.GPS_time))) = lay.layerData{1}.value{2}.data;
end
% Convert surface to rows/bins and smooth
lay_all.surf_twtt_filtered = fir_dec(lay_all.surf_twtt,B,1);

%% Load Data
for frm = frms
  %mdata = load_L1B(sprintf('/cresis/snfs1/dataproducts/ct_data/snow/2011_Greenland_P3/CSARP_post/CSARP_qlook/20110329_01/Data_20110329_01_%03d.mat',frm));
  mdata = load_L1B(sprintf('/cresis/snfs1/dataproducts/ct_data/snow/2012_Greenland_P3/CSARP_post/CSARP_qlook/20120330_04/Data_20120330_04_%03d.mat',frm));
  dt = mdata.Time(2)-mdata.Time(1);
  t0 = mdata.Time(1);
  %lay = load(sprintf('/cresis/snfs1/dataproducts/ct_data/snow/2011_Greenland_P3/CSARP_layerData/20110329_01/Data_20110329_01_%03d.mat',frm));
  lay = load(sprintf('/cresis/snfs1/dataproducts/ct_data/snow/2012_Greenland_P3/CSARP_layerData/20120330_04/Data_20120330_04_%03d.mat',frm));
  
  surf_twtt_filtered_bin = (interp1(lay_all.GPS_time,lay_all.surf_twtt_filtered,mdata.GPS_time)-t0)/dt - 50;
  
  if 0
    clf;
    plot(surf_twtt_filtered_bin);
    hold on;
    plot((lay.layerData{1}.value{2}.data-t0) / dt - 50);
  end
  
  % Flatten data with smoothed surface
  mdata.Data(isnan(mdata.Data)) = 0;
  mdata.Data = fft(mdata.Data);
  Nt = size(mdata.Data,1);
  freq = ifftshift( 1/Nt*(-floor(Nt/2) : floor((Nt-1)/2)) ).';
  for rline=1:size(mdata.Data,2)
    mdata.Data(:,rline) = mdata.Data(:,rline) .* exp(1i*2*pi*freq*surf_twtt_filtered_bin(rline));
  end
  mdata.Data = abs(ifft(mdata.Data));
  
  % Convert layers to rows/bins and flatten with smoothed surface
  layer = [];
  for idx = 1:length(lay.layerData)
    if idx == 1
%       layer(idx,:) = 50*ones(size(lay.layerData{idx}.value{2}.data));
      layer(idx,:) = (interp1(lay.GPS_time,lay.layerData{idx}.value{2}.data,mdata.GPS_time)-t0)/dt - surf_twtt_filtered_bin;
    else
      layer(idx,:) = (interp1(lay.GPS_time,lay.layerData{idx}.value{2}.data,mdata.GPS_time)-t0)/dt - surf_twtt_filtered_bin;
    end
  end

  % Filter data
  mdata.Data = fir_dec(mdata.Data,B,1);
  
  % Sampling
  if isempty(cdata.Data)
    cdata.dt = dt;
    Nt_truncate = round(max_time/cdata.dt) + 50;
  end
  if dt ~= cdata.dt
    keyboard
  end
  
  % Truncate
  mdata.Data = mdata.Data(1:Nt_truncate,:);
  
  % Concatenate results
  cdata.t0 = surf_twtt_filtered_bin - 50*dt;
  cdata.Data = cat(2,cdata.Data,mdata.Data);
  cdata.GPS_time = cat(2,cdata.GPS_time,mdata.GPS_time);
  cdata.layer = cat(2,cdata.layer,layer);
  
%   clf
%   imagesc(10*log10(cdata.Data))
%   hold on;
%   plot(cdata.layer.');
%   keyboard
end


clf
imagesc(10*log10(cdata.Data))
hold on;
plot(cdata.layer.');




