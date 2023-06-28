% script example_accum_reader
%
% Examples of
% 1. How to load OIB accumulation radar data posted at NSIDC
% for 2010 Greenland P-3 field season.
% 
% 2. Load 2011,2012,2013 Greenland P-3 field seasons and plot
% versus two way travel time, WGS-84 elevation, and range.
%
% Author: John Paden

if 0
  % ====================================================================
  % User settings
  % ====================================================================
  fn = 'D:\accum\accum.20100508A.4.151_200.dat';
  year = 2010;
  month = 5;
  day = 8;
  incoh_ave = 1;
  t0 = 2.296e-6; % found by looking at surface multiple
  
  % ====================================================================
  % Automated section
  % ====================================================================
  fid = fopen(fn,'r','ieee-le');
  data = fread(fid,[2549 inf],'single');
  fclose(fid);
  
  utc_time_sod = double(data(1,:));
  utc_time_fractions = double(data(2,:));
  data = data(3:end,:);
  
  utc_time = datenum(year,month,day,0,0,utc_time_sod+utc_time_fractions/1e5);
  fprintf('%s to %s\n', datestr(utc_time(1)), datestr(utc_time(end)));
  
  fs = 125e6;
  dt = 1/fs;
  Nt = size(data,1);
  time = (0:dt:(Nt-1)*dt).' + 2.296e-6;
  
  figure(1); clf;
  imagesc([],time*1e6,10*log10(filter2(ones(1,incoh_ave)/incoh_ave,data)));
  cur_caxis = caxis;
  xlabel('Along-track (range lines)');
  ylabel('Fast time (us)');
  colormap(gray(256));
  colorbar;
end

if 1
  c = 2.997924580003452e+008;
  er_ice = 3.15;
  mdata = load('Data_20111107_02_191.mat')
  mdata = uncompress_echogram(mdata);

  max_elev = max(mdata.Elevation);
  dRange = max_elev - mdata.Elevation;
  dt = mdata.Time(2)-mdata.Time(1);
  dBins = round(dRange / (c/2) / dt);
  dtime = dRange/(c/2);
  zero_pad_len = max(abs(dBins));
  mdata.Data = cat(1,mdata.Data,zeros(zero_pad_len,size(mdata.Data,2)));
  mdata.Time = mdata.Time(1) + (mdata.Time(2)-mdata.Time(1)) * (0:1:size(mdata.Data,1)-1);
  update_surf = true;
  mdata.Surface_Bin = zeros(size(mdata.Surface));
  for rline = 1:size(mdata.Data,2)
    mdata.Data(:,rline) = interp1(mdata.Time, mdata.Data(:,rline), mdata.Time - dtime(rline), 'linear',0);
    mdata.Elevation(rline) = mdata.Elevation(rline) + dRange(rline);
    % Use original surface:
    mdata.Surface(rline) = mdata.Surface(rline) + dtime(rline);
    % Update surface:
    mdata.Surface_Bin(rline) = find(lp(mdata.Data(:,rline)) > -40,1);
  end
  if update_surf
    % Update surface:
    mdata.Surface = interp1(1:length(mdata.Time),mdata.Time,mdata.Surface_Bin);
  else
    % Update surface bins:
    mdata.Surface_Bin = interp1(mdata.Time,1:length(mdata.Time),mdata.Surface_Bin);
  end

  % Set which bins to plot
  good_bins = min(mdata.Surface_Bin)-50 : max(mdata.Surface_Bin)+600;
  
  figure(1); clf;
  imagesc([],mdata.Time(good_bins)*1e6,lp(mdata.Data(good_bins,:)))
  xlabel('Range line');
  ylabel('Two way travel time (us)')
  colormap(1-gray(256))
  hold on
  plot(mdata.Surface*1e6);
  hold off;
  
  % Re-interpolate to WGS-84 elevation y-axis
  max_surf_bin = max(mdata.Surface_Bin);
  max_range = mdata.Time(max_surf_bin) * c/2 + (mdata.Time(end)-mdata.Time(max_surf_bin)) * c/2 / sqrt(er_ice);
  dz = dt * c / 2 / sqrt(er_ice);
  range_yaxis = 0:dz:max_range;
  mdata.Data_WGS84 = zeros(length(range_yaxis),size(mdata.Data,2));
  for rline = 1:size(mdata.Data,2)
    orig_range_yaxis = [mdata.Time(1:mdata.Surface_Bin(rline))*c/2  ...
      mdata.Time(mdata.Surface_Bin(rline))*c/2 + (mdata.Time(mdata.Surface_Bin(rline)+1:end) - mdata.Time(mdata.Surface_Bin(rline)))*c/2/sqrt(er_ice)];
    mdata.Data_WGS84(:,rline) = interp1(orig_range_yaxis, mdata.Data(:,rline), range_yaxis);
  end
  mdata.Surface_Range = mdata.Surface * c/2;
  
  % Set which ranges to plot
  good_bins = round(interp1(range_yaxis, 1:length(range_yaxis), [min(mdata.Surface_Range)-50 max(mdata.Surface_Range)+200]));
  good_bins = good_bins(1):good_bins(2);
  
  figure(2); clf;
  imagesc([],range_yaxis(good_bins),lp(mdata.Data_WGS84(good_bins,:)));
  xlabel('Range line');
  ylabel('Range (m)')
  colormap(1-gray(256))
  hold on
  plot(mdata.Surface*c/2);
  hold off;
  
  figure(3); clf;
  imagesc([],mdata.Elevation(1) - range_yaxis(good_bins),lp(mdata.Data_WGS84(good_bins,:)));
  xlabel('Range line');
  ylabel('WGS-84 (m)')
  set(gca,'YDir','normal')
  colormap(1-gray(256))
  hold on
  plot(mdata.Elevation(1) - mdata.Surface*c/2);
  hold off;
  
end