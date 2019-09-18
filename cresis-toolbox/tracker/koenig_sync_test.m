if 1
  dk=load('/cresis/snfs1/dataproducts/public/data/temp/internal_layers/NASA_OIB_test_files/image_files/greenland_picks_final_2009-2012_20140602/2012/data_mat/layers_20120330_04_dec09.mat');
  new_layer = [];
  for lay_idx = 1:30
    [~,y] = max(dk.arr_layers==lay_idx);
    y(y==1) = NaN; % Assume 1-index max index implies no layer in this column
    new_layer(lay_idx,:) = y;
  end
  dk.time = dk.time_trace(1)*(0:size(dk.data_out,1)-1).';
  dk.time = dk.time - dk.time(205);
  new_layer_twtt = interp1(1:size(dk.data_out,1),dk.time,new_layer);
  
%   echo_fn = '/cresis/snfs1/dataproducts/public/data/snow/2012_Greenland_P3/CSARP_qlook/20120330_04/Data_20120330_04_174.mat';
%   lay_fn = '/cresis/snfs1/dataproducts/public/data/temp/internal_layers/NASA_OIB_test_files/snow/2012_Greenland_P3/CSARP_layerData/20120330_04/Data_20120330_04_174.mat';
  echo_fn = '/cresis/snfs1/dataproducts/public/data/snow/2012_Greenland_P3/CSARP_qlook/20120330_04/Data_20120330_04_190.mat';
  lay_fn = '/cresis/snfs1/dataproducts/public/data/temp/internal_layers/NASA_OIB_test_files/snow/2012_Greenland_P3/CSARP_layerData/20120330_04/Data_20120330_04_190.mat';
  mdata = load_L1B_L2(echo_fn,lay_fn);
  mdata.lay = mdata.lay([1 3:end],:);
  
  % Truncate datasets to the overlapping region
  gps_mask = dk.time_gps >= mdata.GPS_time(1) & dk.time_gps <= mdata.GPS_time(end);
  if ~any(gps_mask)
    gps_mask(:) = true;
  end
  dk.lat = dk.lat(gps_mask);
  dk.lon = dk.lon(gps_mask);
  dk.data_out = dk.data_out(:,gps_mask);
  dk.time_trace = dk.time_trace(gps_mask);
  dk.time_gps = dk.time_gps(gps_mask);
  dk.arr_segs = dk.arr_segs(:,gps_mask);
  dk.arr_layers = dk.arr_layers(:,gps_mask);
  new_layer_twtt = new_layer_twtt(:,gps_mask);
  
  if any(dk.time_trace~=dk.time_trace(1))
    error('Script does not handle changes in time_trace.');
  end
  mdata.Time_Surf = mdata.Time - mdata.Time(55);
  dk.data_out = interp1(dk.time, dk.data_out, mdata.Time_Surf);
  time_nan_mask = ~all(isnan(dk.data_out).').';
  if any(time_nan_mask)
    time_nan_idxs = find(time_nan_mask);
    dk.time = mdata.Time_Surf(time_nan_mask);
    dk.data_out = dk.data_out(time_nan_idxs,:);
    mdata.Data = mdata.Data(time_nan_idxs,:);
    mdata.Time = mdata.Time_Surf(time_nan_idxs);
    new_layer = interp1(dk.time,1:size(dk.data_out,1),new_layer_twtt);
  end
  
  figure(1); clf;
  plot(dk.lon,dk.lat,'LineWidth',2);
  hold on;
  grid on;
  
  figure(1);
  plot(mdata.Longitude,mdata.Latitude);
  hold on;
  
  figure(2); clf;
  imagesc(10*log10(mdata.Data));
  colormap(1-gray(256));
  hold on;
  %   plot(mdata.lay.');
  for lay_idx = 1:size(mdata.lay,1)
    plot(mdata.lay(lay_idx,:));
  end
  
  figure(3); clf;
  imagesc(dk.data_out);
  colormap(1-gray(256));
  hold on;
  %   plot(new_layer.');
  for lay_idx = 1:size(new_layer,1)
    plot(new_layer(lay_idx,:));
  end
  
  link_figures([2 3],'y')
  
  
elseif 0
  %Data_20120404_01_072
  dk=load('/cresis/snfs1/dataproducts/public/data/temp/internal_layers/NASA_OIB_test_files/image_files/greenland_picks_final_2009-2012_20140602/2012/data_mat/layers_20120404_01_dec08.mat');
  new_layer = [];
  for lay_idx = 1:30
    [~,y] = max(dk.arr_layers==lay_idx);
    y(y==1) = NaN; % Assume 1-index max index implies no layer in this column
    new_layer(lay_idx,:) = y;
  end
  dk.time = dk.time_trace(1)*(0:size(dk.data_out,1)-1).';
  dk.time = dk.time - dk.time(200);
  new_layer_twtt = interp1(1:size(dk.data_out,1),dk.time,new_layer);
  
  echo_fn = '/cresis/snfs1/dataproducts/public/data/snow/2012_Greenland_P3/CSARP_qlook/20120404_01/Data_20120404_01_072.mat';
  lay_fn = '/cresis/snfs1/dataproducts/public/data/temp/internal_layers/NASA_OIB_test_files/snow/2012_Greenland_P3/CSARP_layerData/20120404_01/Data_20120404_01_072.mat';
  mdata = load_L1B_L2(echo_fn,lay_fn);
  mdata.lay = mdata.lay([1 3:end],:);
  
  % Truncate datasets to the overlapping region
  gps_mask = dk.time_gps >= mdata.GPS_time(1) & dk.time_gps <= mdata.GPS_time(end);
  if ~any(gps_mask)
    gps_mask(:) = true;
  end
  dk.lat = dk.lat(gps_mask);
  dk.lon = dk.lon(gps_mask);
  dk.data_out = dk.data_out(:,gps_mask);
  dk.time_trace = dk.time_trace(gps_mask);
  dk.time_gps = dk.time_gps(gps_mask);
  dk.arr_segs = dk.arr_segs(:,gps_mask);
  dk.arr_layers = dk.arr_layers(:,gps_mask);
  new_layer_twtt = new_layer_twtt(:,gps_mask);
  
  if any(dk.time_trace~=dk.time_trace(1))
    error('Script does not handle changes in time_trace.');
  end
  mdata.Time_Surf = mdata.Time - mdata.Time(55);
  dk.data_out = interp1(dk.time, dk.data_out, mdata.Time_Surf);
  time_nan_mask = ~all(isnan(dk.data_out).').';
  time_nan_idxs = find(time_nan_mask);
  dk.time = mdata.Time_Surf(time_nan_mask);
  dk.data_out = dk.data_out(time_nan_idxs,:);
  mdata.Data = mdata.Data(time_nan_idxs,:);
  mdata.Time = mdata.Time_Surf(time_nan_idxs);
  new_layer = interp1(dk.time,1:size(dk.data_out,1),new_layer_twtt);
  
  
  figure(1); clf;
  plot(dk.lon,dk.lat,'LineWidth',2);
  hold on;
  grid on;
  
  figure(1);
  plot(mdata.Longitude,mdata.Latitude);
  hold on;
  
  figure(2); clf;
  imagesc(10*log10(mdata.Data));
  colormap(1-gray(256));
  hold on;
  %   plot(mdata.lay.');
  for lay_idx = 1:size(mdata.lay,1)
    plot(mdata.lay(lay_idx,:));
  end
  
  figure(3); clf;
  imagesc(dk.data_out);
  colormap(1-gray(256));
  hold on;
  %   plot(new_layer.');
  for lay_idx = 1:size(new_layer,1)
    plot(new_layer(lay_idx,:));
  end
  
  %   link_figures([2 3],'y')
  
elseif 0
  dk=load('/cresis/snfs1/dataproducts/public/data/temp/internal_layers/NASA_OIB_test_files/image_files/greenland_picks_final_2009-2012_20140602/2012/data_mat/layers_20120330_04_dec08.mat');
  new_layer = [];
  for lay_idx = 1:30
    [~,y] = max(dk.arr_layers==lay_idx);
    y(y==1) = NaN; % Assume 1-index max index implies no layer in this column
    new_layer(lay_idx,:) = y;
  end
  dk.time = dk.time_trace(1)*(0:size(dk.data_out,1)-1).';
  dk.time = dk.time - dk.time(200);
  new_layer_twtt = interp1(1:size(dk.data_out,1),dk.time,new_layer);
  
  echo_fn = '/cresis/snfs1/dataproducts/public/data/snow/2012_Greenland_P3/CSARP_qlook/20120330_04/Data_20120330_04_153.mat';
  lay_fn = '/cresis/snfs1/dataproducts/public/data/temp/internal_layers/NASA_OIB_test_files/snow/2012_Greenland_P3/CSARP_layerData/20120330_04/Data_20120330_04_153.mat';
  mdata = load_L1B_L2(echo_fn,lay_fn);
  mdata.lay = mdata.lay([1 3:end],:);
  
  % Truncate datasets to the overlapping region
  gps_mask = dk.time_gps >= mdata.GPS_time(1) & dk.time_gps <= mdata.GPS_time(end);
  dk.lat = dk.lat(gps_mask);
  dk.lon = dk.lon(gps_mask);
  dk.data_out = dk.data_out(:,gps_mask);
  dk.time_trace = dk.time_trace(gps_mask);
  dk.time_gps = dk.time_gps(gps_mask);
  dk.arr_segs = dk.arr_segs(:,gps_mask);
  dk.arr_layers = dk.arr_layers(:,gps_mask);
  new_layer_twtt = new_layer_twtt(:,gps_mask);
  
  if any(dk.time_trace~=dk.time_trace(1))
    error('Script does not handle changes in time_trace.');
  end
  mdata.Time_Surf = mdata.Time - mdata.Time(55);
  dk.data_out = interp1(dk.time, dk.data_out, mdata.Time_Surf);
  time_nan_mask = ~all(isnan(dk.data_out).').';
  time_nan_idxs = find(time_nan_mask);
  dk.time = mdata.Time_Surf(time_nan_mask);
  dk.data_out = dk.data_out(time_nan_idxs,:);
  mdata.Data = mdata.Data(time_nan_idxs,:);
  mdata.Time = mdata.Time_Surf(time_nan_idxs);
  new_layer = interp1(dk.time,1:size(dk.data_out,1),new_layer_twtt);
  
  
  figure(1); clf;
  plot(dk.lon,dk.lat,'LineWidth',2);
  hold on;
  grid on;
  
  figure(1);
  plot(mdata.Longitude,mdata.Latitude);
  hold on;
  
  figure(2); clf;
  imagesc(10*log10(mdata.Data));
  colormap(1-gray(256));
  hold on;
  %   plot(mdata.lay.');
  for lay_idx = 1:size(mdata.lay,1)
    plot(mdata.lay(lay_idx,:));
  end
  
  figure(3); clf;
  imagesc(dk.data_out);
  colormap(1-gray(256));
  hold on;
  %   plot(new_layer.');
  for lay_idx = 1:size(new_layer,1)
    plot(new_layer(lay_idx,:));
  end
  
  link_figures([2 3],'y')
  
  
end