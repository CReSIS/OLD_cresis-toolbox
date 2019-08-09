fn = '';
if 0
  %% 2011 and 2014 Comparison
  insar_mode = 3; % 1 to find equalization coefficients, 2 to process data, 3 to differential SAR
  
  equalization_2011 = 10.^(zeros(1,15)/20) .* exp(1i*([-17.8 -14.4 -16.3 -15.2 -20.0 -19.4 -16.0 0.0 -10.4 -9.5 -12.2 -28.3 -31.6 -32.0 -19.4])/180*pi);
  equalization_2014 = 10.^(zeros(1,15)/20) .* exp(1i*([122.5 121.6 126.7 105.5 130.1 120.1 128.7 -0.0 -134.2 121.1 36.4 125.8 -171.5 -1.0 128.4]-[42.8 45.5 48.2 50.9 53.6 56.3 58.9 0.0 2.7 5.4 8.1 93.6 96.3 99.0 101.7]/2)/180*pi);
  % Enable for waveform 1
%   rbins = [];
%   if ispc
%     fn = fullfile('X:/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_2011_2014_wf1'));
%   else
%     fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2011_Greenland_P3/CSARP_insar/',sprintf('rds_thule_2011_2014_wf1'));
%   end
  % Enable for waveform 2
  rbins = 220:420;
  if ispc
    fn = fullfile('X:/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_2011_2014'));
  else
    fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2011_Greenland_P3/CSARP_insar/',sprintf('rds_thule_2011_2014'));
  end
  % Enable for waveform 3
%   rbins = [];
%   if ispc
%     fn = fullfile('X:/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_2011_2014_wf3'));
%   else
%     fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2011_Greenland_P3/CSARP_insar/',sprintf('rds_thule_2011_2014_wf3'));
%   end
  equalization = [equalization_2014 equalization_2011];
  
  master_idx = 8;

elseif 0
  %% 2012 and 2014 Comparison
  insar_mode = 3; % 1 to find equalization coefficients, 2 to process data, 3 to differential SAR
  
  equalization_2012 = 10.^(zeros(1,15)/20) .* exp(1i*([-8.0 -7.3 0.0 -10.8 -0.2 -8.4 1.5 -22.7 -26.3 -21.8 -13.9 11.6 12.5 14.1 14.0])/180*pi);
  equalization_2014 = 10.^(zeros(1,15)/20) .* exp(1i*([122.5 121.6 126.7 105.5 130.1 120.1 128.7 -0.0 -134.2 121.1 36.4 125.8 -171.5 -1.0 128.4]-[42.8 45.5 48.2 50.9 53.6 56.3 58.9 0.0 2.7 5.4 8.1 93.6 96.3 99.0 101.7]/2)/180*pi);
  % Enable for waveform 1
%   rbins = [];
%   if ispc
%     fn = fullfile('X:/ct_data/rds/2012_Greenland_P3/CSARP_insar/',sprintf('rds_thule_2012_2014_wf1'));
%   else
%     fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2012_Greenland_P3/CSARP_insar/',sprintf('rds_thule_2012_2014_wf1'));
%   end
  % Enable for waveform 2
  rbins = 220:420;
  if ispc
    fn = fullfile('X:/ct_data/rds/2012_Greenland_P3/CSARP_insar/',sprintf('rds_thule_2012_2014_wf2'));
  else
    fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2012_Greenland_P3/CSARP_insar/',sprintf('rds_thule_2012_2014_wf2'));
  end
  % Enable for waveform 3
%   rbins = [];
%   if ispc
%     fn = fullfile('X:/ct_data/rds/2012_Greenland_P3/CSARP_insar/',sprintf('rds_thule_2012_2014_wf3'));
%   else
%     fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2012_Greenland_P3/CSARP_insar/',sprintf('rds_thule_2012_2014_wf3'));
%   end
  equalization = [equalization_2014 equalization_2012];
  
  master_idx = 8;

elseif 1
  %% 2013 and 2014 Comparison
  insar_mode = 3; % 1 to find equalization coefficients, 2 to process data, 3 to differential SAR
  
    coregistration_time_shift = [zeros(1,15) 1.43*ones(1,7)];
  equalization_2013 = 10.^(zeros(1,7)/20) .* exp(1i*([-3.9 -5.0 0.0 -8.1 -0.5 -4.3 1.1])/180*pi);
  equalization_2014 = 10.^(zeros(1,15)/20) .* exp(1i*([122.5 121.6 126.7 105.5 130.1 120.1 128.7 -0.0 -134.2 121.1 36.4 125.8 -171.5 -1.0 128.4]-[42.8 45.5 48.2 50.9 53.6 56.3 58.9 0.0 2.7 5.4 8.1 93.6 96.3 99.0 101.7]/2)/180*pi);
  % Enable for waveform 1
%   rbins = [];
%   if ispc
%     fn = fullfile('X:/ct_data/rds/2013_Greenland_P3/CSARP_insar/',sprintf('rds_thule_2013_2014_wf1'));
%   else
%     fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2013_Greenland_P3/CSARP_insar/',sprintf('rds_thule_2013_2014_wf1'));
%   end
  % Enable for waveform 2
    rbins = [200:420];
  if ispc
    fn = fullfile('X:/ct_data/rds/2013_Greenland_P3/CSARP_insar/',sprintf('rds_thule_2013_2014_wf2'));
  else
    fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2013_Greenland_P3/CSARP_insar/',sprintf('rds_thule_2013_2014_wf2'));
  end
  % Enable for waveform 3
%   rbins = [];
%   if ispc
%     fn = fullfile('X:/ct_data/rds/2013_Greenland_P3/CSARP_insar/',sprintf('rds_thule_2013_2014_wf3'));
%   else
%     fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2013_Greenland_P3/CSARP_insar/',sprintf('rds_thule_2013_2014_wf3'));
%   end
  equalization = [equalization_2014 equalization_2013];
  
  baseline_master_idx = 8;
  master_idx = 8;
  if 1
    output_fn_midfix = '';
    pass_en_mask = true(1,7+15);
  end

elseif 0
  %% 2011
  insar_mode = 2; % 1 to find equalization coefficients, 2 to process data, 3 to differential SAR

  wf = 1;
  if wf == 1
    rbins = [200:420];
    equalization = 10.^(zeros(1,15)/20) .* exp(1i*([-17.8 -14.4 -16.3 -15.2 -20.0 -19.4 -16.0 0.0 -10.4 -9.5 -12.2 -28.3 -31.6 -32.0 -19.4])/180*pi);
  end
  
  if ispc
    fn = fullfile('X:/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20110502_02_032_wf%d.mat',wf));
  else
    fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20110502_02_032_wf%d.mat',wf));
  end
%   if ispc
%     fn = fullfile('X:/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20110506_01_004_wf%d.mat',wf));
%   else
%     fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20110506_01_004_wf%d.mat',wf));
%   end
%   if ispc
%     fn = fullfile('X:/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20110509_01_004_wf%d.mat',wf));
%   else
%     fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20110509_01_004_wf%d.mat',wf));
%   end
%   if ispc
%     fn = fullfile('X:/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20110509_02_034_wf%d.mat',wf));
%   else
%     fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20110509_02_034_wf%d.mat',wf));
%   end

  baseline_master_idx = 16;
  master_idx = 8;
  if 0
    output_fn_midfix = '';
    pass_en_mask = true(1,16);
    pass_en_mask(16) = false; % Disable the 2014 pass
  elseif 1
    baseline_master_idx = 8;
    output_fn_midfix = '_ref2011';
    pass_en_mask = true(1,16);
    pass_en_mask(16) = false; % Disable the 2014 pass
  elseif 0
    output_fn_midfix = '_master';
    pass_en_mask = true(1,16);
    pass_en_mask([1:7 9:15]) = false; % Disable all but the master indices
  else
    output_fn_midfix = '_center';
    pass_en_mask = true(1,16);
    pass_en_mask(8:16) = false; % Disable all but the center elements
  end
  
elseif 0
  %% 2012
  wf = 1;
  insar_mode = 1; % 1 to find equalization coefficients, 2 to process data, 3 to differential SAR

  if wf == 1
%     rbins = [200:420];
%     equalization = ones(1,15);
    equalization = 10.^(zeros(1,15)/20) .* exp(1i*([-8.0 -7.3 0.0 -10.8 -0.2 -8.4 1.5 -22.7 -26.3 -21.8 -13.9 11.6 12.5 14.1 14.0])/180*pi);
    % -11.6 -10.6 0.0 -14.0 0.5 -13.0 -4.2 -33.6 -36.4 -27.0 -23.3 6.6 5.9 6.9 7.2
  end
  
%   if ispc
%     fn = fullfile('X:/ct_data/rds/2012_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20120503_03_067_wf%d.mat',wf));
%   else
%     fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2012_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20120503_03_067_wf%d.mat',wf));
%   end
  if ispc
    fn = fullfile('X:/ct_data/rds/2012_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20120516_01_089_wf%d.mat',wf));
  else
    fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2012_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20120516_01_089_wf%d.mat',wf));
  end

  master_idx = 3;
  % For combined file

elseif 0
  %% 2013
  wf = 1;
  insar_mode = 2; % 1 to find equalization coefficients, 2 to process data, 3 to differential SAR, 4 to check coregistration

  if wf == 1
    rbins = [200:420];
%     equalization = ones(1,7);
    equalization = 10.^(zeros(1,7)/20) .* exp(1i*([-3.9 -5.0 0.0 -8.1 -0.5 -4.3 1.1])/180*pi);
    coregistration_time_shift = [1.43 1.43 1.43 1.43 1.43 1.43 1.43 0];
%     coregistration_time_shift = zeros(1,8);
  end
  
%   if ispc
%     fn = fullfile('X:/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20130419_01_004_wf%d.mat',wf));
%   else
%     fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20130419_01_004_wf%d.mat',wf));
%   end
  if ispc
    fn = fullfile('X:/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20130426_01_004_wf%d.mat',wf));
  else
    fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20130426_01_004_wf%d.mat',wf));
  end

  baseline_master_idx = 8;
  master_idx = 3;
  if 1
    output_fn_midfix = '';
    pass_en_mask = true(1,8);
    pass_en_mask(8) = false; % Disable the 2014 pass
  elseif 0
    baseline_master_idx = 3;
    output_fn_midfix = '_ref2013';
    pass_en_mask = true(1,8);
    pass_en_mask(8) = false; % Disable the 2014 pass
  elseif 0
    output_fn_midfix = '_master';
    pass_en_mask = true(1,8);
    pass_en_mask([1:2 4:7]) = false; % Disable all but the master indices
  end
  
elseif 0
  %% 2014
  wf = 3;
  insar_mode = 2; % 1 to find equalization coefficients, 2 to process data, 3 to differential SAR

  if wf == 1
    rbins = [130:220];
    equalization = 10.^(zeros(1,15)/20) .* exp(1i*([122.5 121.6 126.7 105.5 130.1 120.1 128.7 -0.0 -134.2 121.1 36.4 125.8 -171.5 -1.0 128.4]-[42.8 45.5 48.2 50.9 53.6 56.3 58.9 0.0 2.7 5.4 8.1 93.6 96.3 99.0 101.7]/2)/180*pi);
  elseif wf == 2
    rbins = 280:420;
    equalization = 10.^(zeros(1,15)/20) .* exp(1i*([122.5 121.6 126.7 105.5 130.1 120.1 128.7 -0.0 -134.2 121.1 36.4 125.8 -171.5 -1.0 128.4]-[42.8 45.5 48.2 50.9 53.6 56.3 58.9 0.0 2.7 5.4 8.1 93.6 96.3 99.0 101.7]/2)/180*pi);
  elseif wf == 3
    equalization = 10.^(zeros(1,15)/20) .* exp(1i*([122.5 121.6 126.7 105.5 130.1 120.1 128.7 -0.0 -134.2 121.1 36.4 125.8 -171.5 -1.0 128.4]-[42.8 45.5 48.2 50.9 53.6 56.3 58.9 0.0 2.7 5.4 8.1 93.6 96.3 99.0 101.7]/2)/180*pi);
    rbins = 420:500;
    rbins = [];
  end
%   if ispc
%     fn = fullfile('X:/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20140429_01_005_wf%d.mat',wf));
%   else
%     fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20140429_01_005_wf%d.mat',wf));
%   end
  if ispc
    fn = fullfile('X:/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20140429_01_067_wf%d.mat',wf));
  else
    fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20140429_01_067_wf%d.mat',wf));
  end
%   if ispc
%     fn = fullfile('X:/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_combine_wf%d.mat',wf));
%     equalization = [equalization equalization equalization equalization];
%   else
%     fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_combine_wf%d.mat',wf));
%     equalization = [equalization equalization equalization equalization];
%   end
  master_idx = 8;

elseif 0
  insar_mode = 2; % 1 to find equalization coefficients, 2 to process data, 3 to differential SAR
  if ispc
    fn = 'X:/ct_data/rds/2016_Greenland_TOdtu/CSARP_insar/north.mat';
  else
    fn = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Greenland_TOdtu/CSARP_insar/north.mat';
  end
  master_idx = 1;
  rbins = 20:300;
elseif 0
  insar_mode = 2; % 1 to find equalization coefficients, 2 to process data, 3 to differential SAR
  if ispc
    fn = 'X:/ct_data/rds/2016_Greenland_TOdtu/CSARP_insar/middle.mat';
  else
    fn = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Greenland_TOdtu/CSARP_insar/middle.mat';
  end
  master_idx = 10;
  rbins = 20:300;
elseif 0
  insar_mode = 2; % 1 to find equalization coefficients, 2 to process data, 3 to differential SAR
  if ispc
    fn = 'X:/ct_data/rds/2016_Greenland_TOdtu/CSARP_insar/south.mat';
  else
    fn = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Greenland_TOdtu/CSARP_insar/south.mat';
  end
  master_idx = 1;
  rbins = 20:300;
elseif 0
  insar_mode = 2; % 1 to find equalization coefficients, 2 to process data, 3 to differential SAR
  if ispc
    fn = 'X:/ct_data/rds/2016_Greenland_TOdtu/CSARP_insar/iceland_south.mat';
  else
    fn = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Greenland_TOdtu/CSARP_insar/iceland_south.mat';
  end
  master_idx = 1;
  rbins = 20:300;
elseif 0
  insar_mode = 2; % 1 to find equalization coefficients, 2 to process data, 3 to differential SAR
  if ispc
    fn = 'X:/ct_data/rds/2016_Greenland_TOdtu/CSARP_insar/iceland_north.mat';
  else
    fn = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Greenland_TOdtu/CSARP_insar/iceland_north.mat';
  end
  master_idx = 1;
  rbins = 20:300;
elseif 0
  insar_mode = 2; % 1 to find equalization coefficients, 2 to process data, 3 to differential SAR
  if ispc
    fn = 'X:/ct_data/rds/2016_Greenland_G1XB/CSARP_insar/good_line.mat';
  else
    fn = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Greenland_G1XB/CSARP_insar/good_line.mat';
  end
  master_idx = 1;
  rbins = 20:100;
elseif 0
  insar_mode = 2; % 1 to find equalization coefficients, 2 to process data, 3 to differential SAR
  if ispc
    fn = 'X:/ct_data/rds/2016_Greenland_G1XB/CSARP_insar/medium_line.mat';
  else
    fn = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Greenland_G1XB/CSARP_insar/medium_line.mat';
  end
  master_idx = 2;
  rbins = 20:100;
elseif 0
  insar_mode = 2; % 1 to find equalization coefficients, 2 to process data, 3 to differential SAR
  if ispc
    fn = 'X:/ct_data/rds/2016_Greenland_G1XB/CSARP_insar/bad_line.mat';
  else
    fn = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Greenland_G1XB/CSARP_insar/bad_line.mat';
  end
  master_idx = 2;
  rbins = 20:100;
end

%% Setup
% =========================================================================

proj = geotiffinfo(ct_filename_gis([],fullfile('greenland','Landsat-7','Greenland_natural_90m.tif')));

physical_constants;

load(fn);

% Input check
% -----------------------
if ~exist('pass_en_mask','var')
  pass_en_mask = [];
end
% Assume any undefined passes are enabled
pass_en_mask(end+1:length(pass)) = true;
pass_en_idxs = find(pass_en_mask);

if isempty(rbins)
  rbins = 1:size(pass(pass_idx).data,1);
end
for pass_idx = 1:length(pass)
  if pass_en_mask(pass_idx)
    rbins = intersect(rbins,1:size(pass(pass_idx).data,1));
  end
end

if ~exist('output_fn_midfix','var') || isempty(output_fn_midfix)
  output_fn_midfix = '';
end

% For passes with undefined equalization, set to 1
if ~exist('equalization','var') || isempty(equalization)
  equalization = ones(1,length(pass));
elseif numel(equalization) < length(pass)
  equalization(end+1:length(pass)) = 1;
end

if ~exist('coregistration_time_shift','var') || isempty(coregistration_time_shift)
  coregistration_time_shift = zeros(1,length(pass));;
end

%% Plot Results
% =========================================================================
h_fig_map = figure(100); clf;
h_plot_map = [];
h_legend_map = {};
hold on;
axis('equal');
h_fig_elev = figure(101); clf;
h_plot_elev = [];
h_legend_elev = {};
hold on;
xlabel('Range line');
ylabel('WGS-84 elevation (m)');
grid on;

h_data_axes = [];
for pass_idx = 1:length(pass)
  if pass_en_mask(pass_idx)
    figure(pass_idx); clf;
    set(pass_idx,'WindowStyle','docked')
    imagesc(lp(pass(pass_idx).data(rbins,:)))
    colormap(1-gray(256));
    h_data_axes(end+1) = gca;
  end
  
  % Apply GPS time offset
  if 0
    pass_idx = 5;
    time_offset = -5;
    pass(pass_idx).ecef = interp1(pass(pass_idx).gps_time,pass(pass_idx).ecef.',pass(pass_idx).gps_time+time_offset,'linear','extrap').';
    pass(pass_idx).x = interp1(pass(pass_idx).gps_time,pass(pass_idx).x.',pass(pass_idx).gps_time+time_offset,'linear','extrap').';
  end
  
  pass(pass_idx).ecef = pass(pass_idx).origin;
  for rline = 1:size(pass(pass_idx).origin,2)
    pass(pass_idx).ecef(:,rline) = pass(pass_idx).ecef(:,rline) ...
      + [pass(pass_idx).x(:,rline) pass(pass_idx).y(:,rline) pass(pass_idx).z(:,rline)]*pass(pass_idx).pos(:,rline);
  end
  [pass(pass_idx).lat,pass(pass_idx).lon,pass(pass_idx).elev] = ecef2geodetic(referenceEllipsoid('wgs84'), ...
    pass(pass_idx).ecef(1,:), pass(pass_idx).ecef(2,:), pass(pass_idx).ecef(3,:));
  
  figure(h_fig_map);
  if 0
    h_plot = plot(pass(pass_idx).lon, pass(pass_idx).lat,'.');
    color = get(h_plot,'Color');
    h_text = text(pass(pass_idx).lon(1), pass(pass_idx).lat(1), sprintf('%d', pass_idx), 'Color', color);
  else
    [pass(pass_idx).proj_x,pass(pass_idx).proj_y] = projfwd(proj,pass(pass_idx).lat,pass(pass_idx).lon);
    h_plot_map(end+1) = plot(pass(pass_idx).proj_x/1e3, pass(pass_idx).proj_y/1e3,'.');
    h_legend_map{end+1} = sprintf('%d',pass_idx);
    color = get(h_plot_map(end),'Color');
    h_text = text(pass(pass_idx).proj_x(1)/1e3, pass(pass_idx).proj_y(1)/1e3, sprintf('%d', pass_idx), 'Color', color);
    xlabel('X (km)');
    ylabel('Y (km)');
    grid on;
  end
  
  figure(h_fig_elev);
  if baseline_master_idx == pass_idx
    h_plot_elev(end+1) = plot(pass(pass_idx).elev,'LineWidth',2,'UserData',pass_idx);
  else
    h_plot_elev(end+1) = plot(pass(pass_idx).elev,'UserData',pass_idx);
  end
  h_legend_elev{end+1} = sprintf('%d',pass_idx');
end
linkaxes(h_data_axes,'xy');
legend(h_plot_map,h_legend_map);

%% Co-Register Results
% =========================================================================


if 1
  % Option 1: Use a single pass as the reference.
  ref = pass(baseline_master_idx);

else
  % Option 2: Take a single pass as master and construct along track vectors relative
  % to it. Fit a polynomial to all the data using the along track vector
  % from step 1. Use this as the reference. This may be useful to do if no
  % single pass follows the center of the tube of passes.
  
  % With reference, run SAR coord system in a special mode where the mean is
  % not taken across Lsar, but instead each point is directly passed to the
  % output.
end

along_track = geodetic_to_along_track(ref.lat,ref.lon,ref.elev);
ref.surface_bin = interp1(ref.wfs(ref.wf).time, 1:length(ref.wfs(ref.wf).time), ref.surface);

if 0
  h_fig_ref_idx = figure(102); clf;
  hold on;
end

data = [];
for pass_idx = 1:length(pass)
    
  pass(pass_idx).ref_idx = zeros(1,size(pass(pass_idx).origin,2));
  last_idx = 0;
  for rline = 1:size(pass(pass_idx).ecef,2)
    % offset: offset from position rline of pass pass_idx from every point
    %         on the ref line
    offset = bsxfun(@minus, pass(pass_idx).ecef(:,rline), ref.ecef);
    % dist: converts offset into along-track distance of the slave line
    dist = offset.'*pass(pass_idx).x(:,rline);
    % min_idx: finds the point on the reference line that is closest to
    %          this point
    [min_dist,min_idx] = min(abs(dist));
    %       if min_idx == last_idx
    %         keyboard
    %       end
    last_idx = min_idx;
    pass(pass_idx).ref_idx(rline) = min_idx;
    
    % x_offset: along-track offset on master line of the closest point
    x_offset = offset(:,min_idx).'*ref.x(:,min_idx);
    % Compute FCS of slave point in master line coordinate system
    %   FCS: flight (aka SAR) coordinate system
    pass(pass_idx).along_track(rline) = along_track(min_idx) + x_offset;
    pass(pass_idx).ref_y(rline) = offset(:,min_idx).'*ref.y(:,min_idx);
    pass(pass_idx).ref_z(rline) = offset(:,min_idx).'*ref.z(:,min_idx);

    if 0
      % Compute the location of all pixels from this range line in ECEF
      pass(pass_idx).wfs(pass(pass_idx).wf).time;
      time = pass(pass_idx).time(2)-pass(pass_idx).time(1);
      
      range = time * c/2;
      range(time>pass(pass_idx).surface(rline)) = pass(pass_idx).surface(rline)*c/2 ...
        + (time(time>pass(pass_idx).surface(rline)) - pass(pass_idx).surface(rline))*c/2/sqrt(er_ice);
      
      pixels = pass(pass_idx).ecef(:,rline) + pass(pass_idx).z(:,rline)*range;
      
      % Compute the location of all pixels from this range line in the master
      % line coordinate system FCS.
      Tfcs = [ref.x, ref.y, ref.z];
      pass(pass_idx).pixels(:,:,rline) = (pixels - ref.origin) / Tfcs;
    end
  end
  
  pass(pass_idx).along_track_slave = geodetic_to_along_track(pass(pass_idx).lat,pass(pass_idx).lon,pass(pass_idx).elev);
  
  if 0
    % Debug plot showing indexes for alignment of passes
    figure(h_fig_ref_idx);
    plot(pass(pass_idx).ref_idx)
    drawnow;
  end
  
  % Resample images and position vectors onto a common along-track axes
  % 1. Oversample slave data by 10x in along track
  Mx = 10;
  Nx = size(pass(pass_idx).data,2);
  data_oversample = interpft(pass(pass_idx).data.',Mx*Nx);
  % 2. Interpolate to find the oversampled slave axes
  along_track_oversample = interp1(0:Nx-1, ...
    pass(pass_idx).along_track, (0:Nx*Mx-1)/Mx,'linear','extrap');
  % 3. Interpolate oversampled slave data onto master along track axes
  pass(pass_idx).ref_data = interp1(along_track_oversample, ...
    data_oversample, along_track,'linear','extrap').';
  
  pass(pass_idx).ref_y = interp1(pass(pass_idx).along_track, ...
    pass(pass_idx).ref_y, along_track,'linear','extrap').';
  pass(pass_idx).ref_z = interp1(pass(pass_idx).along_track, ...
    pass(pass_idx).ref_z, along_track,'linear','extrap').';

  % Apply fixed coregistration time shift
  freq = pass(pass_idx).wfs(pass(pass_idx).wf).freq;
  freq = freq - freq(1); % Remove center frequency offset
  dt = coregistration_time_shift(pass_idx) * (pass(pass_idx).wfs(pass(pass_idx).wf).time(2)-pass(pass_idx).wfs(pass(pass_idx).wf).time(1));
  pass(pass_idx).ref_data = ifft(bsxfun(@times,fft(pass(pass_idx).ref_data),exp(-1i*2*pi*freq*dt)));
  
  if insar_mode == 1
    % Motion compensation of FCS z-motion
    for rline = 1:size(pass(pass_idx).ref_data,2)
      % Convert z-offset into time-offset assuming nadir DOA
      dt = pass(pass_idx).ref_z(rline)/(c/2);
      pass(pass_idx).ref_data(:,rline) = ifft(fft(pass(pass_idx).ref_data(:,rline)) ...
        .*exp(1i*2*pi*pass(pass_idx).wfs(pass(pass_idx).wf).freq*dt) );
    end
  elseif insar_mode == 2 || insar_mode == 4
    % Co-register images using GPS and nadir squint angle assumption
    %
    % Motion compensation of FCS z-motion without center frequency so there
    % is no phase shift.
    freq = pass(pass_idx).wfs(pass(pass_idx).wf).freq;
    freq = freq - freq(1); % Remove center frequency offset
    for rline = 1:size(pass(pass_idx).ref_data,2)
      % Convert z-offset into time-offset assuming nadir DOA
      dt = pass(pass_idx).ref_z(rline)/(c/2);
      pass(pass_idx).ref_data(:,rline) = ifft(fft(pass(pass_idx).ref_data(:,rline)) ...
        .*exp(1i*2*pi*freq*dt) );
    end
  elseif insar_mode == 3
    % Motion compensation of FCS z-motion and slope compensation

    % True time delay shift for z-offset
    for rline = 1:size(pass(pass_idx).ref_data,2)
      % Convert z-offset into time-offset assuming nadir DOA
      dt = pass(pass_idx).ref_z(rline)/(c/2);
      pass(pass_idx).ref_data(:,rline) = ifft(fft(pass(pass_idx).ref_data(:,rline)) ...
        .*exp(1i*2*pi*pass(pass_idx).wfs(pass(pass_idx).wf).freq*dt) );
    end

    % Phase only correction for slope
    if 1
      % Using file generated from this dataset
      [fn_dir,fn_name] = fileparts(fn);
      fn_slope = fullfile(fn_dir,[fn_name '_slope.mat']);
      load(fn_slope,'slope','GPS_time','Latitude','Longitude','Elevation','Time','Surface');
      slope = interp1(GPS_time,slope.',pass(baseline_master_idx).gps_time).';
      slope = interp_finite(slope.').';
      slope = interp1(Time,slope,pass(pass_idx).wfs(pass(pass_idx).wf).time);
      slope = interp_finite(slope);
      
      pass(pass_idx).ref_data = pass(pass_idx).ref_data .* exp(-1i*4*pi*pass(pass_idx).wfs(pass(pass_idx).wf).freq(1)/c *bsxfun(@times,sin(slope),pass(pass_idx).ref_y(:).'));
      
    elseif 1
      % Using file generated from another dataset
      fn = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/rds_thule_20140429_01_067_wf2_slope.mat';
      
      % TBD
      
    end
    
  elseif 0
    % Co-register images using cross-correlation
    keyboard
  end

  % Match time axis to baseline_master_idx
  if 0
    pass(pass_idx).ref_data = interp1(pass(pass_idx).wfs(pass(pass_idx).wf).time, pass(pass_idx).ref_data, pass(baseline_master_idx).wfs(pass(baseline_master_idx).wf).time, 'linear', 0);
  else
    Mt = 4;
    Nt = length(pass(pass_idx).wfs(pass(pass_idx).wf).time);
    pass(pass_idx).ref_data = interpft(pass(pass_idx).ref_data,Mt*Nt);
    dt = pass(pass_idx).wfs(pass(pass_idx).wf).time(2)-pass(pass_idx).wfs(pass(pass_idx).wf).time(1);
    time_Mt = pass(pass_idx).wfs(pass(pass_idx).wf).time(1) + dt/Mt*(0:Mt*Nt-1);
    pass(pass_idx).ref_data = interp1(time_Mt, pass(pass_idx).ref_data, pass(baseline_master_idx).wfs(pass(baseline_master_idx).wf).time, 'linear', 0);
  end
  
  if 0
    % Normalize surface phase
    Nt = size(pass(pass_idx).ref_data,1);
    Nx = size(pass(pass_idx).ref_data,2);
    H = pass(baseline_master_idx).ref_data(round(ref.surface_bin)+(0:Nx-1)*Nt) .* conj(pass(pass_idx).ref_data(round(ref.surface_bin)+(0:Nx-1)*Nt));
    H = exp(1i*angle(H));
    pass(pass_idx).ref_data = bsxfun(@times,pass(pass_idx).ref_data,H);
  end
  
  % Concatenate data into a single matrix
  data = cat(3,data,pass(pass_idx).ref_data);
end

%% Apply equalization
% -----------------------
if insar_mode == 2 || insar_mode == 3 || insar_mode == 4
  equalization = reshape(equalization,[1 1 numel(equalization)]);
  data(:,:,pass_en_idxs) = bsxfun(@times,data(:,:,pass_en_idxs),1./equalization(:,:,pass_en_idxs));
end

if 0
  %% Coregister: Data Dependent method to estimate System Time Delay
  % Apply fixed coregistration time shift
  for pass_out_idx = 2%1:length(pass_en_idxs)
    pass_idx = pass_en_idxs(pass_out_idx);
    freq = pass(pass_idx).wfs(pass(pass_idx).wf).freq;
    freq = freq - freq(1); % Remove center frequency offset
    coregistration_time_shifts = -2:0.05:2;
    coregistration_time_shifts = -1.6:0.01:-1.3;
%     coregistration_time_shifts = -0.2:0.01:0.2;
    coherence_sum = [];
    for coregistration_time_shift_idx = 1:length(coregistration_time_shifts)
      coregistration_time_shift = coregistration_time_shifts(coregistration_time_shift_idx);
      dt = coregistration_time_shift * (pass(pass_idx).wfs(pass(pass_idx).wf).time(2)-pass(pass_idx).wfs(pass(pass_idx).wf).time(1));
      adjusted = ifft(bsxfun(@times,fft(data(:,:,pass_idx)),exp(-1i*2*pi*freq*dt)));
      coherence = fir_dec(adjusted(rbins,:) .* conj(data(rbins,:,master_idx)) ./ abs(adjusted(rbins,:) .* data(rbins,:,master_idx)),ones(1,7)/7,1);
      coherence = fir_dec(coherence.',ones(1,3)/3,1).';
      coherence = abs(coherence);
      %     coherence_sum(coregistration_time_shift_idx) = sum(coherence(coherence>0.5));
      coherence_sum(coregistration_time_shift_idx) = sum(coherence(coherence>0));
      %    TriangleRayIntersection fprintf('%g %.2f\n', coregistration_time_shift, coherence_sum(coregistration_time_shift_idx));
      %     imagesc(coherence); colormap(1-gray(256));
      %     pause
    end
  end
  figure(1000); clf;
  plot(coregistration_time_shifts,coherence_sum)
  [~,coregistration_time_shift_idx] = max(coherence_sum);
  coregistration_time_shift = coregistration_time_shifts(coregistration_time_shift_idx)
  return
end

%% Plot interferograms
h_data_axes = [];
new_equalization = [];
for pass_out_idx = 1:length(pass_en_idxs)
  pass_idx = pass_en_idxs(pass_out_idx);
  
  figure(pass_idx); clf;
  set(pass_idx,'WindowStyle','docked')
  if 0
    imagesc(lp(data(rbins,:,pass_idx)))
    colormap(1-gray(256));
    ylabel('Range bin');
    xlabel('Range line');
    title(sprintf('%s_%03d %d',pass(pass_idx).param_sar.day_seg,pass(pass_idx).param_sar.load.frm,pass(pass_idx).direction),'interpreter','none')
    caxis([-90 8]);
  else
    % Form interferogram (couple options)
    complex_data = fir_dec(data(rbins,:,pass_idx) .* conj(data(rbins,:,master_idx)),ones(1,11)/11,1);
    if ~exist('equalization_rlines','var') || isempty(equalization_rlines)
      new_equalization(pass_idx) = mean(complex_data(:)); % equalization only valid when motion compensation with phase is used
    else
      new_equalization(pass_idx) = mean(mean(complex_data(:,equalization_rlines))); % equalization only valid when motion compensation with phase is used
    end
    if insar_mode == 1
      complex_data = complex_data ./ new_equalization(pass_idx);
    end
    % Plot interferogram
    if insar_mode == 4
      imagesc(lp(data(rbins,:,pass_idx)));
      colormap(1-gray(256));
      h_colorbar = colorbar;
      set(get(h_colorbar,'ylabel'),'string','Relative power (dB)');
    else
      coherence = abs(fir_dec(data(rbins,:,pass_idx) .* conj(data(rbins,:,master_idx)) ./ abs(data(rbins,:,pass_idx) .* data(rbins,:,master_idx)),ones(1,11)/11,1)) ...
        .* exp(1i*angle(complex_data));
      imagesc(hsv_plot_coherence(coherence,[0 1]));
      colormap(hsv(256))
      h_colorbar = colorbar;
      caxis([-pi pi])
      set(get(h_colorbar,'ylabel'),'string','Angle (radians)');
    end
    ylabel('Range bin');
    xlabel('Range line');
    title(sprintf('%s_%03d %d',pass(pass_idx).param_sar.day_seg,pass(pass_idx).param_sar.load.frm,pass(pass_idx).direction),'interpreter','none')
    
    if 0
      imagesc(abs(coherence));
      colormap(1-gray(256))
      h_colorbar = colorbar;
      caxis([0 1])
      set(get(h_colorbar,'ylabel'),'string','Coherence');
      ylabel('Range bin');
      xlabel('Range line');
      title(sprintf('%s_%03d %d',pass(pass_idx).param_sar.day_seg,pass(pass_idx).param_sar.load.frm,pass(pass_idx).direction),'interpreter','none')
    end
    
  end
  h_data_axes(end+1) = gca;
  
end
fprintf('=============================================\n');
fprintf('New equalization\n');
fprintf('%.1f ', lp(new_equalization)-mean(lp(new_equalization(pass_en_idxs))));
fprintf('\n');
fprintf('%.1f ', angle(new_equalization)*180/pi)
fprintf('\n');
fprintf('=============================================\n');
linkaxes(h_data_axes,'xy');
if insar_mode == 1 || insar_mode == 4
  return
end
if insar_mode == 3
  [fn_dir,fn_name] = fileparts(fn);
  fn_insar = fullfile(fn_dir,[fn_name '_insar.mat']);
  param_sar = pass(master_idx).param_sar;
  param_records = pass(master_idx).param_records;
  save(fn_insar,'-v7.3','data','ref','param_sar','param_records');
  return
end

%% Array Processing

% Package data to call array_proc.m
% 1. Data
% 2. Trajectory and attitude
% 3. Array processing parameters
data = {permute(data,[1 2 4 5 3])};

param.array = [];
param.array.method = 1;
% param.array.Nsv = 256; param.array = rmfield(param.array,'theta'); % Forces default theta
% param.array.theta = linspace(-20,20,256);
param.array.theta = linspace(-6,6,256);
param.array.Nsrc = 2;
param.array.bin_rng = [-4:4];
param.array.line_rng = [-20:20];
param.array.dbin = 1;
param.array.dline = 11;
param.array.freq_rng = 1;
h_fig_baseline = figure(200); clf;
h_plot_baseline = [];
h_legend_baseline = {};
for pass_out_idx = 1:length(pass_en_idxs)
  pass_idx = pass_en_idxs(pass_out_idx);
  
  param.array.fcs{1}{pass_out_idx}.pos = along_track;
  param.array.fcs{1}{pass_out_idx}.pos(2,:) = pass(pass_idx).ref_y;
  param.array.fcs{1}{pass_out_idx}.pos(3,:) = pass(pass_idx).ref_z;
  param.array.fcs{1}{pass_out_idx}.base_line ...        
    = sqrt( (pass(pass_idx).ref_z - pass(master_idx).ref_z).^2 ...
      + (pass(pass_idx).ref_y - pass(master_idx).ref_y).^2 );
    
  h_plot_baseline(end+1) = plot(param.array.fcs{1}{pass_out_idx}.base_line);
  h_legend_baseline{end+1} = sprintf('%d',pass_idx);
  hold on;

  param.array.fcs{1}{pass_out_idx}.surface = ref.surface;
end
xlabel('Range line');
ylabel('Baseline (m)');
grid on;
legend(h_plot_baseline,h_legend_baseline);

param.array.wfs.time = ref.wfs(ref.wf).time;
dt = param.array.wfs.time(2)-param.array.wfs.time(1);
param.array_proc.bin0 = param.array.wfs.time/dt;
param.array.sv_fh = @array_proc_sv;
param.array.wfs.fc = ref.wfs(ref.wf).fc;
param.array.imgs = {[ones(length(pass_en_idxs),1), (1:length(pass_en_idxs)).']};
param.array.tomo_en = true;

%%
array_proc_methods;
param = array_proc(param);
param.array.method = STANDARD_METHOD;
[param_array0,result0] = array_proc(param,{data{1}(:,:,:,:,pass_en_idxs)});
% param.array.method = MVDR_METHOD;
% [param_array1,result1] = array_proc(param,{data{1}(:,:,:,:,pass_en_idxs)});
param.array.method = MUSIC_METHOD;
[param_array2,result2] = array_proc(param,{data{1}(:,:,:,:,pass_en_idxs)});

figure(101); clf;
imagesc(lp(result0.img))
title('Periodogram')
colormap(1-gray(256))
h_axes = gca;

% figure(102); clf;
% imagesc(lp(result1.img))
% title('MVDR')
% colormap(1-gray(256))
% h_axes(end+1) = gca;

figure(103); clf;
imagesc(lp(result2.img))
title('MUSIC')
colormap(1-gray(256))
h_axes(end+1) = gca;

figure(104); clf;
imagesc(lp(fir_dec(abs(data{1}(:,:,master_idx)).^2, ones(size(param.array.line_rng)), param.array.dline, ...
  1-param.array.line_rng(1), size(data{1}(:,:,master_idx),2)-length(param.array.line_rng)+1)))
title('Single Channel');
colormap(1-gray(256))
h_axes(end+1) = gca;

linkaxes(h_axes,'xy');

%% Save Results
[fn_dir,fn_name] = fileparts(fn);

Tomo = result0.tomo;
Data = result0.img;
GPS_time = ref.gps_time(param_array0.array_proc.lines);
Latitude = ref.lat(param_array0.array_proc.lines);
Longitude = ref.lon(param_array0.array_proc.lines);
Elevation = ref.elev(param_array0.array_proc.lines);
Roll = ref.roll(param_array0.array_proc.lines);
Pitch = ref.pitch(param_array0.array_proc.lines);
Heading = ref.heading(param_array0.array_proc.lines);
Surface = ref.surface(param_array0.array_proc.lines);
Bottom = nan(size(Surface));
param_sar = pass(master_idx).param_sar;
param_records = pass(master_idx).param_records;
param_array = param_array0;
Time = param.array.wfs.time(param_array0.array_proc.bins);
file_version = '1';
fn_mat = fullfile(fn_dir,[fn_name output_fn_midfix '_standard.mat']);
save('-v7.3',fn_mat,'Tomo','Data','Latitude','Longitude','Elevation','GPS_time', ...
  'Surface','Bottom','Time','param_array','param_records', ...
  'param_sar', 'Roll', 'Pitch', 'Heading', 'file_version');

Tomo = result2.tomo;
Data = result2.img;
GPS_time = ref.gps_time(param_array2.array_proc.lines);
Latitude = ref.lat(param_array2.array_proc.lines);
Longitude = ref.lon(param_array2.array_proc.lines);
Elevation = ref.elev(param_array2.array_proc.lines);
Roll = ref.roll(param_array2.array_proc.lines);
Pitch = ref.pitch(param_array2.array_proc.lines);
Heading = ref.heading(param_array2.array_proc.lines);
Surface = ref.surface(param_array2.array_proc.lines);
Bottom = nan(size(Surface));
param_sar = pass(master_idx).param_sar;
param_records = pass(master_idx).param_records;
param_array = param_array2;
Time = param.array.wfs.time(param_array2.array_proc.bins);
file_version = '1';
fn_mat = fullfile(fn_dir,[fn_name output_fn_midfix '_music.mat']);
save('-v7.3',fn_mat,'Tomo','Data','Latitude','Longitude','Elevation','GPS_time', ...
  'Surface','Bottom','Time','param_array','param_records', ...
  'param_sar', 'Roll', 'Pitch', 'Heading', 'file_version');

fn_map = fullfile(fn_dir,[fn_name '_map.fig']);
saveas(h_fig_map,fn_map);
fn_elev = fullfile(fn_dir,[fn_name '_elev.fig']);
saveas(h_fig_elev,fn_elev);
fn_baseline = fullfile(fn_dir,[fn_name '_baseline.fig']);
saveas(h_fig_baseline,fn_baseline);
fn_standard = fullfile(fn_dir,[fn_name '_standard.fig']);
saveas(101,fn_standard);
% fn_mvdr = fullfile(fn_dir,[fn_name '_mvdr.fig']);
% saveas(102,fn_mvdr);
fn_music = fullfile(fn_dir,[fn_name '_music.fig']);
saveas(103,fn_music);
fn_single = fullfile(fn_dir,[fn_name '_single.fig']);
saveas(104,fn_single);

return



%% Other plots and plot setup


figure(500); clf;
for idx = 1800:-1:600%1:50:size(result1.img,3)
  imagesc(lp(result0.img(20:120,:,idx)))
  colormap(jet(256));
  caxis([-80 -20])
  idx
  pause
end

%%

figure; plot(array_param2.theta*180/pi,lp(result1.img(20+54-1,:,1291)))
grid on
xlim([-90 90])
xlabel('Direction of arrival (deg)');
ylabel('Relative power (dB)');
title('Bad line, MVDR passes [2 4 6], rline 1291');

figure(h_fig_ref_idx); h_axes = gca;
%figure(h_fig_map); h_axes(end+1) = gca;
figure(h_fig_elev); h_axes(end+1) = gca;
linkaxes(h_axes,'x');

h_fig_combined = figure(103); clf;
imagesc(lp(mean(data(rbins,:,[1 6]),3)));
colormap(1-gray(256));
set(h_fig_combined,'WindowStyle','Docked');
hold on;
plot(ref.surface_bin-rbins(1)+1);

surf_data = zeros(length(pass),size(data,2));
for pass_idx = 1:length(pass)
  for rline = 1:size(data,2)
    surf_data(pass_idx,rline) = data(round(ref.surface_bin(rline)),rline,pass_idx);
  end
end

h_fig_angle = figure(104); clf;
complex_data = surf_data(6,:) .* conj(surf_data(1,:));
complex_data(lp(complex_data) < -96) = NaN;
plot(180/pi*angle(complex_data),'.')
xlabel('Range line');
ylabel('Angle (deg)');
grid on;
