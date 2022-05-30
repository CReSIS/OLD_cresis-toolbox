%% User Settings
% =========================================================================

param = [];
vel_fn = {};
vel_mult = [];

if ispc
  vel_fn_dir = 'Y:\cbarnett\year_greenland_vv\VX_vel_components';
  % FOR X-COMPONENT VELOCITIES USE: 
  % 'Y:\cbarnett\year_greenland_vv\VX_vel_components'
else
  vel_fn_dir = '/cresis/snfs1/scratch/cbarnett/year_greenland_vv/';
end

example_str = 'Petermann_line1_2011_2014_2015_2017_2018_2019_X';

if strcmpi(example_str, 'Petermann_line1_2011_2014_2015_2017_2018_2019_X')
  %% Petermann Line 1 - 2011, 2014, 2015, 2017, 2018, 2019
  if ispc
    fn = fullfile('X:\ct_data\rds\2014_Greenland_P3\CSARP_multipass\',sprintf('Petermann_line1_2011_2014_2015_2017_2018_2019_multipass02.mat'));
  else
    fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_multipass/',sprintf('Petermann_line1_2011_2014_2015_2017_2018_2019_multipass02'));
  end
  vel_fn{1}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2012_2013_vx_v02.1.tif');
  vel_mult{1} = [3];
  vel_fn{2}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2014_2015_vx_v02.1.tif');
  vel_mult{2} = [1];
  vel_fn{3}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2015_2016_vx_v02.1.tif');
  vel_fn{3}{2} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2016_2017_vx_v02.1.tif');
  vel_mult{3} = [1 1];
  vel_fn{4}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2017_2018_vx_v02.1.tif');
  vel_mult{4} = [1];
  % velocity correction for 2018-2019 season uses 2017-2018 data since no
  % NCSIDC velocity map exists for this period of time
  vel_fn{5}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2017_2018_vx_v02.1.tif'); 
  vel_mult{5} = [1];
end

if strcmpi(example_str, 'Petermann_line2_2007_2013_2014_2017')
  %% Petermann Line 2 - 2007, 2013, 2014, 2017
  if ispc
    fn = fullfile('X:\ct_data\rds\2014_Greenland_P3\CSARP_multipass\',sprintf('Petermann_line2_2007_2013_2014_2017_multipass.mat'));
  else
    fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_multipass/',sprintf('Petermann_line2_2007_2013_2014_2017_multipass'));
  end
  % velocity correction for 2007-2013 uses 2009-2010 and 2012-2013
  % velocities twice to correct for missing data between 2010-2011 and
  % 2011-2012
%   vel_fn{1}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2007_2008_vv_v02.1.tif');
%   vel_fn{1}{2} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2008_2009_vv_v02.1.tif');
%   vel_fn{1}{3} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2009_2010_vv_v02.1.tif');
%   vel_fn{1}{4} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2012_2013_vv_v02.1.tif');
%   vel_mult{1} = [1 1 2 2];
  vel_fn{2}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2013_2014_vv_v02.1.tif');
  vel_mult{2} = [1];
  vel_fn{3}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2014_2015_vv_v02.1.tif');
  vel_fn{3}{2} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2015_2016_vv_v02.1.tif');
  vel_fn{3}{3} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2016_2017_vv_v02.1.tif');
  vel_mult{3} = [1 1 1];
end

% if strcmpi(example_str, 'Petermann_line3_2007_2010A_2010B_2017')
%   % Petermann Line 3 - 2007, 2010A, 2010B, 2017 
%   if ispc
%     fn = fullfile('X:\ct_data\rds\2014_Greenland_P3\CSARP_multipass\',sprintf('Petermann_line3_2007_2010_2010DC8_2017_multipass.mat'));
%   else
%     fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_multipass/',sprintf('Petermann_line3_2007_2010_2010CD8_2017_multipass'));
%   end
% end

if strcmpi(example_str, 'Petermann_line4_2010_2011_2013_2014_2017_X')
  %% Petermann Line 4 - 2010, 2011, 2013, 2014, 2017 x COMPONENT
  if ispc
    fn = fullfile('X:\ct_data\rds\2014_Greenland_P3\CSARP_multipass\',sprintf('Petermann_line4_2010_2011_2013_2014_2017_multipass02.mat'));
  else
    fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_multipass/',sprintf('Petermann_line4_2010_2011_2013_2014_2017_multipass02'));
  end
  % 2010-2011 and 2011-2013 use NSIDC years 2012-2013 correction since none
  % exists for the covered time period.
  vel_fn{1}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2009_2010_vx_v02.1.tif');
  vel_mult{1} = [1];
  vel_fn{2}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2009_2010_vx_v02.1.tif');
  vel_fn{2}{2} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2012_2013_vx_v02.1.tif');
  vel_mult{2} = [1 1];
  vel_fn{3}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2012_2013_vx_v02.1.tif');
  vel_mult{3} = [3];
  vel_fn{4}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2014_2015_vx_v02.1.tif');
  vel_fn{4}{2} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2015_2016_vx_v02.1.tif');
  vel_fn{4}{3} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2016_2017_vx_v02.1.tif');
  vel_mult{4} = [1 1 1];
end

if strcmpi(example_str, 'Petermann_oblique_2014_2017_2018')
  %% Petermann Oblique 2014, 2017, 2018 
  % NOT WORKING CURRENTLY, CANNOT FIND multipass files
  if ispc
    fn = fullfile('X:\ct_data\rds\2014_Greenland_P3\CSARP_multipass\',sprintf('Petermann_oblique_2014_2017_2018_multipass02.mat'));
  else
    fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_multipass/',sprintf('Petermann_oblique_2014_2017_2018_multipass02'));
  end
  vel_fn{1}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2014_2015_vv_v02.1.tif');
  vel_mult{1} = [1];
  vel_fn{2}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2015_2016_vv_v02.1.tif');
  vel_fn{2}{2} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2016_2017_vv_v02.1.tif');
  vel_mult{2} = [1 1];
  vel_fn{3}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2017_2018_vv_v02.1.tif');
  vel_mult{3} = [1];
end

if strcmpi(example_str, '79N_line1_2010_2014_2016_2018_2019')
  %% 79N Line 1 - 2010, 2014, 2016, 2018, 2019
  if ispc
    fn = fullfile('X:\ct_data\rds\2014_Greenland_P3\CSARP_multipass\',sprintf('79N_line1_2010_2014_2016_2018_2019_multipass.mat'));
  else
    fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_multipass/',sprintf('79N_line1_2010_2014_2016_2018_2019_multipass'));
  end
  
  vel_fn{1}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2012_2013_vv_v02.1.tif');
  vel_fn{1}{2} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2013_2014_vv_v02.1.tif');
  vel_mult{1} = [3 1];
  vel_fn{2}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2014_2015_vv_v02.1.tif');
  vel_fn{2}{2} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2015_2016_vv_v02.1.tif');
  vel_mult{2} = [1 1];
  vel_fn{3}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2016_2017_vv_v02.1.tif');
  vel_fn{3}{2} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2017_2018_vv_v02.1.tif');
  vel_mult{3} = [1 1];
  vel_fn{4}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2017_2018_vv_v02.1.tif');
  vel_mult{4} = [1];
end

if strcmpi(example_str, 'Humboldt_line1_2012_2013_2014_2017')
    %% Humboldt Glacier - 2012, 2013, 2014, 2017
    if ispc
      fn = fullfile('X:\ct_data\rds\2013_Greenland_P3\CSARP_multipass\',sprintf('Humboldt_line1_2012_2013_2014_2017_multipass.mat'));
    else
      fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2013_Greenland_P3/CSARP_multipass/',sprintf('Humboldt_line1_2012_2013_2014_2017_multipass'));
    end
    vel_fn{1}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2012_2013_vv_v02.1.tif');
    vel_mult{1} = [1];
    vel_fn{2}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2013_2014_vv_v02.1.tif');
    vel_mult{2} = [1];
    vel_fn{3}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2014_2015_vv_v02.1.tif');
    vel_fn{3}{2} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2015_2016_vv_v02.1.tif');
    vel_fn{3}{3} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2016_2017_vv_v02.1.tif');
    vel_mult{3} = [1 1 1];
end

if strcmpi(example_str, 'Ryder_line1_2011_2013_2015_2019')
    %% Ryder Glacier - 2011, 2013, 2015, 2019
    if ispc
      fn = fullfile('X:\ct_data\rds\2013_Greenland_P3\CSARP_multipass\',sprintf('Ryder_line1_2011_2013_2015_2019_multipass02.mat'));
    else
      fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2013_Greenland_P3/CSARP_multipass/',sprintf('Ryder_line1_2011_2013_2015_2019_multipass02'));
    end
    vel_fn{1}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2012_2013_vv_v02.1.tif');
    vel_mult{1} = [2];
    vel_fn{2}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2014_2015_vv_v02.1.tif');
    vel_mult{2} = [2];
    vel_fn{3}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2015_2016_vv_v02.1.tif');
    vel_fn{3}{2} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2016_2017_vv_v02.1.tif');
    vel_fn{3}{3} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2017_2018_vv_v02.1.tif');
    vel_mult{3} = [1 1 2];
    % Use the NSIDC 2017-2018 velocity for the 2018-2019 season shift
end

if strcmpi(example_str, 'Steensby_line1_2011_2013_2015_2019')
    %% Steensby Glacier - 2011, 2013, 2015, 2019
    if ispc
      fn = fullfile('X:\ct_data\rds\2013_Greenland_P3\CSARP_multipass\',sprintf('Steensby_line1_2011_2013_2015_2019_multipass02.mat'));
    else
      fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2013_Greenland_P3/CSARP_multipass/',sprintf('Steensby_line1_2011_2013_2015_2019_multipass02'));
    end
        vel_fn{1}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2012_2013_vv_v02.1.tif');
    vel_mult{1} = [2];
    vel_fn{2}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2014_2015_vv_v02.1.tif');
    vel_mult{2} = [2];
    vel_fn{3}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2015_2016_vv_v02.1.tif');
    vel_fn{3}{2} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2016_2017_vv_v02.1.tif');
    vel_fn{3}{3} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2017_2018_vv_v02.1.tif');
    vel_mult{3} = [1 1 2];
end

if strcmpi(example_str, 'ZI_line1_2010_2014_2016_2017_2018_2019')
    %% Zachariae Isstrom - 2010, 2014_A, 2016, 2017, 2018, 2019 DC8 TEST
    if ispc
      fn = fullfile('X:\ct_data\rds\2014_Greenland_P3\CSARP_multipass\',sprintf('ZI_line1_2010_2014_2016_2017_2018_2019_multipass02.mat'));
    else
      fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_multipass/',sprintf('ZI_line1_2010_2014_2016_2017_2018_2019_multipass02'));
    end
    vel_fn{1}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2012_2013_vv_v02.1.tif');
    vel_mult{1} = [4];
    vel_fn{2}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic500_2014_2015_vv_v02.1.tif');
    vel_fn{2}{2} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2015_2016_vv_v02.1.tif');
    vel_mult{2} = [1 1];
    vel_fn{3}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2016_2017_vv_v02.1.tif');
    vel_mult{3} = [1];
    vel_fn{4}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2017_2018_vv_v02.1.tif');
    vel_mult{4} = [1];
    vel_fn{5}{1} = fullfile(vel_fn_dir,'greenland_vel_mosaic200_2017_2018_vv_v02.1.tif');
    vel_mult{5} = [1];
    % Use the NSIDC 2017-2018 velocity for the 2018-2019 season shift
end

%%
load(fn);

% Load multipass.multipass output
load(fn,'pass','param_multipass');

vel_param = []; vel_param.velocity_coregister = [];
if ~isfield(vel_param.velocity_coregister,'debug_plots') || isempty(vel_param.velocity_coregister.debug_plots)
  vel_param.velocity_coregister.debug_plots = {'visible'};
end
enable_visible_plot = any(strcmp('visible',vel_param.velocity_coregister.debug_plots));
if ~isempty(vel_param.velocity_coregister.debug_plots)
  h_fig = get_figures(4,enable_visible_plot);
end

%% Interpolate velocity
% =========================================================================
%  Interpolate the surface velocity maps for each pass. We do not
%  interpolate the last pass in the list because this function assumes it
%  to be the last pass in time and therefore does not need any velocity
%  correction. All other passes should have a list of velocity maps that
%  allow the velocity to be integrated over time to correct for velocity
%  from when the particular pass was taken to the last pass.

figure(h_fig(1)); clf;
baseline_master_idx = param_multipass.multipass.baseline_master_idx;
for pass_idx = 1:length(pass)
  
  if pass_idx == length(pass)
    pass(pass_idx).vel = zeros(1,length(pass(baseline_master_idx).along_track));
    break;
  end
  
  for sub_pass_idx = pass_idx:length(pass)-1
    for vel_idx = 1:length(vel_fn{sub_pass_idx})
      vel_proj = geotiffinfo(vel_fn{sub_pass_idx}{vel_idx});
      
      % Read the image
      [vel, R, tmp] = geotiffread(vel_fn{sub_pass_idx}{vel_idx});
      vel = vel*vel_mult{sub_pass_idx}(vel_idx);
      
      xaxis = R(3,1) + R(2,1)*(1:size(vel,2));
      yaxis = R(3,2) + R(1,2)*(1:size(vel,1));
      
      [pass(pass_idx).x,pass(pass_idx).y] = projfwd(vel_proj,pass(baseline_master_idx).lat,pass(baseline_master_idx).lon);
      
      % CHANGE THE pass(pass_idx).vel TO vel_x OR vel_y FOR X AND Y VELS
      if vel_idx == 1 && sub_pass_idx == pass_idx
        pass(pass_idx).vel = interp2(xaxis,yaxis,vel,pass(pass_idx).x,pass(pass_idx).y);
        pass(pass_idx).velx = pass(pass_idx).vel;
      else
        pass(pass_idx).vel = pass(pass_idx).vel + interp2(xaxis,yaxis,vel,pass(pass_idx).x,pass(pass_idx).y);
        pass(pass_idx).velx = pass(pass_idx).vel; 
      end
      
      % Debug: plot velocity geotiff
      if 0
        figure(h_fig(1));
        imagesc(xaxis/1e3, yaxis/1e3, vel);
        set(gca,'YDir','normal');
        colorbar;
        caxis([0 1500*vel_mult{sub_pass_idx}(vel_idx)]);
        axis([357         532       -1244        -975]);
        % axis([-299        -144       -1140        -933]);
        % caxis([0 5000]);
        hold on
        plot(pass(pass_idx).x/1e3,pass(pass_idx).y/1e3,'k');
      end
    end
  end
end

%% Plot velocity corrected layers
% =========================================================================
%  Plot the integrated velocity and the surface velocity corrected ice surface
%  and bottom layers for each pass.
% figure(h_fig(2)); clf;
% figure(h_fig(3)); clf;
% figure(h_fig(4)); clf;
% leg_str = {};
% h_vel = [];
% h_plot = [];
% for pass_idx = 1:length(pass)
%   [year month day] = datevec(epoch_to_datenum(pass(pass_idx).gps_time(1)));
%   
%   figure(h_fig(2));
%   h_vel(pass_idx) = plot(pass(baseline_master_idx).along_track/1e3, pass(pass_idx).vel);
%   grid on;
%   hold on;
%   xlabel('Along-track (km)');
%   ylabel('Horizontal movement (m)');
%   
%   figure(h_fig(3));
%   lay_idx = 1;
%   alongtrack = pass(baseline_master_idx).along_track + pass(pass_idx).vel;
%   h_plot(pass_idx) = plot(alongtrack/1e3, pass(pass_idx).layers(lay_idx).layer_elev);
%   grid on;
%   leg_str{pass_idx} = num2str(year);
%   hold on;
%   for lay_idx = 2:length(pass(pass_idx).layers)
%     h_plot2 = plot(alongtrack/1e3, pass(pass_idx).layers(lay_idx).layer_elev);
%     set(h_plot2,'Color',get(h_plot(pass_idx),'Color'));
%   end
%   xlabel('Along-track (km)');
%   ylabel('Elevation (m,WGS-84)');
%   
%   figure(h_fig(4));
%   lay_idx = 1;
%   alongtrack = pass(baseline_master_idx).along_track;
%   h_plot(pass_idx) = plot(alongtrack/1e3, pass(pass_idx).layers(lay_idx).layer_elev);
%   grid on;
%   hold on;
%   for lay_idx = 2:length(pass(pass_idx).layers)
%     h_plot2 = plot(alongtrack/1e3, pass(pass_idx).layers(lay_idx).layer_elev);
%     set(h_plot2,'Color',get(h_plot(pass_idx),'Color'));
%   end
%   xlabel('Along-track (km)');
%   ylabel('Elevation (m,WGS-84)');
%   
% end
% legend(h_plot,leg_str,'location','best');
% legend(h_vel,leg_str,'location','best');
% 
% [fn_dir,fn_name] = fileparts(fn);
% fn_motion = fullfile(fn_dir,[fn_name '_motion.fig']);
% fn_layer_comp = fullfile(fn_dir,[fn_name '_layer_comp.fig']);
% fn_layer = fullfile(fn_dir,[fn_name '_layer.fig']);
% fprintf('Saving %s (%s)\n', fn_motion, datestr(now));
% saveas(2,fn_motion);
% fprintf('Saving %s (%s)\n', fn_layer_comp, datestr(now));
% saveas(3,fn_layer_comp);
% fprintf('Saving %s (%s)\n', fn_layer, datestr(now));
% saveas(4,fn_layer);
% 
