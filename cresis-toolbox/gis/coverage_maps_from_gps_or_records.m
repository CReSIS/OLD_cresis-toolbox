% script coverage_maps_from_gps_or_records
%
% Example:
%  coverage_maps_from_gps_or_records % after modifying user settings
%  saveas(fig_h,'my_map.jpg')
%  % Add some other graphics
%  fid = fopen(ct_filename_gis([],'antarctica\GroundingLine\InSAR_Ant_GroundingLine\InSAR_GL_Antarctica.txt'));
%  GL = textscan(fid,'%f%f%*s%*s%*s%*s%*s%*s%*s%*s%*s');
%  fclose(fid);
%  [x,y] = projfwd(proj,GL{1},GL{2});
%  hold on;
%  plot(x/1e3,y/1e3,'k.');
%
% Author: John Paden

%% User Settings

% fig_h = figure handle to use for plotting
fig_h = 1;

% dt = time spacing (sec) between each plotted point
dt = 1;

% data_source: string 'gps' or 'records' (which files in csarp_support to use)
data_source = 'gps'; 

% Greenland
% geotiff_fn = '/cresis/snfs1/dataproducts/GIS_data/arctic/NaturalEarth_Data/Arctic_NaturalEarth.tif';
geotiff_fn = 'X:\GIS_data/arctic/NaturalEarth_Data/Arctic_NaturalEarth.tif';

% Antarctica
geotiff_fn = ct_filename_gis(gRadar,'antarctica/NaturalEarth_Data/Antarctica_NaturalEarth.tif');
axis_limits = [-2800    500       -1500       1800];

% Cell vector of param spreadsheet filenames
param_fns = {ct_filename_param('snow_param_2014_Greenland_P3.xls')};

% plot_color: string argument to plot function (e.g. 'b.' for blue dots)
plot_color = 'b.';

%% Automated Section
figure(fig_h); clf;
proj = plot_geotiff(geotiff_fn,[],[],fig_h);
axis normal; axis equal;
if ~isempty(axis_limits)
  axis(axis_limits);
end
xlabel('X (km)');
ylabel('Y (km)');

for param_fn_idx = 1:length(param_fns)
  param_fn = param_fns{param_fn_idx};
  fprintf('Parameters %s\n', param_fn);
  
  params = read_param_xls(param_fn);
  
  yyyymmdd_list = {};
  for param_idx = 1:length(params)
    param = params(param_idx);
    
    param.yyyymmdd = param.day_seg(1:8);
    if strcmp(data_source,'gps') && any(strcmp(param.yyyymmdd,yyyymmdd_list))
      % Already done this segment
      continue;
    end
    
    if param.cmd.generic && isempty(regexpi(param.cmd.notes,'do not process'))
      yyyymmdd_list{end+1} = param.yyyymmdd;
      if strcmp(data_source,'gps')
        gps_fn = ct_filename_support(param,[],'gps',1);
        fprintf('  Processing %s: %s (%s)\n', param.day_seg, gps_fn, datestr(now,'HH:MM:SS'));
        gps = load(gps_fn,'gps_time','lat','lon');
      elseif strcmp(data_source,'records')
        records_fn = ct_filename_support(param,[],'records');
        fprintf('  Processing %s: %s (%s)\n', param.day_seg, records_fn, datestr(now,'HH:MM:SS'));
        gps = load(records_fn,'gps_time','lat','lon');
      end
      
      idxs = logical(zeros(size(gps.gps_time)));
      idxs(1) = 1;
      cur_time = gps.gps_time(1);
      for idx = 2:length(gps.gps_time)
        if gps.gps_time(idx) > cur_time + dt
          idxs(idx) = true;
          cur_time = gps.gps_time(idx);
        end
      end
      [x,y] = projfwd(proj,gps.lat(idxs),gps.lon(idxs));
      hold on;
      h = plot(x/1e3,y/1e3,plot_color,'UserData',param.day_seg);
      drawnow;
      hold off;
    end
    
  end
  
end
