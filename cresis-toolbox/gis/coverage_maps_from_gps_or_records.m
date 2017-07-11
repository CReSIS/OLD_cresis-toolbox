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

global gRadar
%% User Settings

% fig_h = figure handle to use for plotting
fig_h = 1;

% dt = time spacing (sec) between each plotted point
dt = 1;

% data_source: string 'gps' or 'records' (which files in csarp_support to use)
data_source = 'records';

% Greenland
% geotiff_fn = ct_filename_gis(gRadar,fullfile('arctic','NaturalEarth_Data','Arctic_NaturalEarth.tif'));
% geotiff_fn = ct_filename_gis(gRadar,fullfile('arctic','Landsat-7','arctic_natural_90m.tif'));
geotiff_fn = ct_filename_gis(gRadar,fullfile('greenland','Landsat-7','Greenland_natural.tif'));
% geotiff_fn = ct_filename_gis(gRadar,fullfile('canada','Landsat-7','Canada_90m.tif'));
axis_limits = [-1024        -207       -1363        -560];

% Antarctica
% geotiff_fn = ct_filename_gis(gRadar,fullfile('antarctica','NaturalEarth_Data','Antarctica_NaturalEarth.tif'));
% axis_limits = [-2800    500       -1500       1800];

% Cell vector of param spreadsheet filenames
param_fns = {ct_filename_param('rds_param_2014_Greenland_P3.xls')};

% plot_args: cell array of argument to plot function (e.g. 'b.' for blue dots)
plot_args = {'b-','LineWidth',1};

% along_track_sampling: desired along track sampling of plot (m)
along_track_sampling = 50;

% label: Add text labels to the plots for specific frames
%  .day_seg: cell list of day segments to match
%  .frm: equal length cell vector to .day_seg with specific frames to label
%  .text: equal length cell vector to .day_seg with string for text label
label.day_seg = {};
label.frm = {};
label.text = {};
% label.day_seg = {'20091102_02','20091016_01','20091016_01','20091102_02','20091102_02'};
% label.frm = {8, 21, 26, 23, 32};
% label.text = {' 1',' 2',' 3',' 4',' 5'};

% figure_position: Set figure Position property (leave empty to use default/current figure size)
figure_position = [];
% figure_position = [221   -60   654   704];

%% Automated Section
figure(fig_h); clf;
if exist(figure_position,'var') && ~isempty(figure_position)
  set(fig_h,'Position',figure_position);
end
proj = plot_geotiff(geotiff_fn,[],[],fig_h);
axis normal; axis equal;
if ~isempty(axis_limits)
  axis(axis_limits);
end
xlabel('X (km)');
ylabel('Y (km)');

label.h_text = [];
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
    
    if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
      continue;
    end
    
    if isempty(regexpi(param.cmd.notes,'do not process'))
      yyyymmdd_list{end+1} = param.yyyymmdd;
      if strcmp(data_source,'gps')
        gps_fn = ct_filename_support(param,[],'gps',1);
        fprintf('  Processing %s: %s (%s)\n', param.day_seg, gps_fn, datestr(now,'HH:MM:SS'));
        gps = load(gps_fn,'gps_time','lat','lon');
      elseif strcmp(data_source,'records')
        records_fn = ct_filename_support(param,[],'records');
        if ~exist(records_fn,'file')
          warning('No records file: %s\n', records_fn);
          continue;
        end
        frames_fn = ct_filename_support(param,param.records.frames_fn,'frames');
        if ~exist(frames_fn,'file')
          warning('No frames file: %s\n', frames_fn);
          continue;
        end
        
        fprintf('  Processing %s: %s (%s)\n', param.day_seg, records_fn, datestr(now,'HH:MM:SS'));
        gps = load(records_fn,'gps_time','lat','lon');
        
        if ~isfield(param.records,'frames_fn')
          param.records.frames_fn = '';
        end
        
        % Load frames file
        load(frames_fn);
        
        if isempty(param.cmd.frms)
          param.cmd.frms = 1:length(frames.frame_idxs);
        end
        % Remove frames that do not exist from param.cmd.frms list
        [valid_frms,keep_idxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
        if length(valid_frms) ~= length(param.cmd.frms)
          bad_mask = ones(size(param.cmd.frms));
          bad_mask(keep_idxs) = 0;
          warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
            param.cmd.frms(find(bad_mask,1)));
          param.cmd.frms = valid_frms;
        end
        % Force into a row vector
        param.cmd.frms = param.cmd.frms(:).';
      end

      if strcmp(data_source,'gps')
        % Force GPS time to be monotonic (only necessary for corrupt GPS data)
        idxs = logical(zeros(size(gps.gps_time)));
        idxs(1) = 1;
        cur_time = gps.gps_time(1);
        for idx = 2:length(gps.gps_time)
          if gps.gps_time(idx) > cur_time + dt
            idxs(idx) = true;
            cur_time = gps.gps_time(idx);
          end
        end
        
        % Plot
        [x,y] = projfwd(proj,gps.lat(idxs),gps.lon(idxs));
        hold on;
        h = plot(x/1e3,y/1e3,plot_args{:},'UserData',param.day_seg);
        drawnow;
        hold off;
        
      elseif strcmp(data_source,'records')
        lat = NaN*zeros(size(gps.lat));
        for frm = param.cmd.frms
          first_idx = frames.frame_idxs(frm);
          if frm ~= length(frames.frame_idxs)
            last_idx = frames.frame_idxs(frm+1)-1;
          else
            last_idx = length(gps.gps_time);
          end
          lat(first_idx:last_idx) = gps.lat(first_idx:last_idx);
        end
        
        along_track = geodetic_to_along_track(gps.lat,gps.lon);
        decim_idxs = get_equal_alongtrack_spacing_idxs(along_track,along_track_sampling);
        
        [x,y] = projfwd(proj,lat(decim_idxs),gps.lon(decim_idxs));
        hold on;
        h = plot(x/1e3,y/1e3,plot_args{:},'UserData',param.day_seg);
        drawnow;
        
        if isfield('label','day_seg')
          match_idxs = find(strcmpi(param.day_seg,label.day_seg(:).'));
          for match_idx = match_idxs
            for frm = label.frm{match_idx}
              first_idx = frames.frame_idxs(frm);
              [x,y] = projfwd(proj,gps.lat(first_idx),gps.lon(first_idx));
              plot(x/1e3,y/1e3,'.','MarkerSize',12,'Color','Black');
              label.h_text(match_idx) = text(x/1e3,y/1e3,label.text{match_idx},'Color','black','FontWeight','bold','FontSize',14);
            end
          end
        end
        
      end
    end
    
  end
  
end

return

%% Below are common steps to customize the figure
axis normal;
axis_limits = [-3190        -428       -1284        1478];
axis(axis_limits);
axis equal;
axis([-2562        -950       -1615        1893]);

pos = get(label.h_text(2),'Position');
set(label.h_text(2),'Position',pos+[-50 75 0]);

pos = get(label.h_text(3),'Position');
set(label.h_text(3),'Position',pos-[10 40 0]);

pos = get(label.h_text(5),'Position');
set(label.h_text(5),'Position',pos-[-15 -25 0]);

h_children = get(gca,'Children');
h_img = h_children(end);
clip_and_resample_image(h_img,gca,10);


