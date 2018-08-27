% script make_coverage_maps.m
% 
% Author: Rohan Choudhari

function test_make_coverage_maps(params, user_variables)

% This function creates the coverage map depending on the user_variables and
% params provided.
% 
%   user_variables: = struct whose fields are used to set the user settings
%     .data_source: string
%       'gps': Loads the data from gps files
%       'records': Loads the data from records files
%       'layer': Loads the layerData using opsLoadLayers
%     .proj: projection information returned by plot_geotiff.m
%       - eg. geotiff_fn = ct_filename_gis(gRadar,fullfile('arctic','NaturalEarth_Data','Arctic_NaturalEarth.tif'));
%             user_variables.proj = plot_geotiff(geotiff_fn,[],[], fig_h_idx);
%     .season_idx: integer
%       - index of the current season from the 
%     .save_img: 1 or 0
%       1: to save
%       2: to not save
%     .ext: string
%       '.m': MATLAB FIG-file and MATLAB code that opens figure 
%       '.jpg': JPEG image
%       '.fig': MATLAB® FIG-file
%     '.axis_limits': vector
%       - a vector of axis limits for the geotiff
%     '.along_track_sampling': integer
%       - desired along track sampling of plot (m)
%       - default - 100
%     '.plot_args': cell array
%       - cell array of argument to plot function (e.g. 'b.' for blue dots)
%       - default - {'LineWidth',1}
%     label: Add text labels to the plots for specific frames
%       .day_seg: cell list of day segments to match
%       .frm: equal length cell vector to .day_seg with specific frames to label
%       .text: equal length cell vector to .day_seg with string for text label
%     '.figure_position': vector
%       - to set the figure position
%       eg. [221   -60   654   704]


global gRadar;
yyyymmdd_list = {};


%Set the title for the figure, initialize required variables
if strcmp(user_variables.data_source,'gps') || strcmp(user_variables.data_source,'records')
  title(strcat(params(1).season_name,'_coverage_map'),'Interpreter','none');
  fprintf('  Creating coverage map for %s: \n', params(1).season_name);
elseif strcmp(user_variables.data_source,'layer')
  fprintf('  Creating coverage map for %s: \n', user_variables.season_name);
  
  % for plotting the quality of data
  bad_data_handles = [];
  moderate_data_handles = [];
  good_data_handles = [];
  
  bottom = [];
  quality = [];
  x = [];
  y = [];
end

%Set custom axis limits if provided by the user
if isfield(user_variables, 'axis_limits')
  if ~isempty(user_variables.axis_limits)
    axis(user_variables.axis_limits);
  end
end

%Set the axes labels
xlabel('X (km)');
ylabel('Y (km)');
user_variables.label.h_text = [];

%Looping through all the params for the season
if strcmp(user_variables.data_source,'gps') || strcmp(user_variables.data_source,'records')
  for param_idx = 1:length(params)
    
    parameter = params(param_idx);
    
    %Skipping segments that are already covered
    parameter.yyyymmdd = parameter.day_seg(1:8);
    if strcmp(user_variables.data_source,'gps') && any(strcmp(parameter.yyyymmdd,yyyymmdd_list))
      % Already done this segment
      return;
    end
    
    %Plotting if not 'do not process'
    if isempty(regexpi(parameter.cmd.notes,'do not process'))
      yyyymmdd_list{end+1} = parameter.yyyymmdd; %Adding segment to the list of segments already covered
      
      %Load lat, lon depending on data source
      if strcmp(user_variables.data_source,'gps')
        gps_fn = ct_filename_support(parameter,[],'gps',1);
        gps = load(gps_fn,'gps_time','lat','lon');
      elseif strcmp(user_variables.data_source,'records')
        records_fn = ct_filename_support(parameter,[],'records');
        if ~exist(records_fn,'file')
          warning('No records file: %s\n', records_fn);
          return;
        end
        frames_fn = ct_filename_support(parameter,'','frames');
        if ~exist(frames_fn,'file')
          warning('No frames file: %s\n', frames_fn);
          return;
        end
        gps = load(records_fn,'gps_time','lat','lon');
        
        % Load frames file
        load(frames_fn);
        
        if isempty(parameter.cmd.frms)
          parameter.cmd.frms = 1:length(frames.frame_idxs);
        end
        
        % Remove frames that do not exist from param.cmd.frms list
        [valid_frms,keep_idxs] = intersect(parameter.cmd.frms, 1:length(frames.frame_idxs));
        if length(valid_frms) ~= length(parameter.cmd.frms)
          bad_mask = ones(size(parameter.cmd.frms));
          bad_mask(keep_idxs) = 0;
          warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
            parameter.cmd.frms(find(bad_mask,1)));
          parameter.cmd.frms = valid_frms;
        end
        % Force into a row vector
        parameter.cmd.frms = parameter.cmd.frms(:).';
      end
      
      if strcmp(user_variables.data_source,'gps')
        % Force GPS time to be monotonic (only necessary for corrupt GPS data)
        idxs = logical(zeros(size(gps.gps_time)));
        idxs(1) = 1;
        cur_time = gps.gps_time(1);
        for idx = 2:length(gps.gps_time)
          if gps.gps_time(idx) > cur_time + user_variables.dt
            idxs(idx) = true;
            cur_time = gps.gps_time(idx);
          end
        end
        
        % Plot
        [x,y] = projfwd(user_variables.proj,gps.lat(idxs),gps.lon(idxs));
        hold on;
        h = plot(x/1e3,y/1e3,user_variables.plot_args{:},'UserData',parameter.day_seg);
        drawnow;
        hold off;
        
      elseif strcmp(user_variables.data_source,'records')
        lat = NaN*zeros(size(gps.lat));
        for frm = parameter.cmd.frms
          first_idx = frames.frame_idxs(frm);
          if frm ~= length(frames.frame_idxs)
            last_idx = frames.frame_idxs(frm+1)-1;
          else
            last_idx = length(gps.gps_time);
          end
          lat(first_idx:last_idx) = gps.lat(first_idx:last_idx);
        end
        
        along_track = geodetic_to_along_track(gps.lat,gps.lon);
        decim_idxs = get_equal_alongtrack_spacing_idxs(along_track,user_variables.along_track_sampling);        
        [x,y] = projfwd(user_variables.proj,lat(decim_idxs),gps.lon(decim_idxs));
        hold on;
        h = plot(x/1e3,y/1e3,user_variables.plot_args{:},'UserData',parameter.day_seg);
        drawnow;
        
        if isfield(user_variables.label,'day_seg')
          match_idxs = find(strcmpi(parameter.day_seg,user_variables.label.day_seg(:).'));
          for match_idx = match_idxs
            for frm = user_variables.label.frm{match_idx}
              first_idx = frames.frame_idxs(frm);
              [x,y] = projfwd(user_variables.proj,gps.lat(first_idx),gps.lon(first_idx));
              plot(x/1e3,y/1e3,'.','MarkerSize',12,'Color','Black');
              user_variables.label.h_text(match_idx) = text(x/1e3,y/1e3,user_variables.label.text{match_idx},'Color','black','FontWeight','bold','FontSize',14);
            end
          end
        end
        
      end
    end
    % saving the coverage map
    if user_variables.save_img == 1
      if strcmp(user_variables.data_source,'gps') || strcmp(user_variables.data_source,'records')
        image_fn_name = sprintf('%s_coverage_map_%s%s', params(1).season_name, date, user_variables.ext);
        image_fn = fullfile(user_variables.out_dir,image_fn_name);
        saveas(figure(user_variables.season_idx), image_fn);
      elseif strcmp(user_variables.data_source,'layer')
        image_fn = fullfile(user_variables.out_dir,'Master_coverage_map');
        saveas(figure(1), image_fn);
      end
    end
    
    
  end

elseif strcmp(user_variables.data_source,'layer')
  layer_params = user_variables.layer_params;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  x_red = [];
  x_magenta = [];
  x_yellow = [];
  x_green = [];
  y_red = [];
  y_yellow = [];
  y_green = [];
  y_magenta = [];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        
  for season_idx = 1:length(params)
    fprintf('Loading season %s \n',  user_variables.season_stats{season_idx}.details.season);
        
    points_before = 0;
    points_after = 0;
    user_variables.season_stats{season_idx}.seg = {};
    
    segs = {};
    for param_idx = 1:length(params{season_idx})
      param = params{season_idx}(param_idx);
      param = merge_structs(param,gRadar);
      seg = {};
      seg.name = param.day_seg;
      try
        layer = opsLoadLayers(param,layer_params);
      catch
%         fprintf('%s not found\n', param.day_seg);
        continue;
      end
      
%       fprintf('opsLoadLayers %s\n', param.day_seg);
      
      if isempty(layer)
        continue;
      end;
      for layer_idx = 1:length(layer)
        
        if length(layer(layer_idx).lat) ~= length(layer(layer_idx).type)
          fprintf('    Bad file\n');
          continue;
        end
        along_track = geodetic_to_along_track(layer(layer_idx).lat,layer(layer_idx).lon,zeros(size(layer(layer_idx).lat)));
        idxs = get_equal_alongtrack_spacing_idxs(along_track,500);
        points_before = points_before + length(along_track);
        points_after = points_after + length(idxs);
        
        for idx = 1:(length(idxs)-1)
          bottom_sec = isfinite(layer(layer_idx).twtt(idxs(idx):idxs(idx+1)-1));
          quality_sec = layer(layer_idx).quality(idxs(idx):idxs(idx+1)-1);
          [x_sec,y_sec] = projfwd(user_variables.proj,layer(layer_idx).lat(idxs(idx):idxs(idx+1)-1),layer(layer_idx).lon(idxs(idx):idxs(idx+1)-1));
          x_sec = x_sec/1e3;
          y_sec = y_sec/1e3;
          
          red = isnan(quality_sec);
          magenta = quality_sec==3;
          yellow = bottom_sec==1 & (quality_sec==1 | quality_sec==2);
          green = bottom_sec==1 & quality_sec==1;
          
          if(all(green==1))
            x_green = cat(2, x_green,x_sec(round(length(green)/2)));
            y_green = cat(2, y_green,y_sec(round(length(green)/2)));
          elseif(all(yellow==1))
            x_yellow = cat(2, x_yellow,x_sec(round(length(yellow)/2)));
            y_yellow = cat(2, y_yellow,y_sec(round(length(yellow)/2)));
          elseif(all(red==1))
            x_red = cat(2, x_red,x_sec(round(length(red)/2)));
            y_red = cat(2, y_red,y_sec(round(length(red)/2)));
          end
          if(any(magenta==1))
            x_magenta = cat(2, x_magenta,x_sec(round(length(magenta)/2)));
            y_magenta = cat(2, y_magenta,y_sec(round(length(magenta)/2)));
          end
          
        end
        
        
        
        %         bottom = cat(2,bottom,isfinite(layer(layer_idx).twtt(idxs)));
        %         quality = cat(2,quality,layer(layer_idx).quality(idxs));
        %         [x_tmp,y_tmp] = projfwd(user_variables.proj,layer(layer_idx).lat(idxs),layer(layer_idx).lon(idxs));
        %         x = cat(2,x,x_tmp/1e3);
        %         y = cat(2,y,y_tmp/1e3);
      end
      %       seg.before = length(along_track);
      %       seg.after = length(idxs);
      %       user_variables.season_stats{season_idx}.seg{end+1} = seg;
      %       fprintf('\t%s\t%s\t%d\t%d\n', user_variables.season_stats{season_idx}.details.season, param.day_seg, seg.before, seg.after);

    end
      plot(x_green,y_green,'g.');
      plot(x_yellow,y_yellow,'y.');
      plot(x_magenta,y_magenta,'m.');
      plot(x_red,y_red,'r.');

    
%     user_variables.season_stats{season_idx}.details.before = points_before;
%     user_variables.season_stats{season_idx}.details.after = points_after;
%     fprintf('\t%s\t%d\t%d\n', user_variables.season_stats{season_idx}.details.season, points_after, points_before);
     









    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AGGREGATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     no_data = bottom==0;
% %     no_data = isnan(quality);
% %     good_or_moderate = bottom==1 & (quality==1 | quality==2);
% %     good = bottom==1 & quality==1;
% %     
%     
%     length_quality = length(quality);
%     club = 15;
%     parts = fix(length_quality/club);
%     leftover = length_quality - (parts*club);
%     
%     dim(1:1, 1:parts) = club;
%     dim = cat(2, dim, leftover);
%     grouped = mat2cell(quality, 1, dim);
%     grouped_bottom = mat2cell(bottom, 1, dim);
%     grouped_x = mat2cell(x, 1, dim);
%     grouped_y = mat2cell(y, 1, dim);
%     
%     for row_idx = 1:length(grouped)
% %       no_data = grouped_bottom{row_idx}==0;
%       red = isnan(grouped{row_idx});
%       magenta = grouped{row_idx}==3;
%       yellow = grouped_bottom{row_idx}==1 & (grouped{row_idx}==1 | grouped{row_idx}==2);
%       green = grouped_bottom{row_idx}==1 & grouped{row_idx}==1;
%       
%       x_plot = grouped_x{row_idx};
%       y_plot = grouped_y{row_idx};
%       
%       if(all(green==1))
%         x_green = cat(2, x_green,x_plot(round(length(green)/2))); 
%         y_green = cat(2, y_green,y_plot(round(length(green)/2))); 
% %         plot(x_plot(round(length(green)/2)),y_plot(round(length(green)/2)),'g.');
%       elseif(all(yellow==1))        
%         x_yellow = cat(2, x_yellow,x_plot(round(length(yellow)/2))); 
%         y_yellow = cat(2, y_yellow,y_plot(round(length(yellow)/2))); 
% %         plot(x_plot(round(length(yellow)/2)),y_plot(round(length(yellow)/2)),'y.');
%       elseif(all(red==1))
%         x_red = cat(2, x_red,x_plot(round(length(red)/2))); 
%         y_red = cat(2, y_red,y_plot(round(length(red)/2)));
% %         plot(x_plot(round(length(red)/2)),y_plot(round(length(red)/2)),'r.');
%       end
%       if(any(magenta==1))
%         x_magenta = cat(2, x_magenta,x_plot(round(length(magenta)/2)));
%         y_magenta = cat(2, y_magenta,y_plot(round(length(magenta)/2)));
%       end
%     end
%     
%   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLAIN DECIMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     no_bottom_mask = bottom==0;
%     moderate_mask = bottom==1 & quality > 1;
%     bottom_mask = bottom==1 & ~moderate_mask;
%     
%     if any(no_bottom_mask)
%       bad_data_handles(season_idx) = plot(x(no_bottom_mask),y(no_bottom_mask),'r.');
%     else
%       bad_data_handles(season_idx) = plot(1,NaN,'r.');
%     end
%     if any(moderate_mask)
%       moderate_data_handles(season_idx) = plot(x(moderate_mask),y(moderate_mask),'y.');
%     else
%       moderate_data_handles(season_idx) = plot(1,NaN,'y.');
%     end
%     if any(bottom_mask)
%       good_data_handles(season_idx) = plot(x(bottom_mask),y(bottom_mask),'g.');
%     else
%       good_data_handles(season_idx) = plot(1,NaN,'g.');
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   plot(x_green,y_green,'g.');
%   plot(x_yellow,y_yellow,'y.');
%   plot(x_red,y_red,'r.');
%   plot(x_magenta,y_magenta,'m.');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % saving the coverage map
    if user_variables.save_img == 1
      for ext_idx = 1:length(user_variables.ext)
        if strcmp(user_variables.ext{ext_idx}, '.fig')
%           image_fn_name = sprintf('saveas_2014_P3_coverage_map_%s', date);
%           image_fn = fullfile(user_variables.out_dir,image_fn_name);
%           saveas(figure(1), image_fn);

          image_fn_name = sprintf('aggregate_greenland_master_coverage_map_%s', date);
          image_fn = fullfile(user_variables.out_dir,image_fn_name);
          savefig(figure(1), image_fn, 'compact');

%           image_fn_name = sprintf('hgsave_2014_P3_coverage_map_%s', date);
%           image_fn = fullfile(user_variables.out_dir,image_fn_name);
%           hgsave(1, strcat('',image_fn), '-v7.3');
        else
%           image_fn_name = sprintf('Master_fig_coverage_map_%s%s', date, user_variables.ext{ext_idx});
%           image_fn = fullfile(user_variables.out_dir,image_fn_name);
%           saveas(figure(1), image_fn);

        end
        %         saveas(image_fn, user_variables.ext{ext_idx},'-v7.3');
      end
    end
    
end
return

