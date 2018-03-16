Author: Mohanad Al-Ibadi 
% clear all;clc;startup

%Group 1:
% frm_inters = {{'20140401_03_009','20140401_03_024'},{'20140401_03_008','20140401_03_025'},{'20140401_03_010','20140401_03_017'}...
%   {'20140401_03_010','20140401_03_021'},{'20140401_03_028','20140401_03_028'},...
%   {'20140401_03_025','20140401_03_032'},{'20140401_03_026','20140401_03_031'}};

% Group 2:
% frm_inters = {{'20140401_03_007','20140401_03_033'},{'20140401_03_012','20140401_03_013'},{'20140506_01_007','20140506_01_008'},...
% {'20140506_01_004','20140506_01_008'},{'20140506_01_018','20140506_01_029'}};

% Group 3:
% frm_inters = {{'20140401_03_035','20140401_03_041'},{'20140506_01_037','20140506_01_040'},{'20140506_01_035','20140506_01_036'},...
%   {'20140506_01_030','20140506_01_030'},{'20140506_01_028','20140506_01_021'},...
%   {'20140506_01_028','20140506_01_024'},{'20140506_01_012','20140506_01_015'},{'20140506_01_007','20140506_01_009'},...
%   {'20140506_01_007','20140506_01_007'},{'20140506_01_005','20140506_01_007'},{'20140506_01_005','20140506_01_008'}};

% I didn't find any intersection between these frames even though they
% appear to be intersected from the picker.
% frm_inters = {{'20140401_03_035','20140401_03_042'},{'20140401_03_045','20140401_03_046'},...
%   {'20140401_03_041','20140401_03_038'},{'20140401_03_045','20140401_03_039'},{'20140401_03_012','20140401_03_015'}};

% Bad intersections:
% {{'20140506_01_004','20140506_01_007'},{'20140401_03_017','20140401_03_021'},{'20140506_01_034','20140506_01_034'},...
% {'20140401_03_028','20140401_03_028'},{'20140506_01_030','20140506_01_030'},{'20140506_01_007','20140506_01_007'}

% Couldn't handle Slice position problem:
% {{'20140506_01_039','20140506_01_041'}}}

% All groups (the good ones):
% frm_inters = {{'20140401_03_009','20140401_03_024'},{'20140401_03_008','20140401_03_025'},{'20140401_03_010','20140401_03_017'}...
%   {'20140401_03_010','20140401_03_021'},{'20140401_03_025','20140401_03_032'},{'20140401_03_026','20140401_03_031'},...
%   {'20140401_03_007','20140401_03_033'},{'20140401_03_012','20140401_03_013'},{'20140506_01_007','20140506_01_008'},...
%   {'20140506_01_004','20140506_01_008'},{'20140506_01_018','20140506_01_029'},{'20140506_01_037','20140506_01_040'},...
%   {'20140506_01_035','20140506_01_036'},{'20140506_01_028','20140506_01_021'},{'20140506_01_028','20140506_01_024'},...
%   {'20140506_01_012','20140506_01_015'},{'20140506_01_007','20140506_01_009'},{'20140506_01_005','20140506_01_007'},...
%   {'20140506_01_005','20140506_01_008'},{'20140401_03_035','20140401_03_041'}};

RadConf18_flag = 0;
% For RadConf18 ===. Remember to set doa_trim = 8 when you generate these DEMs
% frm_inters = {{'20140401_03_009','20140401_03_024'},{'20140401_03_010','20140401_03_017'},{'20140401_03_010','20140401_03_021'},...
%   {'20140401_03_007','20140401_03_033'},{'20140506_01_004','20140506_01_008'},{'20140506_01_018','20140506_01_029'},...
%   {'20140506_01_037','20140506_01_040'},{'20140506_01_035','20140506_01_036'},{'20140506_01_028','20140506_01_021'},...
%   {'20140506_01_028','20140506_01_024'},{'20140506_01_005','20140506_01_007'},{'20140506_01_005','20140506_01_008'}};

frm_inters = {{'20140506_01_004','20140506_01_007'}};

if RadConf18_flag
  % For RadConf18
  out_dir = '/users/mohanad/IceSheetProject/Conferences/Radar conference/figures/for the paper';
else
  % General output directory
  out_dir = '/users/mohanad/IceSheetProject/Conferences/Radar conference/figures';
end

param.radar_name   = 'rds';
param.season_name  = '2014_Greenland_P3';
surfdata_source    = 'surfData';
DEM_source         = 'CSA_DEM';
data_source        = 'music3D';

layer = 'bottom';
% layer = 'ice surface';
layers_plot = {'ice surface','bottom'};
save_flag     = 0;
save_all_flag = 0;
plot_flag     = 1;
colorbar_fig_flag = 0;

root_mean_squared_err_trim_tmp = 0;
failed_intersections_idx = 0;

for frm_inter_idx = 1:length(frm_inters)
  fprintf('\n Working on the intersection of frame %s and frame %s\n',frm_inters{frm_inter_idx}{1},frm_inters{frm_inter_idx}{2});
  day_seg_frm_cell = frm_inters{frm_inter_idx};
  day_seg_frm_cell_1 = day_seg_frm_cell(1);
  day_seg_frm_cell_2 = day_seg_frm_cell(2);
  param.day_seg = day_seg_frm_cell_1{1}(1:end-4);
  frms = [str2num(fliplr(day_seg_frm_cell_1{1}(end:-1:end-2))), str2num(fliplr(day_seg_frm_cell_2{1}(end:-1:end-2)))];
  
  %   param = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),param.day_seg,'post');
  for frm_idx = 1:length(frms)
    %% Call the DEMs (i.e. geotiff) and plot them separately
    frm = frms(frm_idx);
    surfdata_fn = fullfile(ct_filename_out(param,surfdata_source),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
    data_fn = fullfile(ct_filename_out(param,data_source),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
    load(surfdata_fn)
    surf_names = {surf.name};
    bot_idx = strmatch(layer,surf_names,'exact');
    layer_idx = bot_idx;
    
    surfdata_all{frm_idx} = surf;
    surfdata{frm_idx} = surf(layer_idx).y;
    data_tmp = load(data_fn);
    data{frm_idx} = data_tmp.Topography.img;
    
    
    surf_str = lower(surf(layer_idx).name);
    idx = regexp([' ' surf_str],'(?<=\s+)\S','start')-1;
    surf_str(idx) = upper(surf_str(idx));
    if ~strncmp(surf_str,'Ice',3)
      surf_str = ['Ice ',surf_str];
    end
    
    surf_str_lower = lower(surf_str);
    layer_name = strcat(surf_str_lower(1:3),'_',surf_str_lower(5:end),'_','DEM');
    geotiff_fn = fullfile(ct_filename_out(param,DEM_source),sprintf('%s_%03d_%s.tif',param.day_seg,frm,layer_name));
    
    surf_str_save = strrep(surf_str,' ','_');
    
    mat_file_path = fullfile(ct_filename_out(param,DEM_source),sprintf('%s_%03d_%s.mat',param.day_seg,frm,'ice_surface'));
    mat_file_name = load(mat_file_path);
    doa_points_coordinates{frm_idx} = mat_file_name.doa_points_coordinates;
    
    points{frm_idx}.x = (doa_points_coordinates{frm_idx}.x(33,:))./1e3;
    points{frm_idx}.y = (doa_points_coordinates{frm_idx}.y(33,:))./1e3;
    
    %     lat = data_tmp.Latitude;
    %     lon = data_tmp.Longitude;
    %     proj = geotiffinfo(geotiff_fn);
    %     [points{frm_idx}.x,points{frm_idx}.y] = projfwd(proj,lat,lon);
    %     points{frm_idx}.x = (points{frm_idx}.x)./1e3;
    %     points{frm_idx}.y = (points{frm_idx}.y)./1e3;
    
    [DEM{frm_idx}, R, tmp] = geotiffread(geotiff_fn);
    %   DEM_both{frm_idx} = DEM;
    
    surf_str_save = strrep(surf_str,' ','_');
    if plot_flag
      f = figure(frm_idx); clf;
      h_img = imagesc( (R(3,1) + R(2,1)*(0:size(DEM{frm_idx},2)-1))/1e3, (R(3,2) + R(1,2)*(0:size(DEM{frm_idx},1)-1))/1e3, DEM{frm_idx});
      hold on
      plot(points{frm_idx}.x,points{frm_idx}.y)
      set(gca,'YDir','normal');
      xlabel('X (km)');
      ylabel('Y (km)');
      Title = strcat('Frame #',num2str(frm),': Original DEM');
      title(Title)
      hcolor = colorbar;
      set(get(hcolor,'YLabel'),'String',sprintf('%s Elevation (WGS-84,m)',surf_str))
    end
    if save_flag
      % Save
      out_fn_name = sprintf('%s_%03d_%s_%s',param.day_seg,frm,surf_str_save,'Original');
      saveas(f,[fullfile(out_dir,out_fn_name),'.fig']);
      saveas(f,[fullfile(out_dir,out_fn_name),'.jpg']);
    end
    
    % x_rng/y_rng: the range of the DEM along the x/y-axis.
    x_rng{frm_idx} = (R(3,1) + R(2,1)*(0:size(DEM{frm_idx},2)-1))/1e3;
    y_rng{frm_idx} = (R(3,2) + R(1,2)*(0:size(DEM{frm_idx},1)-1))/1e3;
    % x_rng_lims/y_rng_lims: used later to extract the crossovers
    x_rng_lims{frm_idx} = [x_rng{frm_idx}(1) x_rng{frm_idx}(end)];
    y_rng_lims{frm_idx} = [y_rng{frm_idx}(1) y_rng{frm_idx}(end)];
    
    if 0
      %For debugging
      figure(1000);clf;
      plot(doa_points_coordinates{frm_idx}.x./1e3,doa_points_coordinates{frm_idx}.y./1e3,'b')
      hold on
      plot(doa_points_coordinates{frm_idx}.x(33,1)./1e3,doa_points_coordinates{frm_idx}.y(33,1)./1e3,'rx','MarkerSize',10)
      plot(doa_points_coordinates{frm_idx}.x(33,3332/2)./1e3,doa_points_coordinates{frm_idx}.y(33,3332/2)./1e3,'rx','MarkerSize',10)
      xlabel('Along-track (Km)')
      ylabel('Across-track (Km)')
      grid on
    end
  end
  
  %% Draw the intersecting areas between the two frames separately
  % Extract the limits of the matched DEMs matrices
  st_x  = max(x_rng_lims{1}(1),x_rng_lims{2}(1));
  end_x = min(x_rng_lims{1}(2),x_rng_lims{2}(2));
  st_y  = max(y_rng_lims{1}(2),y_rng_lims{2}(2));
  end_y = min(y_rng_lims{1}(1),y_rng_lims{2}(1));
  
  lim_x = [st_x end_x];
  lim_y = [st_y end_y];
  for frm_idx = 1:length(frms)
    frm = frms(frm_idx);
    DEM_match = [];
    x_DEM = find(x_rng{frm_idx}>lim_x(1) & x_rng{frm_idx}<lim_x(2));
    y_DEM = find(y_rng{frm_idx}>lim_y(1) & y_rng{frm_idx}<lim_y(2));
    
    for idx_x=1:length(y_DEM)
      for idx_y = 1:length(x_DEM)
        DEM_match(idx_x,idx_y) = DEM{frm_idx}(y_DEM(idx_x),x_DEM(idx_y));
      end
    end
    DEM_match_tmp{frm_idx} = DEM_match;
    
    if 0
      % Plot the matched area
      f = figure(frm_idx+200);clf
      h_img = imagesc(x_rng{frm_idx}(x_DEM),y_rng{frm_idx}(y_DEM),DEM_match);
      set(gca,'YDir','normal');
      xlabel('X (km)');
      ylabel('Y (km)');
      Title = strcat('Frame #',num2str(frm),': DEM of the matched area');
      title(Title)
      hcolor = colorbar;
      set(get(hcolor,'YLabel'),'String',sprintf('%s Elevation (WGS-84,m)',surf_str))
      
      if save_flag
        % Save
        out_fn_name = sprintf('%s_%03d_%s_%s',param.day_seg,frm,surf_str_save,'Matched');
        saveas(f,[fullfile(out_dir,out_fn_name),'.fig']);
        saveas(f,[fullfile(out_dir,out_fn_name),'.jpg']);
      end
    end
  end
  % Make sure the sizes match
  if size(DEM_match_tmp{1},1) ~= size(DEM_match_tmp{2},1)
    DEM_match_tmp{1} = DEM_match_tmp{1}(1:min(size(DEM_match_tmp{1},1),size(DEM_match_tmp{2},1)),:);
    DEM_match_tmp{2} = DEM_match_tmp{2}(1:min(size(DEM_match_tmp{1},1),size(DEM_match_tmp{2},1)),:);
  end
  if size(DEM_match_tmp{1},2) ~= size(DEM_match_tmp{2},2)
    DEM_match_tmp{1} = DEM_match_tmp{1}(:,1:min(size(DEM_match_tmp{1},2),size(DEM_match_tmp{2},2)));
    DEM_match_tmp{2} = DEM_match_tmp{2}(:,1:min(size(DEM_match_tmp{1},2),size(DEM_match_tmp{2},2)));
  end
  
  %   DEM_match_tmp{1} = double(DEM_match_tmp{1});
  %   DEM_match_tmp{2} = double(DEM_match_tmp{2});
  
  DEM_nan_frm_1_idx = DEM_match_tmp{1}<=-9999;
  DEM_nan_frm_2_idx = DEM_match_tmp{2}<=-9999;
  DEM_match_tmp{1}(DEM_nan_frm_1_idx) = NaN;
  DEM_match_tmp{2}(DEM_nan_frm_2_idx) = NaN;
  DEM_match_sum = DEM_match_tmp{1}+DEM_match_tmp{2};
  
  %   DEM_match_tmp{1} = int16(DEM_match_tmp{1});
  %   DEM_match_tmp{2} = int16(DEM_match_tmp{2});
  
  x_lim = x_rng{2}((x_DEM));
  y_lim = y_rng{2}((y_DEM));
  
  [r,c] = find(~isnan(DEM_match_sum));
  
  if isempty(r) || isempty(c)
    failed_intersections_idx = failed_intersections_idx+1;
    empty_intersections{failed_intersections_idx} = frm_inters{frm_inter_idx};
    fprintf('\nWarning: No intersection found between frame %s and frame %s\n',...
      frm_inters{frm_inter_idx}{1},frm_inters{frm_inter_idx}{2})
    fprintf('\n Skipping ...\n')
    continue
  end
  
  r_min = min(max(r(find(c==c(end)))),max(r(find(c==c(1)))));
  r_max = max(min(r(find(c==c(end)))),min(r(find(c==c(1)))));
  %   r_min = min(r(find(c==c(end))),r(find(c==c(1))));
  %   r_min = max(r_min);
  %   r_max = max(r(find(c==c(end))),r(find(c==c(1))));
  %   r_max = min(r_max);
  
  c_min = find(c==c(1));
  c_min = min(c(c_min));
  c_max = find(c==c(end));
  c_max = min(c(c_max));
  
  x_min = x_lim(c_min);
  y_min = y_lim(r_max);
  x_max = x_lim(c_max);
  y_max = y_lim(r_min);
  
  DEM_overlap_center = [(x_max+x_min)/2 (y_max+y_min)/2];
  
  if 0
    figure(frm_idx+300);clf
    h_img = imagesc(x_rng{frm_idx}(x_DEM),y_rng{frm_idx}(y_DEM),DEM_match_sum);
    set(gca,'YDir','normal');
    xlabel('X (km)');
    ylabel('Y (km)');
    Title = strcat('Frame #',num2str(frm),': DEM of the sum area');
    title(Title)
    hcolor = colorbar;
    set(get(hcolor,'YLabel'),'String',sprintf('%s Elevation (WGS-84,m)',surf_str))
  end
  
  for frm_idx = 1:length(frms)
    x_DEM_overlap = find(x_rng{frm_idx}>x_min & x_rng{frm_idx}<x_max);
    y_DEM_overlap = find(y_rng{frm_idx}>y_min & y_rng{frm_idx}<y_max);
    
    DEM_overlap{frm_idx} = DEM_match_tmp{frm_idx}(r_min:r_max,c(1):c(end));
    %   DEM_overlap{frm_idx} = DEM_match_tmp{frm_idx}(r(1):length(y_DEM_overlap),c(1):length(x_DEM_overlap));
    x_rng_DEM_overlap{frm_idx} = x_rng{frm_idx}(x_DEM_overlap);
    y_rng_DEM_overlap{frm_idx} = y_rng{frm_idx}(y_DEM_overlap);
    
    points_x_lim = find(points{frm_idx}.x>x_min & points{frm_idx}.x<x_max);
    points_y_lim = find(points{frm_idx}.y>y_min & points{frm_idx}.y<y_max);
    %     figure;plot(points{frm_idx}.x(points_x_lim(1:length(points_y_lim))),points{frm_idx}.y(points_y_lim))
    %     min_val = min(length(x_DEM_overlap),length(y_DEM_overlap));
    %     if min_val ~= length(x_DEM_overlap)
    %       [min_pt_val_x,min_pt_idx_x] = min(abs(points{frm_idx}.x-DEM_overlap_center(1)));
    %     else
    %       [min_pt_val_y,min_pt_idx_y] = min(abs(points{frm_idx}.y-DEM_overlap_center(2)));
    %     end
    
    [~,nearest_pt_idx(frm_idx)] = min(sqrt((points{frm_idx}.x-DEM_overlap_center(1)).^2+(points{frm_idx}.y-DEM_overlap_center(2)).^2));
    nearest_pt{frm_idx} = [points{frm_idx}.x(nearest_pt_idx(frm_idx)),points{frm_idx}.y(nearest_pt_idx(frm_idx))] ;
    
    %     nearest_pt_A = [points{frm_idx}.x(nearest_pt_idx(frm_idx)-5),points{frm_idx}.y(nearest_pt_idx(frm_idx)-5)];
    %     nearest_pt_B = [points{frm_idx}.x(nearest_pt_idx(frm_idx)+5),points{frm_idx}.y(nearest_pt_idx(frm_idx)+5)];
    %     nearest_pt_diff = nearest_pt_B-nearest_pt_A;
    
    %     min_len = min(length(x_rng_DEM_overlap{frm_idx}),length(y_rng_DEM_overlap{frm_idx}));
    
    %     null_line = 1.4*null(nearest_pt_diff)';
    %
    %     normal_up{frm_idx} = [nearest_pt{frm_idx}(1),nearest_pt{frm_idx}(2)] + null_line;
    %     normal_dn{frm_idx} = [nearest_pt{frm_idx}(1),nearest_pt{frm_idx}(2)] - null_line;
  end
  
  % set all non overlapped parts to NaN in case the DEMs are not rectangles
  if ~isempty(find(isnan(DEM_overlap{1})))
    nan_idx = find(isnan(DEM_overlap{1}));
    DEM_overlap{2}(nan_idx) = NaN;
  end
  
  if ~isempty(find(isnan(DEM_overlap{2})))
    nan_idx = find(isnan(DEM_overlap{2}));
    DEM_overlap{1}(nan_idx) = NaN;
  end
  
  
  for frm_idx = 1:length(frms)
    if plot_flag
      frm = frms(frm_idx);
      nearest_x = nearest_pt{frm_idx}(1);
      nearest_y = nearest_pt{frm_idx}(2);
      
      f = figure(frm_idx+400);clf
      h_img = imagesc(x_rng_DEM_overlap{frm_idx},y_rng_DEM_overlap{frm_idx},DEM_overlap{frm_idx});
      hold on
      plot(points{frm_idx}.x,points{frm_idx}.y,'Linewidth',1.5)
      %       plot(DEM_overlap_center(1),DEM_overlap_center(2),'k','Linewidth',10)
      scatter(DEM_overlap_center(1),DEM_overlap_center(2),[],'k','filled')
      %      scatter(nearest_x,nearest_y,[],'r','filled')
      plot(nearest_x,nearest_y,'rx','Linewidth',2,'MarkerSize',10)
      set(gca,'YDir','normal');
      xlabel('X (km)');
      ylabel('Y (km)');
      Title = strcat('Frame #',num2str(frm),': DEM of the overlapped area');
      title(Title)
      hcolor = colorbar;
      set(get(hcolor,'YLabel'),'String',sprintf('%s Elevation (WGS-84,m)',surf_str))
    end
    if save_flag
      % Save
      out_fn_name = sprintf('%s_%03d_%s_%s',param.day_seg,frm,surf_str_save,'Overlapped');
      saveas(f,[fullfile(out_dir,out_fn_name),'.fig']);
      saveas(f,[fullfile(out_dir,out_fn_name),'.jpg']);
    end
  end
  
  %% Magnitude height errors
  % Make sure the sizes match
  if size(DEM_overlap{1},1) ~= size(DEM_overlap{2},1)
    DEM_overlap{1} = DEM_overlap{1}(1:min(size(DEM_overlap{1},1),size(DEM_overlap{2},1)),:);
    DEM_overlap{2} = DEM_overlap{2}(1:min(size(DEM_overlap{1},1),size(DEM_overlap{2},1)),:);
  end
  if size(DEM_overlap{1},2) ~= size(DEM_overlap{2},2)
    DEM_overlap{1} = DEM_overlap{1}(:,1:min(size(DEM_overlap{1},2),size(DEM_overlap{2},2)));
    DEM_overlap{2} = DEM_overlap{2}(:,1:min(size(DEM_overlap{1},2),size(DEM_overlap{2},2)));
  end
  DEM_overlap_mag_error = abs(DEM_overlap{1}-DEM_overlap{2});
  Errors = sort(DEM_overlap_mag_error(:),'ascend');
  
  if plot_flag
    % plot: the x_rng and y_rng are from the second frame
    f = figure(frm_idx+500);clf
    h_img = imagesc(x_rng_DEM_overlap{2},y_rng_DEM_overlap{2},DEM_overlap_mag_error);
    set(gca,'YDir','normal');
    xlabel('X (km)');
    ylabel('Y (km)');
    Title = strcat('DEM-overlap magnitude error');
    title(Title)
    hcolor = colorbar;
    if strcmp(layer ,'bottom')
      layer_text = 'Bed Height Error (m)';
    elseif strcmp(layer ,'ice surface')
      layer_text = 'Surface Height Error (m)';
    end
    set(get(hcolor,'YLabel'),'String',layer_text)
  end
  if save_flag
    % Save
    out_fn_name = sprintf('%s_%03d_%s_%03d_%s_%s',param.day_seg,frms(1),'&',frms(2),surf_str_save,'Error');
    saveas(f,[fullfile(out_dir,out_fn_name),'.fig']);
    saveas(f,[fullfile(out_dir,out_fn_name),'.jpg']);
  end
  
  if plot_flag
    % Histogram plot
    f= figure(frm_idx+600);clc
    histogram(DEM_overlap_mag_error(:));
    Title = 'Histogram of the DEM overlap magnitude error';
    title(Title)
  end
  if save_flag
    % Save
    out_fn_name = sprintf('%s_%03d_%s_%03d_%s_%s',param.day_seg,frms(1),'&',frms(2),surf_str_save,'Hist');
    saveas(f,[fullfile(out_dir,out_fn_name),'.fig']);
    saveas(f,[fullfile(out_dir,out_fn_name),'.jpg']);
  end
  
  if 0
    %surf plot
    [x,y] = meshgrid(x_rng(x_DEM),y_rng(y_DEM));
    figure(400);
    mesh(x,y,DEM_overlap_both{1})
    colorbar
    figure(401);
    mesh(x,y,DEM_overlap_both{2})
    colorbar
    xlabel('along-track (Km)')
    xlabel('Cross-track (Km)')
    zlabel('Bed height (m)')
    
    figure(500);
    surf(x,y,DEM_overlap_both{1})
    colorbar
    figure(502);
    surf(x,y,DEM_overlap_both{2})
    colorbar
  end
  
  
  %% Plot all figures together
  slice = nearest_pt_idx;
  surf_layer_idx = strmatch(layers_plot(1),{surfdata_all{1}.name},'exact');
  if isempty(surf_layer_idx)
    surf_layer_idx = strmatch('top',{surfdata_all{1}.name},'exact');
  end
  bott_layer_idx = strmatch(layers_plot(2),{surfdata_all{1}.name},'exact');
  doa_trim = [3 4];
  center_doa_bin = size(doa_points_coordinates{2}.x,1)/2 + 1;
  
  % Flight line coordinates for frame 1
  flight_line_x1 = doa_points_coordinates{1}.x(center_doa_bin,:)./1e3;
  flight_line_y1 = doa_points_coordinates{1}.y(center_doa_bin,:)./1e3;
  
  %   flight_line_tmp_x1 = flight_line_x1(flight_line_x1>=x_rng_DEM_overlap{1}(1) & flight_line_x1<=x_rng_DEM_overlap{1}(end));
  %   flight_line_tmp_y1 = flight_line_y1(flight_line_y1<=y_rng_DEM_overlap{1}(1) & flight_line_y1>=y_rng_DEM_overlap{1}(end));
  %
  %   flight_line_x1 = flight_line_x1(1:min(length(flight_line_x1),length(flight_line_y1)));
  %   flight_line_y1 = fliplr(flight_line_y1(1:min(length(flight_line_x1),length(flight_line_y1))));
  
  % Flight line coordinates for frame 2
  flight_line_x2 = doa_points_coordinates{2}.x(center_doa_bin,:)./1e3;
  flight_line_y2 = doa_points_coordinates{2}.y(center_doa_bin,:)./1e3;
  
  %   flight_line_tmp_x2 = flight_line_x2(flight_line_x2>=x_rng_DEM_overlap{2}(1) & flight_line_x2<=x_rng_DEM_overlap{2}(end));
  %   flight_line_tmp_y2 = flight_line_y2(flight_line_y2<=y_rng_DEM_overlap{2}(1) & flight_line_y2>=y_rng_DEM_overlap{2}(end));
  %
  %   flight_line_x2 = flight_line_x2(1:min(length(flight_line_x2),length(flight_line_y2)));
  %   flight_line_y2 = fliplr(flight_line_y2(1:min(length(flight_line_x2),length(flight_line_y2))));
  
  %   idx = find(abs(flight_line_tmp_x1(1:min(length(flight_line_tmp_x1),length(flight_line_tmp_x2)))...
  %     -flight_line_tmp_x2(1:min(length(flight_line_tmp_x1),length(flight_line_tmp_x2))))<=0.1,1);
  %   px1 = flight_line_x1(idx);
  %   px2 = flight_line_x2(idx);
  %
  %   py1 = flight_line_y1(idx);
  %   py2 = flight_line_y2(idx);
  %    slice = [idx,idx];
  
  
  %   doa_points_coordinates{1}.x([1:doa_trim(1), end-doa_trim(2)+1:end],slice) = NaN;
  %   doa_points_coordinates{1}.y([1:doa_trim(1), end-doa_trim(2)+1:end],slice) = NaN;
  %
  %   doa_points_coordinates{2}.x([1:doa_trim(1), end-doa_trim(2)+1:end],slice) = NaN;
  %   doa_points_coordinates{2}.y([1:doa_trim(1), end-doa_trim(2)+1:end],slice) = NaN;
  
  f_all = figure(22);clf
  % Plot DEM of frm 1
  subplot(231);
  h_img = imagesc(x_rng_DEM_overlap{1},y_rng_DEM_overlap{1},DEM_overlap{1});
  hold on
  % Plot the flight line
  plot(flight_line_x1,flight_line_y1(1:length(flight_line_x1)),'k','LineWidth',1.5);
  %   plot(points{1}.x,points{1}.y,'k','LineWidth',1.5);
  
  % Plot the slice location
  plot(doa_points_coordinates{1}.x(center_doa_bin,slice(1))./1e3,doa_points_coordinates{1}.y(center_doa_bin,slice(1))./1e3,'x','MarkerSize',10,'MarkerEdgeColor',[1/2 0 1/2])
  plot(doa_points_coordinates{1}.x(1:center_doa_bin-1,slice(1))./1e3,doa_points_coordinates{1}.y(1:center_doa_bin-1,slice(1))./1e3,'wx','MarkerSize',10)
  plot(doa_points_coordinates{1}.x(center_doa_bin+1:end,slice(1))./1e3,doa_points_coordinates{1}.y(center_doa_bin+1:end,slice(1))./1e3,'rx','MarkerSize',10)
  %   plot(doa_points_coordinates{1}.x(:,slice(1))./1e3,doa_points_coordinates{1}.y(:,slice(1))./1e3,'wx','MarkerSize',10)
  
  % Plot the right most DoA location
  %   plot(doa_points_coordinates{1}.x(end,slice(1))./1e3,doa_points_coordinates{1}.y(end,slice(1))./1e3,'ro','MarkerSize',10)
  
  legend('Flight path',sprintf('%s %d','Slice', slice(1)))
  %   legend({'Flight path',['Slice location' char(10) '(white is left, red is right)']})
  %   legend(h_plot,'Flight path','Location','northoutside');
  set(gca,'YDir','normal');
  xlabel('X (km)');
  ylabel('Y (km)');
  %   Title = sprintf('Frame %d',frms(1));
  
  day_seg_tmp = param.day_seg;
  day_tmp = day_seg_tmp(1:end-3);
  seg_tmp = day_seg_tmp(end-2:end);
  day_seg = strcat(day_tmp,'\',seg_tmp);
  Title = strcat(day_seg,'\_',sprintf('%03d',frms(1)));
  title(Title)
  %   hcolor = colorbar;
  %   hcolor.Location = 'westoutsid';
  %
  %   set(get(hcolor,'YLabel'),'String',sprintf('%s Elevation (WGS-84,m)',surf_str))
  cmin_frm1 = double(min(DEM_overlap{1}(find(DEM_overlap{1}~=0))));
  cmax_frm1 = double(max(DEM_overlap{1}(find(DEM_overlap{1}~=0))));
  caxis([cmin_frm1 cmax_frm1])
  
  % Plot DEM of frm 2
  subplot(232)
  
  h_img = imagesc(x_rng_DEM_overlap{2},y_rng_DEM_overlap{2},DEM_overlap{2});
  hold on
  % Plot the flight line
  plot(flight_line_x2,flight_line_y2,'k','LineWidth',1.5);
  %   plot(points{2}.x,points{2}.y,'k','LineWidth',1.5);
  
  % Plot the slice location
  % Colors: gray = [1/2 1/2 1/2], orange = [1 165/255 0], purple = [1/2 0
  % 1/2], yellow = [1 1 0]
  plot(doa_points_coordinates{2}.x(center_doa_bin,slice(2))./1e3,doa_points_coordinates{2}.y(center_doa_bin,slice(2))./1e3,'x','MarkerSize',10,'MarkerEdgeColor',[1/2 0 1/2])
  plot(doa_points_coordinates{2}.x(1:center_doa_bin-1,slice(2))./1e3,doa_points_coordinates{2}.y(1:center_doa_bin-1,slice(2))./1e3,'wx')
  plot(doa_points_coordinates{2}.x(center_doa_bin+1:end,slice(2))./1e3,doa_points_coordinates{2}.y(center_doa_bin+1:end,slice(2))./1e3,'rx','MarkerSize',10)
  %   plot(doa_points_coordinates{2}.x(:,slice(2))./1e3,doa_points_coordinates{2}.y(:,slice(2))./1e3,'rx','MarkerSize',10)
  
  % Plot the right most DoA location
  %   plot(doa_points_coordinates{2}.x(end,slice(2))./1e3,doa_points_coordinates{2}.y(end,slice(2))./1e3,'ro','MarkerSize',10)
  
  legend('Flight path',sprintf('%s %d','Slice', slice(2)))
  %   legend({'Flight path',['Slice location' char(10) '(white is left, red is right)']} )
  %   legend('Flight path','Slice location')
  set(gca,'YDir','normal');
  xlabel('X (km)');
  %   ylabel('Y (km)');
  %   Title = sprintf('Frame %d',frms(2));
  
  Title = strcat(day_seg,'\_',sprintf('%03d',frms(2)));
  title(Title)
  
  cmin_frm2 = double(min(DEM_overlap{2}(find(DEM_overlap{2}~=0))));
  cmax_frm2 = double(max(DEM_overlap{2}(find(DEM_overlap{2}~=0))));
  caxis([cmin_frm2 cmax_frm2])
  
  % Plot the error DEM
  subplot(233)
  h_img = imagesc(x_rng_DEM_overlap{2},y_rng_DEM_overlap{2},DEM_overlap_mag_error);
  hold on
  % Plot both flight lines
  plot(flight_line_x1,flight_line_y1,'k','LineWidth',1.5);
  plot(flight_line_x2,flight_line_y2,'k','LineWidth',1.5);
  %   plot(doa_points_coordinates{1}.x(center_doa_bin,:)./1e3,doa_points_coordinates{1}.y(center_doa_bin,:)./1e3,'k','LineWidth',1.5);
  %   plot(doa_points_coordinates{2}.x(center_doa_bin,:)./1e3,doa_points_coordinates{2}.y(center_doa_bin,:)./1e3,'k','LineWidth',1.5);
  %
  set(gca,'YDir','normal');
  xlabel('X (km)');
  %   ylabel('Y (km)');
  title('|Error|')
  if strcmp(layer ,'bottom')
    layer_text = 'Bed Height Error (m)';
  elseif strcmp(layer ,'ice surface')
    layer_text = 'Surface Height Error (m)';
  end
  hcolor = colorbar;
  hcolor.Location = 'eastoutsid';
  
  set(get(hcolor,'YLabel'),'String',layer_text)
  if strcmp(layer ,'bottom')
    cmax = min(100,max(double(DEM_overlap_mag_error(:))));
    if cmax == 0
      cmax = 0.01;
    end
    caxis([0 cmax])
  end
  
  % Plot slice from frm 1
  subplot(234)
  imagesc(10*log10(data{1}(:,:,slice(1))))
  hold on
  
  surf_layer = surfdata_all{1}(surf_layer_idx).y(:,slice(1));
  bott_layer = surfdata_all{1}(bott_layer_idx).y(:,slice(1));
  
  surf_layer([1:doa_trim(1), end-doa_trim(2)+1:end]) = NaN;
  bott_layer([1:doa_trim(1), end-doa_trim(2)+1:end]) = NaN;
  
  surf_layer_left = surf_layer;
  surf_layer_left(35:end) = NaN;
  surf_layer_right = surf_layer;
  surf_layer_right(1:33) = NaN;
  
  bott_layer_left = bott_layer;
  bott_layer_left(35:end) = NaN;
  bott_layer_right = bott_layer;
  bott_layer_right(1:33) = NaN;
  
  min_surf_bin = nanmin(surf_layer);
  max_bott_bin = nanmax(bott_layer);
  
  plot(surf_layer_left,'w','LineWidth',1.5);
  plot(surf_layer_right,'r','LineWidth',1.5);
  
  plot(bott_layer_left,'w','LineWidth',1.5);
  plot(bott_layer_right,'r','LineWidth',1.5);
  
  
  %   plot(surf_layer,'w','LineWidth',1.5);
  %   plot(bott_layer,'r','LineWidth',1.5);
  
  %   plot(length(surf_layer)-doa_trim(2),surf_layer(end-doa_trim(2)),'ro');
  %   plot(length(bott_layer)-doa_trim(2),bott_layer(end-doa_trim(2)),'ro');
  
  %   legend('Surface','Bottom')
  %   legend({['Top curve is ice surface' char(10), 'Bottom curve is ice-bottom']})
  
  set(gca,'YDir','reverse');
  xlabel('Direction of arrival bin');
  ylabel('Range bin');
  Title = strcat(day_seg,'\_',sprintf('%03d: %s %d',frms(1),'Slice', slice(1)));
  %   Title = sprintf('Frame %d: Slice %d',frms(1),slice(1));
  title(Title)
  xlim([1 64])
  ylim([min_surf_bin-50 max_bott_bin+50])
  
  % Plot slice from frm 2
  subplot(235)
  imagesc(10*log10(data{2}(:,:,slice(2))))
  hold on
  surf_layer = surfdata_all{2}(surf_layer_idx).y(:,slice(2));
  bott_layer = surfdata_all{2}(bott_layer_idx).y(:,slice(2));
  
  surf_layer([1:doa_trim(1), end-doa_trim(2)+1:end]) = NaN;
  bott_layer([1:doa_trim(1), end-doa_trim(2)+1:end]) = NaN;
  
  surf_layer_left = surf_layer;
  surf_layer_left(35:end) = NaN;
  surf_layer_right = surf_layer;
  surf_layer_right(1:33) = NaN;
  
  bott_layer_left = bott_layer;
  bott_layer_left(35:end) = NaN;
  bott_layer_right = bott_layer;
  bott_layer_right(1:33) = NaN;
  
  min_surf_bin = nanmin(surf_layer);
  max_bott_bin = nanmax(bott_layer);
  
  plot(surf_layer_left,'w','LineWidth',1.5);
  plot(surf_layer_right,'r','LineWidth',1.5);
  
  plot(bott_layer_left,'w','LineWidth',1.5);
  plot(bott_layer_right,'r','LineWidth',1.5);
  
  %   plot(surf_layer,'w','LineWidth',1.5);
  %   plot(bott_layer,'r','LineWidth',1.5);
  %
  %   plot(length(surf_layer)-doa_trim(2),surf_layer(end-doa_trim(2)),'ro');
  %   plot(length(bott_layer)-doa_trim(2),bott_layer(end-doa_trim(2)),'ro');
  %   legend('Surface','Bottom')
  %   legend({['Top curve is ice surface' char(10), 'Bottom curve is ice-bottom']})
  set(gca,'YDir','reverse');
  xlabel('Direction of arrival bin');
  %   ylabel('Range bin');
  %   Title = sprintf('Frame %d: Slice %d',frms(2),slice(2));
  Title = strcat(day_seg,'\_',sprintf('%03d: %s %d',frms(2),'Slice', slice(2)));
  title(Title)
  xlim([1 64])
  ylim([min_surf_bin-50 max_bott_bin+50])
  
  % Plot the CDF
  subplot(236)
  trim_err_lim = 100;% in meter
  Errors = Errors(1:find(isnan(Errors),1)-1);
  %   err = find(DEM_overlap_mag_error_vec<=trim_err_lim);
  %   [cdf_vals,cdf_pts] = ecdf(DEM_overlap_mag_error_vec(err));
  [cdf_vals,cdf_pts] = ecdf(Errors);
  
  plot(cdf_pts,cdf_vals,'LineWidth',1.5)
  xlabel('|Error|(m)')
  ylabel('F(|Error|)')
  title('cdf(|Error|)')
  grid on
  if max(cdf_pts) ~= 0
    xlim([0 trim_err_lim])
    %     xlim([0 min(max(cdf_pts),trim_err_lim)])
  end
  ylim([0 1])
  
  day_seg_tmp = param.day_seg;
  day_tmp = day_seg_tmp(1:end-3);
  seg_tmp = day_seg_tmp(end-2:end);
  day_seg = strcat(day_tmp,'\',seg_tmp);
  %   Title = strcat(day_seg,'\_',sprintf('%03d%s',frms(1),','),day_seg,'\_',sprintf('%03d',frms(2)));
  %   suptitle(Title)
  
  % Plot the colorbar separately
  f_colorbar = figure(33);
  colorbar
  caxis([nanmin(cmin_frm1,cmin_frm2) nanmax(cmax_frm1,cmax_frm2)])
  
  if save_all_flag
    % Save
    %     fig = gcf;
    set(f_all, 'Units', 'Inches', 'Position', [0, 0, 16 12], 'PaperUnits', 'Inches', 'PaperSize', [16, 12])
    f_all.PaperPositionMode = 'auto';% This property preservs the saved figure from being distorted
    out_fn_name = sprintf('%s_%03d_%s_%03d_%s_%s',param.day_seg,frms(1),'&',frms(2),surf_str_save,'all');
    
    saveas(f_all,[fullfile(out_dir,out_fn_name),'.fig']);
    saveas(f_all,[fullfile(out_dir,out_fn_name),'.jpg']);
  end
  
  if colorbar_fig_flag
    out_fn_name = sprintf('%s_%03d_%s_%03d_%s',param.day_seg,frms(1),'&',frms(2),'Colorbar');
    saveas(f_colorbar,[fullfile(out_dir,out_fn_name),'.fig']);
    saveas(f_colorbar,[fullfile(out_dir,out_fn_name),'.jpg']);
  end
  %% Statistics
  % 1) error percentages
  required_error_in_meter = [10 50 100];
  for error_idx = 1:length(required_error_in_meter)
    required_error = required_error_in_meter(error_idx);
    cdf_pts_sub = find(cdf_pts<=required_error);
    cdf_pts_percent(error_idx) = (length(cdf_pts_sub)/length(cdf_vals))*100;
    fprintf('\n %3.2f percent of the errors is <=%d meters\n',cdf_pts_percent(error_idx),required_error);
  end
  fprintf('\n');
  
  % 2) mean and median
  % case a: using all the errors
  err_min = min(Errors);
  err_max = max(Errors);
  err_mean = mean(Errors);
  err_med = median(Errors);
  err_std = std(Errors);
  %   err_percent_less_than_med = find(Errors==err_med,1)/length(Errors)*100;
  
%     err_trim_percentage = [5:5:95];
  err_trim_percentage = [10];
  % case b: throw away the largest 10 percent of the errors
  for trim_idx = 1:length(err_trim_percentage)
    Errors_trimmed_idx = round((err_trim_percentage(trim_idx)/100)*length(Errors));
    Errors_trimmed = Errors(1:end-Errors_trimmed_idx);
    
    err_trim_min(trim_idx) = min(Errors_trimmed);
    err_trim_max(trim_idx) = max(Errors_trimmed);
    err_trim_mean(trim_idx) = mean(Errors_trimmed);
    err_trim_med(trim_idx) = median(Errors_trimmed);
    err_trim_std(trim_idx) = std(Errors_trimmed);
    
    % 3) root mean/median squared-error
    % case b: throw away the largest 10 percent of the errors
    root_mean_squared_err_trim(trim_idx) = sqrt(mean(Errors_trimmed.^2));
    root_med_squared_err_trim(trim_idx) = sqrt(median(Errors_trimmed.^2));
  end
  
  root_mean_squared_err_trim_tmp = root_mean_squared_err_trim + root_mean_squared_err_trim_tmp;
  
  
  % 3) root mean/median squared-error
  % case a: using all the errors
  root_mean_squared_err = sqrt(mean(Errors.^2));
  root_med_squared_err = sqrt(median(Errors.^2));
  
  if 0
    root_mean_squared_err_trim_tmp_avg = root_mean_squared_err_trim_tmp/length(frm_inters);
    figure(1000);clf
    %   f = plot(err_trim_percentage,root_mean_squared_err_trim,'LineWidth',1.5);
    f = plot(err_trim_percentage,root_mean_squared_err_trim_tmp_avg,'LineWidth',1.5);
    xlim([err_trim_percentage(1),err_trim_percentage(end)])
    ylim([min(root_mean_squared_err_trim_tmp_avg), max(root_mean_squared_err_trim_tmp_avg)])
    xlabel('Percentage of the removed largest errors')
    ylabel('RMSE')
    %   title(Title)
    grid on
    
    saveas(f,[fullfile(out_dir,'RMSE_trimmed_largest_errors'),'.fig']);
    saveas(f,[fullfile(out_dir,'RMSE_trimmed_largest_errors'),'.jpg']);
  end
  
  if 0
    % Statistics per crossover
    stats_results{frm_inter_idx}.name =...
      sprintf('%s_%03d_%s_%s_%03d_%s',param.day_seg,frms(1),'and',param.day_seg,frms(2),surf_str_save);
    stats_results{frm_inter_idx}.err_min  = err_min;
    stats_results{frm_inter_idx}.err_max  = err_max;
    stats_results{frm_inter_idx}.err_mean = err_mean;
    stats_results{frm_inter_idx}.err_med  = err_med;
    stats_results{frm_inter_idx}.err_std  = err_std;
    stats_results{frm_inter_idx}.RMeanSE  = root_mean_squared_err;
    stats_results{frm_inter_idx}.RMedSE   = root_med_squared_err;
    
    stats_results{frm_inter_idx}.err_trim_min  = err_trim_min(1);
    stats_results{frm_inter_idx}.err_trim_max  = err_trim_max(1);
    stats_results{frm_inter_idx}.err_trim_mean = err_trim_mean(1);
    stats_results{frm_inter_idx}.err_trim_med  = err_trim_med(1);
    stats_results{frm_inter_idx}.err_trim_std  = err_trim_std(1);
    stats_results{frm_inter_idx}.RMeanSE_trim  = root_mean_squared_err_trim(1);
    stats_results{frm_inter_idx}.RMedSE_trim   = root_med_squared_err_trim(1);
  end
  
  if 1
    %Statistics as vectors
    stats_results.num_err(frm_inter_idx)  = length(Errors);
    stats_results.err_min(frm_inter_idx)  = err_min;
    stats_results.err_max(frm_inter_idx)  = err_max;
    stats_results.err_mean(frm_inter_idx) = err_mean;
    stats_results.err_med(frm_inter_idx)  = err_med;
    stats_results.err_std(frm_inter_idx)  = err_std;
    stats_results.RMeanSE(frm_inter_idx)  = root_mean_squared_err;
    stats_results.RMedSE(frm_inter_idx)   = root_med_squared_err;
    
    stats_results.err_trim_min(frm_inter_idx)  = err_trim_min(1);
    stats_results.err_trim_max(frm_inter_idx)  = err_trim_max(1);
    stats_results.err_trim_mean(frm_inter_idx) = err_trim_mean(1);
    stats_results.err_trim_med(frm_inter_idx)  = err_trim_med(1);
    stats_results.err_trim_std(frm_inter_idx)  = err_trim_std(1);
    stats_results.RMeanSE_trim(frm_inter_idx)  = root_mean_squared_err_trim(1);
    stats_results.RMedSE_trim(frm_inter_idx)   = root_med_squared_err_trim(1);
  end
  
  fprintf('\n Minimum Error is %3.2f meters \n',err_min);
  fprintf('\n Maximum Error is %3.2f meters \n',err_max);
  fprintf('\n Mean Error is %3.2f meters \n',err_mean);
  fprintf('\n Median Error is %3.2f meters \n',err_med);
  fprintf('\n Error STD is %3.2f meters \n',err_std);
  fprintf('\n RMeanSE is %3.2f meters \n',root_mean_squared_err);
  fprintf('\n RMedSE is %3.2f meters \n',root_med_squared_err);
  
  fprintf('\n\n When the largest %d percent of the error is trimmed:\n',err_trim_percentage);
  fprintf('\n Minimum Error is %3.2f meters \n',err_trim_min);
  fprintf('\n Maximum Error is %3.2f meters \n',err_trim_max);
  fprintf('\n Mean Error is %3.2f meters \n',err_trim_mean);
  fprintf('\n Median Error is %3.2f meters \n',err_trim_med);
  fprintf('\n Error STD is %3.2f meters \n',err_trim_std);
  fprintf('\n RMeanSE is %3.2f meters \n',root_mean_squared_err_trim);
  fprintf('\n RMedSE is %3.2f meters \n',root_med_squared_err_trim);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% 2D-Filtering the DEMs to reduce the errors
  % FILTERING BLURES THE DEMs===>NO NEEED FOR FILTERING
  if 0
    % Hanning window
    H = kron(hanning(5),hanning(2).');
    % Boxcar window
    %   H = ones(5,2);
    %   H = kron(chebwin(15,190),chebwin(30,190)');
    H = H ./ sum(H(:));
    
    for frm_idx = 1:length(frms)
      DEM_overlap_filtered{frm_idx} = filter2(H,DEM_overlap{frm_idx});
      if plot_flag
        frm = frms(frm_idx);
        f = figure(frm_idx+700);clf
        h_img = imagesc(x_rng_DEM_overlap{frm_idx},y_rng_DEM_overlap{frm_idx},DEM_overlap_filtered{frm_idx});
        set(gca,'YDir','normal');
        xlabel('X (km)');
        ylabel('Y (km)');
        Title = strcat('Frame #',num2str(frm),': Filtered DEM of the overlapped area');
        title(Title)
        hcolor = colorbar;
        set(get(hcolor,'YLabel'),'String',sprintf('%s Elevation (WGS-84,m)',surf_str))
      end
      if save_flag
        % Save
        out_fn_name = sprintf('%s_%03d_%s_%s',param.day_seg,frm,surf_str_save,'Overlapped_Filtered');
        saveas(f,[fullfile(out_dir,out_fn_name),'.fig']);
        saveas(f,[fullfile(out_dir,out_fn_name),'.jpg']);
      end
    end
    
    DEM_overlap_mag_error_filtered = abs(DEM_overlap_filtered{1}-DEM_overlap_filtered{2});
    
    if plot_flag
      % plot: the x_rng and y_rng are from the second frame
      f = figure(frm_idx+800);clf
      h_img = imagesc(x_rng_DEM_overlap{2},y_rng_DEM_overlap{2},DEM_overlap_mag_error_filtered);
      set(gca,'YDir','normal');
      xlabel('X (km)');
      ylabel('Y (km)');
      Title = strcat('Filtered DEM-overlap magnitude error');
      title(Title)
      hcolor = colorbar;
      set(get(hcolor,'YLabel'),'String','Bed Height Error (m)')
    end
    if save_flag
      % Save
      out_fn_name = sprintf('%s_%03d_%s_%03d_%s_%s',param.day_seg,frms(1),'&',frms(2),surf_str_save,'Error_Filtered');
      saveas(f,[fullfile(out_dir,out_fn_name),'.fig']);
      saveas(f,[fullfile(out_dir,out_fn_name),'.jpg']);
    end
    
    if plot_flag
      % Histogram plot
      f= figure(frm_idx+900);clc
      histogram(DEM_overlap_mag_error_filtered(:));
      Title = 'Histogram of the filtered DEM overlap magnitude error';
      title(Title)
    end
    if save_flag
      % Save
      out_fn_name = sprintf('%s_%03d_%s_%03d_%s_%s',param.day_seg,frms(1),'&',frms(2),surf_str_save,'Hist_Filtered');
      saveas(f,[fullfile(out_dir,out_fn_name),'.fig']);
      saveas(f,[fullfile(out_dir,out_fn_name),'.jpg']);
    end
  end
end

if (failed_intersections_idx~=0)
  fprintf('\n Intersections not found between the following frames: \n')
  for idx = 1:length(empty_intersections)
    fprintf('%s and %s\n',empty_intersections{idx}{1},empty_intersections{idx}{2})
  end
end
if 0
  %% Saving
  save(fullfile(out_dir,'stats_results'),'stats_results');
end

if 0
  
  for i=1:20
    mean_val(i) = stats_results{i}.err_mean ;
    med_val(i) = stats_results{i}.err_med ;
    rmse_val(i) = stats_results{i}.RMeanSE;
  end
  [min_mean_val, min_mean_idx] = min(mean_val);
  [min_med_val, min_med_idx] = min(med_val);
  [min_rmse_val, min_rmse_idx] = min(rmse_val);
  
  [max_mean_val, max_mean_idx] = max(mean_val);
  [max_med_val, max_med_idx] = max(med_val);
  [max_rmse_val, max_rmse_idx] = max(rmse_val);
end
