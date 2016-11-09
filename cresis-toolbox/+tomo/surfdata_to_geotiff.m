function surfdata_to_geotiff(param)

  physical_constants;

  frms = param.cmd.frms;
  
  active_surfs = param.active_surfs;
  
  out_dir = ct_filename_out(param,'DEM');
  if ~isdir(out_dir)
    mkdir(out_dir);
  end
  
  for frm_idx = 1:length(frms)
    
    frm = frms(frm_idx);
    
    surfdata_fn = fullfile(ct_filename_out(param,'surfData'),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
    data_fn = fullfile(ct_filename_out(param,'CSA_music'),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
    
    %% Load surface and bottom data
    load(surfdata_fn)
    surf_names = {surf.name};
    
    if ischar(active_surfs)
      active_surfs = find(strcmp(surf_names,active_surfs),1);
      if isempty(active_surfs)
        warning('param.active_surf does not contain a valid string for frame %s_%03d',param.day_seg,frm);
        continue;
      end
    elseif isnumeric(active_surfs) && length(active_surfs)==1
      if round(active_surfs)<=length(surf_names) && round(active_surfs)>=1
        active_surfs = round(active_surfs);
      else
        warning('param.active_surf does not contain a valid value for frame %s_%03d',param.day_seg,frm);
        continue;
      end
    elseif isnumeric(active_surfs)
      mask = round(active_surfs)<=length(surf_names) & round(active_surfs)>=1;
      active_surfs = active_surfs(mask);
    elseif iscell(active_surfs)
      act = [];
      for i = 1:numel(active_surfs)
        if ischar(active_surfs{i})
          a = find(strcmp(surf_names,active_surfs{i}),1);
          if ~isempty(a)
            act(end+1) = a;
          end
        end
      end
      active_surfs = act;
    end
    
    ice_surf_idx = find(strcmp(surf_names,'ice surface'),1);
    bot_idx = find(strcmp(surf_names,'bottom'),1);
    
    % Ensure non-negative ice thickness
    surf(bot_idx).y(surf(bot_idx).y<surf(ice_surf_idx).y) = surf(ice_surf_idx).y(surf(bot_idx).y<surf(ice_surf_idx).y);

    %% Load echogram
    mdata = load(data_fn,'theta','Time','GPS_time','param_combine');
    
    if(~isfield(mdata,'theta'))
      global gRadar
      theta = load(fullfile(gRadar.out_path,'rds','2014_Greenland_P3','CSARP_CSA_music','20140401_03','Data_20140401_03_037'),'theta');
      mdata.theta = theta.theta;
    end

    % Convert surface from range bins to twtt
    ice_surface = interp1(1:length(mdata.Time), mdata.Time, surf(ice_surf_idx).y);
    ice_bottom = interp1(1:length(mdata.Time), mdata.Time, surf(bot_idx).y);
    ice_bottom(ice_bottom<ice_surface) = ice_surface(ice_bottom<ice_surface);

    %% Convert doa,twtt to radar FCS
    % FCS: flight coordinate system (aka SAR coordinate system)
    Nx = length(mdata.GPS_time);
    theta = reshape(mdata.theta,[length(mdata.theta) 1]);
    DOA_trim = 0;
    if isfield(param,'DOA_trim')
      DOA_trim = param.DOA_trim;
    end
    if all(all(isnan(ice_surface)))
      ice_surface = zeros(size(ice_surface));
    end
    [y_bot,z_bot] = tomo.twtt_doa_to_yz(repmat(theta,[1 Nx]),theta(DOA_trim+1:end-DOA_trim,:),ice_surface(DOA_trim+1:end-DOA_trim,:),3.15,ice_bottom);
    [y_surf,z_surf] = tomo.twtt_doa_to_yz(repmat(theta,[1 Nx]),theta(DOA_trim+1:end-DOA_trim,:),ice_surface(DOA_trim+1:end-DOA_trim,:),3.15,ice_surface);
    
    for act_idx = 1:numel(active_surfs)
      
      active_surf = active_surfs(act_idx);
    
      if active_surf == ice_surf_idx
        y_active = y_surf;
        z_active = z_surf;
      elseif active_surf == bot_idx
        y_active = y_bot;
        z_active = z_bot;
      end

      %% Convert from radar FCS to ECEF
      x_plane = zeros(size(y_active));
      y_plane = zeros(size(y_active));
      z_plane = zeros(size(y_active));
      for rline = 1:size(y_active,2)
        x_plane(:,rline) = mdata.param_combine.array_param.fcs{1}{1}.origin(1,rline) ...
          + mdata.param_combine.array_param.fcs{1}{1}.y(1,rline) * y_active(:,rline) ...
          + mdata.param_combine.array_param.fcs{1}{1}.z(1,rline) * z_active(:,rline);
        y_plane(:,rline) = mdata.param_combine.array_param.fcs{1}{1}.origin(2,rline) ...
          + mdata.param_combine.array_param.fcs{1}{1}.y(2,rline) * y_active(:,rline) ...
          + mdata.param_combine.array_param.fcs{1}{1}.z(2,rline) * z_active(:,rline);
        z_plane(:,rline) = mdata.param_combine.array_param.fcs{1}{1}.origin(3,rline) ...
          + mdata.param_combine.array_param.fcs{1}{1}.y(3,rline) * y_active(:,rline) ...
          + mdata.param_combine.array_param.fcs{1}{1}.z(3,rline) * z_active(:,rline);
      end

      %% Convert from ECEF to geodetic
      [points.lat,points.lon,points.elev] = ecef2geodetic(x_plane,y_plane,z_plane,WGS84.ellipsoid);
      points.lat = points.lat * 180/pi;
      points.lon = points.lon * 180/pi;

      points.lat = points.lat(DOA_trim+1:end-DOA_trim,:);
      points.lon = points.lon(DOA_trim+1:end-DOA_trim,:);
      points.elev = points.elev(DOA_trim+1:end-DOA_trim,:);

      %% Convert from geodetic to projection

      % Load geotiff projection
      geotiff_fn = param.geotiff_fn;

      proj = geotiffinfo(geotiff_fn);

      [points.x,points.y] = projfwd(proj,points.lat,points.lon);

      if 0
        %% Debug plot
        clf;
        good_mask = isfinite(points.elev);
        scatter(points.lon(good_mask),points.lat(good_mask),[],points.elev(good_mask),'Marker','.')
        hcolor = colorbar;
        set(get(hcolor,'YLabel'),'String','Elevation (WGS-84,m)')
        hold on;
        hplot = plot(points.lon(33,:),points.lat(33,:),'k');
        xlabel('Longitude (deg,E)');
        ylabel('Latitude (deg,N)');
        legend(hplot,'Flight line');
      end

      %% Grid data

      min_x = min(points.x(:));
      max_x = max(points.x(:));
      min_y = min(points.y(:));
      max_y = max(points.y(:));

      if 1
        [DEM, R, tmp] = geotiffread(geotiff_fn);
        figure(1); clf;
        h_img = imagesc( (R(3,1) + R(2,1)*(0:size(DEM,2)-1))/1e3, (R(3,2) + R(1,2)*(0:size(DEM,1)-1))/1e3, DEM);
        set(gca,'YDir','normal');
        hold on;
        good_mask = isfinite(points.elev);
        scatter(points.x(good_mask)/1e3,points.y(good_mask)/1e3,[],points.elev(good_mask),'Marker','.')
        hcolor = colorbar;

        surf_str = lower(surf(active_surf).name);
        idx = regexp([' ' surf_str],'(?<=\s+)\S','start')-1;
        surf_str(idx) = upper(surf_str(idx));
        if ~strncmp(surf_str,'Ice',3)
          surf_str = ['Ice ',surf_str];
        end
        set(get(hcolor,'YLabel'),'String',sprintf('%s Elevation (WGS-84,m)',surf_str))

        hplot = plot(points.x(33,:)/1e3,points.y(33,:)/1e3,'k');
        xlabel('X (km)');
        ylabel('Y (km)');
        legend(hplot,'Flight line');
        axis([min_x max_x min_y max_y]/1e3);
        clip_and_resample_image(h_img,gca,10);
        surf_str = strrep(surf(active_surf).name,' ','_');
        if ~strncmp(surf_str,'ice',3)
          surf_str = ['ice_',surf_str];
        end
        out_fn_name = sprintf('%s_%03d_%s',param.day_seg,frm,surf_str);
        saveas(1,[fullfile(out_dir,out_fn_name),'.fig']);
        saveas(1,[fullfile(out_dir,out_fn_name),'.jpg']);
      end

      % Griddata the result
      grid_spacing = 25;
      xaxis = min_x:grid_spacing:max_x;
      yaxis = (min_y:grid_spacing:max_y).';
      [xmesh,ymesh] = meshgrid(xaxis,yaxis);

      % Create a constrained delaunay triangulization that forces edges
      % along the boundary (concave_hull) of our swath
      % First side
      good_mask = isfinite(points.x) & isfinite(points.y);
      row = find(good_mask,1);
      col = floor(row/size(good_mask,1)) + 1;
      row = row - (col-1)*size(good_mask,1);
      B = bwtraceboundary(good_mask,[row col],'S',8,inf,'counterclockwise');
      B = B(:,1) + (B(:,2)-1)*size(good_mask,1);
      concave_hull = [B(1:end-1) B(2:end)];
      good_idxs = find(good_mask);
      new_good_idxs = 1:length(good_idxs);
      idx_translate = zeros(size(good_mask));
      idx_translate(good_idxs) = new_good_idxs;
      concave_hull = idx_translate(concave_hull);

      pnts = zeros(3,length(good_idxs));
      pnts(1,:) = points.x(good_idxs);
      pnts(2,:) = points.y(good_idxs);
      pnts(3,:) = points.elev(good_idxs);
      gps_time = repmat(mdata.GPS_time,size(good_mask,1),1);
      gps_time = gps_time(good_idxs);
      
      px = pnts(1,concave_hull(:,1));
      py = pnts(2,concave_hull(:,1));
      c1 = concave_hull(:,1);
      [xi,~,segs] = selfintersect(px,py);
      while ~isempty(xi)
        fprintf('Removing intersections from constraint polygon...\n');
        px(segs(:,1):segs(:,2)) = NaN;
        py(segs(:,2):segs(:,2)) = NaN;
        mask = ~isnan(px)&~isnan(py);
        px = px(mask);
        py = py(mask);
        c1 = c1(mask);
        [xi,~,segs] = selfintersect(px,py);
      end
      fprintf('Intersections removed.\n');
      concave_hull = zeros(length(c1),2);
      concave_hull(:,1) = c1;
      concave_hull(:,2) = [c1(2:end);c1(1)];
      

      % concave_hull = [ [(1:size(pnts,2)-1).' (2:size(pnts,2)).'];
      %   size(pnts,2)*[(1:size(pnts,3)-1).' (2:size(pnts,3)).'];
      %   size(pnts,2)*(size(pnts,3)-1) + [(size(pnts,2):-1:2).' (size(pnts,2)-1:-1:1).'];
      %   size(pnts,2)*[(size(pnts,3):-1:2).' (size(pnts,3)-1:-1:1).']-size(pnts,2)+1 ];
      dt = DelaunayTri(pnts(1,:).',pnts(2,:).',concave_hull);

      F = TriScatteredInterp(dt,pnts(3,:).');
      %F = TriScatteredInterp(dt,pnts(3,:).','natural');
      DEM = F(xmesh,ymesh);

      fprintf('  Creating boundary and removing outside points (%.1f sec)\n', toc);
      if 0
        % Slow method using inpolygon
        boundary = [pnts(1:2,:,1) squeeze(pnts(1:2,end,:)) squeeze(pnts(1:2,end:-1:1,end)) squeeze(pnts(1:2,1,end:-1:1))];
        good_mask = ~isnan(meas{measInd}.DEM);
        good_mask_idx = find(good_mask);
        in = inpolygon(eastMesh(good_mask),northMesh(good_mask),boundary(1,:),boundary(2,:));
        good_mask(good_mask_idx(~in)) = 0;
        meas{measInd}.DEM(~good_mask) = NaN;
      else
        % Faster method using inOutStatus/pointLocation with edge constraints
        % Finds the triangles outside the concave hull (the bad ones)
        bad_tri_mask = ~inOutStatus(dt);
        % Selects just the points that were in the convex hull
        good_mask = ~isnan(DEM);
        good_mask_idx = find(good_mask);
        % For each point in the convex hull, find the triangle enclosing it
        tri_list = pointLocation(dt,xmesh(good_mask),ymesh(good_mask));
        % Use the bad triangle list to find the bad points
        bad_mask = bad_tri_mask(tri_list);
        % Set all the points in bad triangles to NaN
        good_mask(good_mask_idx(bad_mask)) = 0;
        DEM(~good_mask) = NaN;
      end

      % DEM surface
      % phong value for the FaceLighting and EdgeLighting
      % =========================================================================
      % =========================================================================
      if 0
        figure(1); clf;
        set(1,'Color',[1 1 1]);

        hA2 = axes;
        %hC = contourf((xaxis-xaxis(1))/1e3,(yaxis-yaxis(1))/1e3,double(DEM),12);
        clear surf;
        hC = surf((xaxis-xaxis(1))/1e3,(yaxis-yaxis(1))/1e3,double(DEM)*0+25,double(DEM));
        set(hC(1),'EdgeAlpha',0); grid off;
        %axis([2 38 0 10 zlims]);
        set(hA2,'Box','off');
        set(hA2,'View',[7.5   74]);
        set(hA2,'Position',[0.025 0.1 0.875 0.53]);
        set(hA2,'YTick',[0 5 10]);
        set(hA2,'ZColor',[1 1 1]);
        colormap(jet(256))
        hc = colorbar;
        zlims = caxis;
        caxis(zlims);
        set(hc,'Position',[0.95 0.1 0.015 0.8]);
        set(get(hc,'YLabel'),'String','Bed height (m)');

        hA = axes;
        hS = surf((xaxis-xaxis(1))/1e3,(yaxis-yaxis(1))/1e3,double(DEM),double(1*DEM));
        hA = gca; grid off;
        %axis([2 38 0 10 zlims]);
        set(hA,'Box','off');
        % View ["" "elevation angle of observer"]
        set(hA,'View',[7.5 74]);
        set(hA,'Position',[0.025 0.45 0.875 0.66]);
        set(hA,'XTick',[]);
        set(hA,'YTick',[]);
        set(hA,'ZTick',[]);
        set(hA,'XColor',[1 1 1]);
        set(hA,'YColor',[1 1 1]);
        set(hA,'ZColor',[1 1 1]);
        set(hS(1),'EdgeAlpha',0.2);

        if 0 % surf with color
          colormap(jet(256))
          hc = colorbar;
          caxis(zlims);
          set(hc,'Position',[0.95 0.1 0.015 0.8]);
          set(get(hc,'YLabel'),'String','Bed height (m)');
        end

        axes(hA2);
        hx = xlabel('X (km)');
        %set(hx,'Position',[19 -1.8 0]);
        hy = ylabel('Y (km)');
        %set(hy,'Position',[40 3.5 0]);
        %set(hy,'Rotation',54);

        % =========================================================================
        % =========================================================================
        figure(3); clf;
        set(3,'Position',[50 50 500 500]);
        set(3,'Color',[1 1 1]);

        hS = surf((xaxis-xaxis(1))/1e3-2,(yaxis-yaxis(1))/1e3-1,double(DEM),double(DEM));
        hA = gca; grid off;
        %axis([0 36 0 10 zlims]);
        set(hA,'Box','off');
        % View ["" "elevation angle of observer"]
        set(hA,'YTick',[0 5 10]);
        set(hA,'ZColor',[1 1 1]);
        set(hA,'View',[7.5 74]);
        set(hA,'Position',[0.025 0.135 0.875 1.1]);
        set(hA,'ZColor',[1 1 1]);
        set(hS(1),'EdgeAlpha',0.2);

        colormap(jet(256))
        hc = colorbar;
        set(hc,'Position',[0.95 0.1 0.015 0.8]);
        set(get(hc,'YLabel'),'String','Bed height (m)');
        caxis(zlims);

        hx = xlabel('X (km)');
        set(hx,'Position',[19 -1.2 0]);
        set(hx,'Rotation',-2);
        hy = ylabel('Y (km)');
        set(hy,'Position',[38 3.5 0]);
        set(hy,'Rotation',68);
      end

      if 1
        % =========================================================================
        % =========================================================================
        zlims = [min(DEM(:)) max(DEM(:))];
        DEM_contour = DEM;
        DEM_contour(isnan(DEM_contour)) = zlims(1)-100;
        figure(2); clf
        contourf((xaxis-xaxis(1))/1e3,(yaxis-yaxis(1))/1e3,DEM_contour,linspace(zlims(1),zlims(end),10));
        % axis([22 32 1.5 9]);
        hx = xlabel('X (km)');
        hy = ylabel('Y (km)');
        colormap(jet(256));
        hcolor = colorbar;
        surf_str = lower(surf(active_surf).name);
        idx = regexp([' ' surf_str],'(?<=\s+)\S','start')-1;
        surf_str(idx) = upper(surf_str(idx));
        if ~strncmp(surf_str,'Ice',3)
          surf_str = ['Ice ',surf_str];
        end
        set(get(hcolor,'YLabel'),'String',sprintf('%s Elevation (WGS-84,m)',surf_str))
        surf_str = strrep(surf(active_surf).name,' ','_');
        if ~strncmp(surf_str,'ice',3)
          surf_str = ['ice_',surf_str];
        end
        out_fn_name = sprintf('%s_%03d_%s_contour',param.day_seg,frm,surf_str);
        saveas(2,[fullfile(out_dir,out_fn_name),'.fig']);
        saveas(2,[fullfile(out_dir,out_fn_name),'.jpg']);
      end


      %% Save Geotiff of Surface


      %% Save Point File of Surface


      %% Save Geotiff of Bottom

      % ProjectedCSTypeGeoKey: The projection type key to use with geotiff_en.
      % The most reliable way to
      % set this is to create a geotiff file with the right projection in
      % another program and read in the key from that program using geotiffinfo.
      % proj = geotiffinfo('/cresis/projects/GIS_data/greenland/Landsat-7/Greenland_natural_150m.tif');
      % proj.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey
      %  ans =
      %       32767
      % Some keys are not fully supported and you'll just have to try different
      % projections until you find one that works.  Test the output file by
      % loading it with Matlab and another program like Arc or Globalmapper
      ProjectedCSTypeGeoKey = 32622;
      ProjectedCSTypeGeoKey = 3031;
      ProjectedCSTypeGeoKey = proj.GeoTIFFTags.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey;

      R = maprasterref;
      R.XLimWorld = [xaxis(1) xaxis(end)];
      R.YLimWorld = [yaxis(1) yaxis(end)];
      R.RasterSize = size(DEM);
      %R.RasterInterpretation = 'cells';
      R.ColumnsStartFrom = 'south';
      R.RowsStartFrom = 'west';
      %R.RasterInterpretation = 'postings';

      key.GTModelTypeGeoKey  = 1;  % Projected Coordinate System (PCS)
      key.GTRasterTypeGeoKey = 2;  % PixelIsPoint
      key.GTRasterTypeGeoKey = 1;  % PixelIsPoint
      key.ProjectedCSTypeGeoKey = ProjectedCSTypeGeoKey;
      DEM(isnan(DEM)) = -9999;

      DEM_geotiff_dir = ct_filename_out(param,'DEM');
      surf_str = strrep(surf(active_surf).name,' ','_');
      if ~strncmp(surf_str,'ice',3)
        surf_str = ['ice_',surf_str];
      end
      DEM_geotiff_fn_name = sprintf('%s_%03d_%s_DEM',param.day_seg,frm,surf_str);
      DEM_geotiff_fn = fullfile(DEM_geotiff_dir,DEM_geotiff_fn_name);
      fprintf('Creating output %s\n', DEM_geotiff_fn);
      geotiffwrite(DEM_geotiff_fn, int16(DEM), R, 'GeoKeyDirectoryTag', key);


      %% Save Point File of Bottom

      [pnts(4,:),pnts(5,:)] = projinv(proj,pnts(1,:),pnts(2,:));
      pnts(6,:) = gps_time;

      csv_dir = ct_filename_out(param,'DEM');
      surf_str = strrep(surf(active_surf).name,' ','_');
      if ~strncmp(surf_str,'ice',3)
        surf_str = ['ice_',surf_str];
      end
      csv_fn_name = sprintf('%s_%03d_%s.csv',param.day_seg,frm,surf_str);
      csv_fn = fullfile(csv_dir,csv_fn_name);
      fid = fopen(csv_fn,'w');
      fprintf(fid,'X,Y,Elevation,Latitude,Longitude,GPS Time\n');
      fprintf(fid,'%f,%f,%f,%f,%f,%f\n',pnts);
      fclose(fid);

      mat_dir = ct_filename_out(param,'DEM');
      surf_str = strrep(surf(active_surf).name,' ','_');
      if ~strncmp(surf_str,'ice',3)
        surf_str = ['ice_',surf_str];
      end
      mat_fn_name = sprintf('%s_%03d_%s',param.day_seg,frm,surf_str);
      mat_fn = fullfile(mat_dir,mat_fn_name);
      
      sw_version = [];
      if isfield(param,'sw_version')
        sw_version = param.sw_version;
      end
      param_combine = [];
      if isfield(mdata,'param_combine')
        param_combine = mdata.param_combine;
      end
      ice_mask_ref = [];
      if isfield(param,'ice_mask_ref')
        ice_mask_ref = param.ice_mask_ref;
      end
      geotiff_ref = [];
      if isfield(param,'geotiff_ref')
        geotiff_ref = param.geotiff_ref;
      end
      DEM_ref = [];
      if isfield(param,'DEM_ref')
        DEM_ref = param.DEM_ref;
      end

      save(mat_fn,'sw_version','param_combine','ice_mask_ref','geotiff_ref','DEM_ref');


      if 0
        %% Create a 3-D surface image (fancy looking figure)
        fig_h = 20;
        figure(fig_h); clf;
        hA2 = axes;
        hC = surf((xaxis-xaxis(1))/1e3,(yaxis-yaxis(1))/1e3,double(medfilt2(DEM,[3 11])));
        set(hC,'EdgeAlpha',0.1); grid off;
        set(hA2,'View',[190 80]);
        hold on;
        % h = plot3((neem.east-xaxis(1))/1e3,(neem.north-yaxis(1))/1e3,neem.DEM,'ko');
        % set(h, 'LineWidth', 3);
        % set(h, 'MarkerSize', 7);
        hold off;
        %   hx = xlabel('X (km)');
        %   set(hx,'Position',[30 15 0]);
        %   set(hx,'Rotation',-3);
        %   hy = ylabel('Y (km)');
        %   set(hy,'Position',[-3 9 -1500]);
        %   set(hy,'Rotation',57);
        zlabel(sprintf('WGS-84 bed\nelevation (m)'));
        grid on;
        hC = colorbar;
        colormap(jet(256))
        set(get(hC,'YLabel'),'String','WGS-84 bed elevation (m)');
        set(hA2,'Position',[0.12 0.11 0.72 0.815])
        set(hC,'Position',[0.9 0.11 0.022 0.815])

        set(fig_h,'PaperOrientation','landscape');
        set(fig_h,'PaperPosition',[0.5 0.5 10 3]);

        %   print(sprintf('-f%d',fig_h),'-djpeg','-r200','C:\Users\radar\Desktop\NEEM_NGRIP_3D_all.jpg');

      end
    end
  end
end