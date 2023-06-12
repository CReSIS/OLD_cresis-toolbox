function [qq] = qloox(xo, ident)

figures_plot = 1;
figures_visible = 0;

global gRadar;
gRadar.verbose_off = {'read_param_xls_radar','read_param_xls_generic'};
[c, er_ice] = physical_constants('c','er_ice');

xo_table_tag = extractBefore(ident,'_xo_');
idx_xo = str2double(extractAfter(ident,'_xo_'));
reuse_loc = fullfile(gRadar.ct_tmp_path, 'radcal_mat', xo_table_tag);
xo_hdr = sprintf('(xo %d) ', idx_xo);
fig_hdr = xo_hdr;

%% Load(param,records,frames)+find(rec,frm)+check(xo)

if figures_plot == 1

  h_fig_instagram = figure('Name','instagram', 'visible', figures_visible);
  hold on;

  h_fig_geoplot = figure('Name','geoplot', 'visible', figures_visible);
  h_geo = geoaxes(h_fig_geoplot, 'Basemap','darkwater');

  h_fig_records = figure('Name','records', 'visible', figures_visible);
  hold on;
  %   plot(xo.cx_lon, xo.cx_lat, '*', 'LineWidth', 4);

  h_fig_waveform_elev = figure('Name','waveform', 'visible', figures_visible);
  hold on;
  h_fig_waveform_time = figure('Name','waveform2', 'visible', figures_visible);
  hold on;
end

for xx=1:2

  % param
  param{xx} = read_param_xls( ...
    fullfile(gRadar.param_path, sprintf('rds_param_%s.xls', ...
    eval(sprintf('xo.pt%d_season{:}',xx))) ), ...
    eval(sprintf('xo.pt%d_segment',xx)), ...
    [] );

  fprintf('(%s) %d >> %s %s ', datestr(now), ...
    xx, param{xx}.season_name, param{xx}.day_seg);

  % records
  try
    fprintf('>>records');
    records{xx} = records_load( param{xx} );
  catch ME
    fprintf('--Failed. '); return;
  end
  fprintf('--Loaded. ');
  % check the case of prev_rec, xo, next_&_closest_rec
  % check XUNZY to match with xo
  % is_xo_next = xo.

  try
    fprintf('>>frames');
    frames{xx} = frames_load( param{xx} );
  catch ME
    fprintf('--Failed. '); return;
  end
  fprintf('--Loaded. \n');

  %% find(rec+frm)

  dx_rec(xx) = mean(distance_geodetic(records{xx}));

  rr{xx} = eval( sprintf('abs(records{xx}.gps_time - xo.pt%d_gps_time);', xx) );
  [rr_val(xx), rec(xx)] = min(rr{xx},[],2);

  record{xx} = struct_select(records{xx}, length(records{xx}.gps_time), rec(xx), 0);

  frm(xx) = find(rec(xx) >= frames{xx}.frame_idxs,1,'last' );
  if ~strcmp(eval(sprintf('xo.pt%d_frame{:}',xx')), ...
      sprintf('%s_%03d', param{xx}.day_seg, frm(xx)))
    warning('frm mismatch !!!\n');
    return;
  end

  local_hdr{xx} = sprintf('[%s %s_%03d] ', param{xx}.season_name, param{xx}.day_seg, frm(xx));

  %% Load data for frm + find(data_rec)

  dd_full = [];
  try
    fn = fullfile(ct_filename_out(param{xx},'qlook'), ...
      sprintf('Data_%s_%03d.mat', param{xx}.day_seg, ...
      frm(xx)));
    dd_full = load(fn);
  catch
    fn = fullfile(ct_filename_out(param{xx},'qlook'), ...
      sprintf('Data_img_01_%s_%03d.mat', param{xx}.day_seg, ...
      frm(xx)));
    dd_full = load(fn);
  end
  %   dd_power_err = analyze_power_levels(dd_full, sprintf('%s (point %d) %s', xo_hdr ,xx, local_hdr) );

  dx_dd(xx) = mean(distance_geodetic(dd_full));

  dd_rr{xx} = eval( sprintf('abs(dd_full.GPS_time - xo.pt%d_gps_time);', xx) );
  [dd_rr_val(xx), rec_dd(xx)] = min(dd_rr{xx},[],2);

  %% insta_gram
  if figures_plot == 1

    ig_bounds_each_side = 10;
    ig_bounds_max = size(dd_full.Data,2);

    idx_ig_min = max([rec_dd(xx)-ig_bounds_each_side, 1]);
    idx_ig_max = min([rec_dd(xx)+ig_bounds_each_side, ig_bounds_max]);
    idxs_ig = idx_ig_min:idx_ig_max;
    tmp = 10*log10(dd_full.Data(:,idxs_ig));

    figure(h_fig_instagram);
    subplot(2,1,xx);
    imagesc(idxs_ig, dd_full.Time/1e-6, tmp);
    cb = colorbar; cb.Label.String = 'Relative Power, dB';
    hold on;
    xline(rec_dd(xx), 'r-');
    axis tight;
    xlabel('rlines'); ylabel('FastTime, us');
    title(sprintf('partial radargram  %s', local_hdr{xx}), 'Interpreter', 'none');

    clear ig_bounds_each_side ig_bounds_max idx_ig_min idx_ig_max idxs_ig
    clear tmp cb
  end

  %% struct_select rec_dd from dd

  dd = struct_select(dd_full, length(dd_full.GPS_time), rec_dd(xx), 0);

  %% Truncate data

  if ndims(dd.Data)>2; continue; end
  [Nt, Nx, ~] = size(dd.Data);
  if Nt<2; continue; end % to support air_idxs, sub_idxs
  if Nx~=1; continue; end % support only one waveform

  dd_dt = dd.Time(2)-dd.Time(1);
  dd_dt_dist = dd_dt * c/2;

  tmp = 10*log10(dd.Data);
  % [tmp_rline_max_vals, tmp_rline_max_idxs] = max(tmp,[], 1, 'omitnan');
  % max is not always the surface

  % compare (bad way)
  [I,~,~] = find_multiple(dd.Time, '>=', dd.Surface, 1, 'first');
  surf_next_idx = [I{:}];

  % distance or closest (better way)
  bb = abs(dd.Time - dd.Surface);
  [~, surf_closest_idx] = min(bb,[],1);


  if 0
    figure('Name', 'debug');
    plot(dd.Surface/1e-6,'.'); hold on; grid on;
    % plot(dd.Time(tmp_rline_max_idxs)/1e-6,'o'); % dont
    plot(dd.Time(surf_next_idx)/1e-6,'rs'); % bad
    plot(dd.Time(surf_closest_idx)/1e-6,'go'); % good
    XMinorGrid ='on'; yticks( dd.Time(min(surf_closest_idx):max(surf_closest_idx))/1e-6 );
    legend('Surface','Next', 'Closest');
    xlabel('rlines'); ylabel('us');
  end

  clear I surf_next_idx bb

  %% detect surface
  % because surf_closest_idxs is not accurate at all times (tracked)

  if 1
    box_bounds_each_side = round(50/dd_dt_dist); % within 50 meter
    idxs_box = boxing_1D(surf_closest_idx, box_bounds_each_side, Nt);
    [~, box_peak_idx] = max(tmp(idxs_box));
    surf_peak_idx = idxs_box(1) + box_peak_idx -1;

    surface_box(xx,:) = idxs_box;
    clear box_bounds_each_side box_bounds_max idx_box_min idx_box_max idxs_box box_peak_idx
  end

  % even this peak in the box could be something other than surface
  % use findpeaks next time with prominence of 20-23 dB
  % pick the first one close to the surf_closest_idxs
  % deal with feedthrough??

  %% adjust surf_idx

  %tt = dd.Time - dd.Time(1);

  if ( ~exist('surf_peak_idx', 'var') || isempty(surf_peak_idx) )
    surf_idx(xx) = surf_closest_idx; %%%#########################
  else
    surf_idx(xx) = surf_peak_idx;
  end

  if xx==2 %coreg second onto first wf
    a1 = ddd{1}.tmp(surface_box(1,:));
    a2 = tmp(surface_box(2,:));

    % Number of indixes to move a2 wrt a1
    lagaan = coreg_1D(a1,a2, ident);

    if isnan(lagaan) || lagaan == 0
      % no new flag necessary
      % lagaan itself is aflag if not zero or NaN
      % do nothing
    else
      % adjust surf_idx
      surf_idx(xx) = surf_idx(xx) + lagaan;
    end
  end

  %% Transform time axes to WGS-84 elevation

  % create the elevation_axis
  [dd.elev_axis, surface_calc(xx), AGL(xx)] = time2elev(dd.Time, dd.Elevation, surf_idx(xx));

  % create a time axis to match surfaces
  dd.time_axis = dd.Time - dd.Time(surf_idx(xx));

  % Calc powers
  P_est(xx) = 10*log10(1/(4*pi*AGL(xx).^2));
  P_act(xx) = tmp(surf_idx(xx));
  P_err(xx) = P_est(xx) - P_act(xx);

  %% figures
  if figures_plot == 1
    % geoplot
    figure(h_fig_geoplot);
    geoplot(h_geo, records{xx}.lat, records{xx}.lon, ':'); hold on;

    % records
    figure(h_fig_records);

    % plot only a few records on either side of xo
    box_bounds_each_side = max(1, round(500/dx_rec(xx))); % within 500 meter
    plot_rec_idxs = boxing_1D(rec(xx), box_bounds_each_side, length(records{xx}.lon));
    %plot(records{xx}.lon, records{xx}.lat, ':');
    plot3(records{xx}.lon(plot_rec_idxs), records{xx}.lat(plot_rec_idxs), records{xx}.elev(plot_rec_idxs), ':');
    plot3(record{xx}.lon, record{xx}.lat, record{xx}.elev, 'o', 'LineWidth', 1);

    % plot only a few dd_records on either side of xo
    box_bounds_each_side = max(1, round(500/dx_dd(xx))); % within 500 meter
    plot_rec_idxs = boxing_1D(rec_dd(xx), box_bounds_each_side, length(dd_full.Latitude));
    %plot(dd_full.Longitude, dd_full.Latitude, '+', 'LineWidth', 2);
    plot3(dd_full.Longitude(plot_rec_idxs), dd_full.Latitude(plot_rec_idxs), dd_full.Elevation(plot_rec_idxs), '+', 'LineWidth', 2);
    % plot(dd.Longitude, dd.Latitude, 's', 'LineWidth', 1);
    plot3(dd.Longitude, dd.Latitude, dd.Elevation, 's', 'LineWidth', 1);

    % pos_vect
    line([1;1].*dd.Longitude, [1;1].*dd.Latitude, ...
      [dd.Elevation; surface_calc(xx)], ...
      'LineStyle','-', 'LineWidth',1, ...
      'Color',[152,251,152]/256);
    if xx==1
      plot3(dd.Longitude, dd.Latitude, surface_calc(xx), 'o','LineWidth',2,'Color','r');
    else
      plot3(dd.Longitude, dd.Latitude, surface_calc(xx), 'o','LineWidth',2,'Color','b');
    end

    figure(h_fig_waveform_elev);
    plot(dd.elev_axis, tmp,'.-');
    if xx==1
      vert_align = 'top';
    else
      vert_align = 'bottom';
    end
    xline(dd.Elevation, ':', ...
      sprintf('[xo\\_%d] Elev. %0.2f m', xx, dd.Elevation), ...
      'LabelVerticalAlignment', vert_align, 'LabelOrientation', 'aligned');
    xline(surface_calc(xx), ':', ...
      sprintf('[xo\\_%d] Calc. surf %0.2f m', xx, surface_calc(xx)), ...
      'LabelVerticalAlignment', vert_align, 'LabelOrientation', 'aligned');

    figure(h_fig_waveform_time);
    subplot(121); hold on;
    plot(dd.time_axis/1e-6, tmp,'.-');
    subplot(122); hold on;
    plot(dd.time_axis/1e-6, tmp,'.-');

  end

  % for xo table
  dd.tmp = tmp; % log scale
  ddd{xx} = dd;
end

fprintf('\n');

%% figures closing
%set(h_fig_waveform_time, 'Position', get(0, 'Screensize'));

fig_hdr = [fig_hdr, local_hdr{1} local_hdr{2}];

if figures_plot == 1

  figure(h_fig_instagram);
  fig_fn = fullfile(reuse_loc, sprintf('instagram_%s.fig', ident));
  fprintf('Saving %s\n', fig_fn);
  ct_saveas(h_fig_instagram,fig_fn);
  fig_fn = fullfile(reuse_loc, sprintf('instagram_%s.png', ident));
  fprintf('Saving %s\n', fig_fn);
  ct_saveas(h_fig_instagram,fig_fn);

  figure(h_fig_geoplot);
  title(fig_hdr, 'Interpreter', 'None');
  fig_fn = fullfile(reuse_loc, sprintf('geoplot_%s.fig', ident));
  fprintf('Saving %s\n', fig_fn);
  ct_saveas(h_fig_geoplot,fig_fn);
  fig_fn = fullfile(reuse_loc, sprintf('geoplot_%s.png', ident));
  fprintf('Saving %s\n', fig_fn);
  ct_saveas(h_fig_geoplot,fig_fn);


  figure(h_fig_records);

  grid on;
  xlabel('Longitude, deg');
  ylabel('Latitude, deg');
  zlabel('WGS-84 Elevation, m');
  rotate3d on;
  view([45,29]);
  h_tmp=gca;
  line([1 1].*xo.cx_lon(1), [1 1].*xo.cx_lat(1), h_tmp.ZLim, ...
    'LineStyle', ':', 'Color', 'k', 'LineWidth', 2);
  legend(...
    'records1','xo rec1','dd recs1', 'xo dd rec1', 'pos_vect1', 'surf1', ...
    'records2','xo rec2','dd recs2', 'xo dd rec2', 'pos_vect2', 'surf2', ...
    'xo');
  title(fig_hdr, 'Interpreter', 'None');
  fig_fn = fullfile(reuse_loc, sprintf('records_%s.fig', ident));
  set(findobj(h_fig_records,'type','axes'),'FontWeight', 'Bold', 'FontSize',12);
  set(h_fig_records, 'Position', get(0, 'Screensize'));
  fprintf('Saving %s\n', fig_fn);
  ct_saveas(h_fig_records,fig_fn);
  fig_fn = fullfile(reuse_loc, sprintf('records_%s.png', ident));
  fprintf('Saving %s\n', fig_fn);
  ct_saveas(h_fig_records,fig_fn);

  figure(h_fig_waveform_elev);
  grid on;
  set(gca, 'Xdir', 'reverse');
  xlabel('elev\_axis, meter');
  ylabel('Magnitude, dB');
  title(fig_hdr, 'Interpreter', 'None');
  fig_fn = fullfile(reuse_loc, sprintf('waveform_elev_%s.fig', ident));
  set(findobj(h_fig_waveform_elev,'type','axes'),'FontWeight', 'Bold', 'FontSize',12);
  set(h_fig_waveform_elev, 'Position', get(0, 'Screensize'));
  fprintf('Saving %s\n', fig_fn);
  ct_saveas(h_fig_waveform_elev,fig_fn);
  fig_fn = fullfile(reuse_loc, sprintf('waveform_elev_%s.png', ident));
  fprintf('Saving %s\n', fig_fn);
  ct_saveas(h_fig_waveform_elev,fig_fn);

  figure(h_fig_waveform_time);
  subplot(121)
  grid on;
  %set(gca, 'Xdir', 'reverse');
  axis tight;
  % ylim([min(P_act)-100 max(P_act)+10]);
  ylim([-180 -40]);
  xlabel('surface-zeroed Fast-time, us');
  ylabel('Magnitude, dB');
  legend(local_hdr, 'Interpreter','none');
  title(sprintf('Est-Act=Err [1]%0.2f-%0.2f=%0.2f] [2]%0.2f-%0.2f=%0.2f', ...
    P_est(1), P_act(1), P_err(1), ...
    P_est(2), P_act(2), P_err(2)))

  subplot(122)
  grid on;
  axis([-3 9 -110 -30]);
  xlabel('surface-zeroed Fast-time, us');
  ylabel('Magnitude, dB');
  legend(local_hdr, 'Interpreter','none');
  title(sprintf('diff(1,2) [Est, Act, Err] = [%0.2f, %0.2f, %0.2f]', ...
    diff(P_est), diff(P_act), diff(P_err) ));

  sgtitle(sprintf('(xo %d) Power levels: Estimated, Actual, Error', idx_xo), 'Interpreter', 'None');

  fig_fn = fullfile(reuse_loc, sprintf('waveform_time_%s.fig', ident));
  set(findobj(h_fig_waveform_time,'type','axes'),'FontWeight', 'Bold', 'FontSize',12);
  set(h_fig_waveform_time, 'Position', get(0, 'Screensize'));
  fprintf('Saving %s\n', fig_fn);
  ct_saveas(h_fig_waveform_time,fig_fn);
  fig_fn = fullfile(reuse_loc, sprintf('waveform_time_%s.png', ident));
  fprintf('Saving %s\n', fig_fn);
  ct_saveas(h_fig_waveform_time,fig_fn);



end

%% xo table
qq = struct();
qq.dx_rec = dx_rec;
qq.dx_dd = dx_dd;

qq.rec = rec;
qq.rec_dd = rec_dd;
qq.frm = frm;

qq.laggan = lagaan;
qq.surface_box = surface_box;
qq.surf_idx = surf_idx;
qq.AGL = AGL;
qq.surface_calc = surface_calc;

qq.P_act = P_act;
qq.P_est = P_est;
qq.P_err = P_err;

qq.record = record;
qq.dd = ddd;

qq.ident = ident;

%% Find the matching records

%% crossover_view

%% close figs
try
  close(h_fig_instagram);
  close(h_fig_geoplot);
  close(h_fig_records);
  close(h_fig_waveform_elev);
  close(h_fig_waveform_time);
end