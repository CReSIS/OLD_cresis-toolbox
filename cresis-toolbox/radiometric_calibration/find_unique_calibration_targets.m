%% OM

% Script find_unique_calibration_targets.m
%
% Top level script that finds:
% 1. natural targets
% 2. target scenes with close flightpaths and crossovers
%
% Authors: John Paden, Hara Madhav Talasila

%% todo

% take two segments for starters (different)
% load records, frames,
% find crossovers (use lat, lon, elev)
%   --> use distance_geodetic.m
%   --> create xunzy and see angle of crossover, z vectors
% find speculars (look in deconv)

hara;

% =========================================================================
fprintf('=====================================================================\n');
fprintf('%s: (%s)\n', mfilename, datestr(now));
fprintf('=====================================================================\n');

%% Load params
% PICK season(s), day_seg(s)

params      = [];
N_seasons   = [];
N_day_segs  = [];


if 0
  % MULTIPLE day_segs, MULTIPLE seasons
  N_day_segs = 2; % for now
  for idx_day_seg = 1: N_day_segs
    params{idx_day_seg} = param_gui;
    % work on using the cell params later
  end
  
elseif 0
  params{end+1} = read_param_xls('/users/htalasila/scripts/ct_params/rds_param_2018_Greenland_P3.xls',{'20180423_01'},[]);
  params{end+1} = read_param_xls('/users/htalasila/scripts/ct_params/rds_param_2014_Greenland_P3.xls',{'20140502_01'},[]);
  
elseif 0
  % SINGLE day_seg, SINGLE season
  fn_param = '/users/htalasila/scripts/ct_params/rds_param_2014_Greenland_P3.xls';
  day_segs = {'20140401_01'};
  sheets = [];
  params = read_param_xls(fn_param,day_segs,sheets);
  
elseif 1
  % MULTIPLE day_segs, SINGLE season
  fn_param = '/users/htalasila/scripts/ct_params/rds_param_2014_Greenland_P3.xls';
  day_segs = {'20140409_01', '20140409_02', '20140419_02', '20140419_03'};
  day_segs = {'20140409_01', '20140409_02'};
  day_segs = {'20140419_02', '20140419_03'};
  sheets = [];
  %   day_segs = '20140419_02|20140419_03';
  params = read_param_xls(fn_param,day_segs,sheets);
  
  
else
  % generic case...
  params = param_gui;
  
end

%% Setup params

if iscell(params)
  % multiple seasons
else
  params = {params};
end

N_seasons = length(params);
N_day_segs  = [];
for idx_season = 1:N_seasons
  N_day_segs{idx_season} = length(params{idx_season});
end

%% Load records, frames
h_gx = geoaxes('Basemap','darkwater');

for idx_season = 1:N_seasons
  for idx_day_seg = 1:N_day_segs{idx_season}
    
    param = params{idx_season}(idx_day_seg);
    fprintf('\n(%s) season (%d of %d) day_seg (%d of %d) >> %s %s ', datestr(now), ...
      idx_season, N_seasons, idx_day_seg, N_day_segs{idx_season}, ...
      param.season_name, param.day_seg);
    
    try
      fprintf('>>records');
      records{idx_season, idx_day_seg} = records_load( param );
      geoplot(h_gx, records{idx_season, idx_day_seg}.lat, records{idx_season, idx_day_seg}.lon); hold on;
      figure(222)
      plot(records{idx_season, idx_day_seg}.lon, records{idx_season, idx_day_seg}.lat, '.-'); hold on;
    catch ME
      fprintf('--Failed. '); return;
    end
    fprintf('--Loaded. ');
    
    try
      fprintf('>>frames',  datestr(now));
      frames{idx_season, idx_day_seg} = frames_load( param );
    catch ME
      fprintf('--Failed. '); return;
    end
    fprintf('--Loaded. ');
    
  end % day_segs
end % season

fprintf('\n');

%%
% How to find the intersection?

if 0
  % grid the entire greenland??
  for idx_season = 1:N_seasons
    for idx_day_seg = 1:N_day_segs{idx_season}
      d_geo{idx_season, idx_day_seg} = distance_geodetic(records{idx_season, idx_day_seg});
      d_geo_median{idx_season, idx_day_seg} = median(d_geo{idx_season, idx_day_seg})
      
    end
  end
end

% e
% return;

%% this is abandoned because John suggested a better approach
%% gen_crossover_database.m is the next one!!

%% meanwhile
% day_segs = {'20140419_02', '20140419_03'};
% cursor_info for latlon = 68.820138409998606 -49.543126087195560  68.820137918218123 -49.543124828631676
% DataIndex = 5172216      564847
xo_rec_db= { 564847+[0:1], 5172216+[0:1] };
%     xo_rec_db= { [5033608]+[0:1], [6546089]+[0:1] };
%     xo_rec_db= { [6546089]+[-1:0]; [5033608]+[0:1] };

for idx_season = 1:N_seasons
  for idx_day_seg = 1:N_day_segs{idx_season}
    xo_recs = xo_rec_db{idx_season, idx_day_seg};
    xo_frms = [];
    for idx_rec = 1:length(xo_recs)
      rec = xo_recs(idx_rec);
      frm = find(rec >= frames{idx_season, idx_day_seg}.frame_idxs,1,'last' );
      plot(records{idx_season, idx_day_seg}.lon(rec), records{idx_season, idx_day_seg}.lat(rec),'o');
      fprintf('%s %s t[%d] geo[%4.8f, %4.8f, %4.4f] frm = %03d\n', ...
        params{idx_season}(idx_day_seg).season_name, ...
        params{idx_season}(idx_day_seg).day_seg, ...
        records{idx_season, idx_day_seg}.gps_time(rec),...
        records{idx_season, idx_day_seg}.lat(rec),...
        records{idx_season, idx_day_seg}.lon(rec),...
        records{idx_season, idx_day_seg}.elev(rec),...
        frm);
      xo_frms = [xo_frms, frm];
    end
    xo_frm_db{idx_season, idx_day_seg} = unique(xo_frms);
    
    fprintf('dx= %f m \n\n', median(distance_geodetic(struct_select(records{idx_season, idx_day_seg}, length(records{idx_season, idx_day_seg}.gps_time), xo_recs,0) )));
  end
end

%%
%     qloox(params, records, frames, xo_rec_db, xo_frm_db)

for idx_season = 1:N_seasons
  for idx_day_seg = 1:N_day_segs{idx_season}
    fn = fullfile(ct_filename_out(params{idx_season}(idx_day_seg),'qlook'), ...
      sprintf('Data_%s_%03d.mat', params{idx_season}(idx_day_seg).day_seg, ...
      xo_frm_db{idx_season, idx_day_seg}));
    dd = load(fn);
    xo_recs = xo_rec_db{idx_season, idx_day_seg};
    xo_rec_gps_times  = records{idx_season, idx_day_seg}.gps_time(xo_recs);
    
    [~, J, ~] = find_multiple(dd.GPS_time, '>=', xo_rec_gps_times, 1, 'first');
    
    xo_dd_rline_idxs = unique([J{:}]);
    
    extra_picks = 0;
    if ~isempty(xo_dd_rline_idxs)
      xo_dd_rline_picks = max(1, min(xo_dd_rline_idxs)-extra_picks) : ...
        min(max(xo_dd_rline_idxs)+extra_picks, length(dd.GPS_time)) ;
    end
    
    dd = struct_select(dd,length(dd.GPS_time), xo_dd_rline_picks, 0);
    
    if ndims(dd.Data)>2
      continue;
    end
    
    [Nt, Nx, ~] = size(dd.Data);
    
    if Nt<2
      continue; % to support air_idxs, sub_idxs
    end
    
    dx_vec = distance_geodetic(dd);
    dt_ddTime = median(diff(dd.Time));
    
    tmp = 10*log10(dd.Data);
    % [tmp_rline_max_vals, tmp_rline_max_idxs] = max(tmp,[], 1, 'omitnan');
    
    % comparitive
    [I,~,~] = find_multiple(dd.Time, '>=', dd.Surface, 1, 'first');
    surf_next_idxs = [I{:}];
    % distance or closest
    bb = abs(dd.Time - dd.Surface);
    [~, surf_closest_idxs] = min(bb,[],1);
    
    
    figure;
    plot(dd.Surface/1e-6,'.'); hold on; grid on;
    % plot(dd.Time(tmp_rline_max_idxs)/1e-6,'o');
    plot(dd.Time(surf_next_idxs)/1e-6,'rs');
    plot(dd.Time(surf_closest_idxs)/1e-6,'go');
    XMinorGrid ='on'; yticks( dd.Time(min(surf_closest_idxs):max(surf_closest_idxs))/1e-6 );
    legend('Surface','Next', 'Closest');
    xlabel('rlines'); ylabel('us');
    
    
    % transform time axes to WGS-84 elevation
    tt = dd.Time - dd.Time(1);
    [c, er_ice] = physical_constants('c','er_ice');
    big_time = repmat(tt,1,Nx);
    %           big_time = big_time(1:surf_closest_idxs);
    
    for idx_x = 1:Nx
      surf_idx = surf_closest_idxs(idx_x); %%%#########################
      
      air_idxs = 1:surf_idx;
      AGL(idx_x) = dd.Time(surf_idx) *c/2;
      sub_idxs = min(surf_idx+1,Nt):Nt;
      range = [ dd.Time(air_idxs) .*c/2 ; ...
        AGL(idx_x)+(dd.Time(sub_idxs)-dd.Time(surf_idx)) .*(c/2/sqrt(er_ice))];
      elev_axis(:,idx_x) = dd.Elevation(idx_x) - range;
      
    end
    
    hf_9 = figure(999);
    plot(elev_axis,tmp,'.-'); hold on; grid on;
    set(gca, 'Xdir', 'reverse');
    xlabel('elev\_axis, meter');
    ylabel('Magnitude, dB');
    
    figure;
    imagesc(dd.GPS_time, dd.Time/1e-6, tmp, [-100 0] ); hold on;
    plot(dd.GPS_time, dd.Surface/1e-6,'.'); hold on;
    % plot(dd.GPS_time, dd.Time(tmp_rline_max_idxs)/1e-6,'o');
    plot(dd.GPS_time, dd.Time(surf_next_idxs)/1e-6,'rs');
    plot(dd.GPS_time, dd.Time(surf_closest_idxs)/1e-6,'go');
    xlabel('rlines'); ylabel('us');
    legend('Surface','Next', 'Closest');
    XMinorGrid ='on'; yticks( dd.Time(min(surf_closest_idxs):max(surf_closest_idxs))/1e-6 );
    
    %%
    if 0
      figure(3333);
      plot3(dd.Longitude, dd.Latitude, elev_axis); hold on; %, tmp); hold on;
      %           mesh(dd.Longitude, dd.Latitude,  elev_axis, tmp); hold on;
      xlabel('Lon');
      ylabel('Lat');
      zlabel('Elev');
    end
  end
end

return;

unique_elev_axis = unique(elev_axis(:));
size(elev_axis)
prod(size(elev_axis))
figure; plot(unique_elev_axis,'.');
