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
  
elseif 1
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

for idx_season = 1:N_seasons
  for idx_day_seg = 1:N_day_segs{idx_season}
    
    param = params{idx_season}(idx_day_seg);
    fprintf('\n(%s) season (%d of %d) day_seg (%d of %d) >> %s %s ', datestr(now), ...
      idx_season, N_seasons, idx_day_seg, N_day_segs{idx_season}, ...
      param.season_name, param.day_seg);
    
    try
      fprintf('>>records');
      records{idx_season, idx_day_seg} = records_load( param );
      plot(records{idx_season, idx_day_seg}.lat, records{idx_season, idx_day_seg}.lon, '.-'); hold on;
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
    
%     return;
    
    %% meanwhile
    % day_segs = {'20140419_02', '20140419_03'};
    % cursor_info for latlon = 68.820138409998606 -49.543126087195560  68.820137918218123 -49.543124828631676
    % DataIndex = 5172216      564847
    xo_rec_db= { 564847+[0:1], 5172216+[0:1] };
    xo_rec_db= { [5033608]+[0:1], [6546089]+[0:1] };
    xo_rec_db= { [6546089]+[-1:0]; [5033608]+[0:1] };
    
    for idx_season = 1:N_seasons
      for idx_day_seg = 1:N_day_segs{idx_season}
        xo_recs = xo_rec_db{idx_season, idx_day_seg};
        xo_frms = [];
        for idx_rec = 1:length(xo_recs)
          rec = xo_recs(idx_rec);
          frm = find(rec >= frames{idx_season, idx_day_seg}.frame_idxs,1,'last' );
          plot(records{idx_season, idx_day_seg}.lat(rec), records{idx_season, idx_day_seg}.lon(rec),'o');
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
    
    figure;
    for idx_season = 1:N_seasons
      for idx_day_seg = 1:N_day_segs{idx_season}
        fn = fullfile(ct_filename_out(params{idx_season}(idx_day_seg),'standard'), ...
          sprintf('Data_%s_%03d.mat', params{idx_season}(idx_day_seg).day_seg, ...
            xo_frm_db{idx_season, idx_day_seg}));
          dd = load(fn);
          
          dx_vec = distance_geodetic(dd);
          tmp = 10*log10(dd.Data);
          [tmp_rline_max_vals, tmp_rline_max_idxs] = max(tmp,[], 1, 'omitnan');
          
          figure; 
          plot(dd.Surface/1e-6,'.'); hold on; 
          plot(dd.Time(tmp_rline_max_idxs)/1e-6,'o'); 
          xlabel('rlines'); ylabel('us');
          
          figure; 
          imagesc(dd.GPS_time, dd.Time/1e-6, tmp ); hold on; 
          plot(dd.GPS_time, dd.Surface/1e-6,'.'); hold on; plot(dd.GPS_time, dd.Time(tmp_rline_max_idxs)/1e-6,'o');
          xlabel('rlines'); ylabel('us');
          
%           plot3(dd.
      end
    end