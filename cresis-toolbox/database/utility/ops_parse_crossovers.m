% ==================================================================
% USER INPUT SECTION
% ==================================================================

% PATH TO CROSSOVERS MAT FILE (FROM run_ops_get_crossovers)
load('Y:\kpurdon\EXTERNAL\paden\crossovers\output\rdsantarcticCrossovers.mat');

% SCATTER PLOT
plotScatter = false;

% MIMIMUM CROSSOVER ERROR TO REPORT (METERS)
reportBasePath = 'C:\Users\tstafford\Documents\PARSE_CROSSES_REPORT\';
minReportError = 0;

% ==================================================================
% AUTOMATED SECTION BELOW
% ==================================================================

% LOAD CONSTANTS
physical_constants;

% CONSTRUCT OUTPUT VECTORS
data.out.surfaceError =[];
data.out.bottomError = [];
data.out.surf.X = [];
data.out.surf.Y = [];
data.out.surf.frame_1 = {};
data.out.surf.frame_2 = {};
data.out.surf.gps_time_1 = [];
data.out.surf.gps_time_2 = [];

data.out.bot.X = [];
data.out.bot.Y = [];
data.out.bot.frame_1 = {};
data.out.bot.frame_2 = {};
data.out.bot.gps_time_1 = [];
data.out.bot.gps_time_2 = [];

% SPLIT CROSSOVERS BY LAYER
% SURFACE
surf_idxs = data.properties.lyr_id == 1;
surf_X = data.properties.X(surf_idxs);
surf_Y = data.properties.Y(surf_idxs);
surf_twtt1 = data.properties.twtt_1(surf_idxs);
surf_twtt2 = data.properties.twtt_2(surf_idxs);
surf_error = data.properties.abs_error(surf_idxs);
surf_gps1 = data.properties.gps_time_1(surf_idxs);
surf_gps2 = data.properties.gps_time_2(surf_idxs);
surf_elev1 = data.properties.ELEV_1(surf_idxs);
surf_elev2 = data.properties.ELEV_2(surf_idxs);
surf_frame1 = data.properties.frame_1(surf_idxs);
surf_frame2 = data.properties.frame_2(surf_idxs);

%BOTTOM
bot_idxs = data.properties.lyr_id == 2;
bot_X = data.properties.X(bot_idxs);
bot_Y = data.properties.Y(bot_idxs);
bot_twtt1 = data.properties.twtt_1(bot_idxs);
bot_twtt2 = data.properties.twtt_2(bot_idxs);
bot_error = data.properties.abs_error(bot_idxs);
bot_gps1 = data.properties.gps_time_1(bot_idxs);
bot_gps2 = data.properties.gps_time_2(bot_idxs);
bot_elev1 = data.properties.ELEV_1(bot_idxs);
bot_elev2 = data.properties.ELEV_2(bot_idxs);
bot_frame1 = data.properties.frame_1(bot_idxs);
bot_frame2 = data.properties.frame_2(bot_idxs);

% WRITE SURFACE CROSSOVERS
for scross_idx = 1:length(surf_X)
    data.out.surfaceError(end+1) = (surf_elev1(scross_idx) - surf_twtt1(scross_idx)*c/2) - (surf_elev2(scross_idx) - surf_twtt2(scross_idx)*c/2);
    data.out.surf.X(end+1) = surf_X(scross_idx);
    data.out.surf.Y(end+1) = surf_Y(scross_idx);
    data.out.surf.frame_1{end+1} = surf_frame1{scross_idx};
    data.out.surf.frame_2{end+1} = surf_frame2{scross_idx};
    data.out.surf.gps_time_1(end+1) = surf_gps1(scross_idx);
    data.out.surf.gps_time_2(end+1) = surf_gps2(scross_idx);
end

% WRITE ALL BOTTOM CROSSOVERS THAT HAVE MATCHING SURFACE POINTS.
for bcross_idx = 1:length(bot_X)
   % IF BOTH BOTTOM POINTS HAVE A CORRESPONDING SURFACE, WRITE
    if sum(surf_gps1 == bot_gps1(bcross_idx)) > 0 && sum(surf_gps2 == bot_gps2(bcross_idx)) > 0
       % Get the matching surface crossover points 
       match_surf_idx1 = surf_gps1 == bot_gps1(bcross_idx);
       match_surf_idx2 = surf_gps2 == bot_gps2(bcross_idx);
       % MAKE SURE THERE ARE NOT MULTIPLE SURF POINTS WITH SAME GPS & DIFF
       % ELEVS
       if sum(match_surf_idx1) > 1 || sum(match_surf_idx2) > 1
           if (sum(match_surf_idx1) > 1 && length(unique(surf_elev1(match_surf_idx1))) > 1) || (sum(match_surf_idx2) > 1 && length(unique(surf_elev2(match_surf_idx2))) > 1)
               error('MULTIPLE SURFACES WITH SAME GPSTIME AND DIFFERENT ELEVS')
           end      
       end
       % EXTRACT NEEDED MATCHING VARIABLES...
       match_surf_elev1 = surf_elev1(match_surf_idx1);
       match_surf_elev2 = surf_elev2(match_surf_idx2);
       match_surf_twtt1 = surf_twtt1(match_surf_idx1);
       match_surf_twtt2 = surf_twtt2(match_surf_idx2);
       % WRITE BOTTOM ERROR
       data.out.bottomError(end+1) = (match_surf_elev1(1) - match_surf_twtt1(1)*c/2 - (bot_twtt1(bcross_idx) - match_surf_twtt1(1))*c/2/sqrt(er_ice)) - (match_surf_elev2(1) - match_surf_twtt2(1)*c/2 - (bot_twtt2(bcross_idx) - match_surf_twtt2(1))*c/2/sqrt(er_ice));
       data.out.bot.X(end+1) = bot_X(bcross_idx);
       data.out.bot.Y(end+1) = bot_Y(bcross_idx);
       data.out.bot.frame_1{end+1} = bot_frame1{bcross_idx};
       data.out.bot.frame_2{end+1} = bot_frame2{bcross_idx};
       data.out.bot.gps_time_1(end+1) = bot_gps1(bcross_idx);
       data.out.bot.gps_time_2(end+1) = bot_gps2(bcross_idx);
    end
   % DON'T WRITE A BOTTOM ERROR IF NO MATCHING SURFACES
end

% CREATE THE SCATTER PLOT IF ASKED
if plotScatter
  fprintf('   Plottin ...\n');
  figure;
  scatter(data.out.X,data.out.Y,50,abs(data.out.surfaceError),'filled');
  title('Absolute Surface Crossover Error');
  colorbar;
  figure;
  scatter(data.out.X,data.out.Y,50,abs(data.out.bottomError),'filled');
  title('Absolute Bottom Crossover Error');
  colorbar;
end

fprintf('   Writing SURFACE output log ...\n');
surfFid = fopen(fullfile(reportBasePath,'surfaceErrorLog.txt'),'w+');
fprintf(surfFid,'SEASON1,FRAME1,SEASON2,FRAME2,ABSERROR,GPS1,GPS2\n');

% CREATE SURFACE OVER LIMIT REPORT
for outIdx = 1:length(data.out.surf.gps_time_1)
  if abs(data.out.surfaceError(outIdx)) > minReportError
    % QUERY FOR SEASON
    query1 = sprintf('SELECT season_name from rds_seasons WHERE season_id=(SELECT season_id from rds_segments where segment_id=(SELECT segment_id from rds_frames WHERE frame_name=''%s''));',data.out.surf.frame_1{outIdx});
    [s,d1] = ops_query(query1);
    if s == 2
      d1{1} = 'UNKNOWN';
    end
    query2 = sprintf('SELECT season_name from rds_seasons WHERE season_id=(SELECT season_id from rds_segments where segment_id=(SELECT segment_id from rds_frames WHERE frame_name=''%s''));',data.out.surf.frame_2{outIdx});
    [s,d2] = ops_query(query2);
    if s == 2
      d2{1} = 'UNKNOWN';
    end
    
    % WRITE ERRORS TO FILE IF THEY ARE > MINREPORTERROR
    if abs(data.out.surfaceError(outIdx)) > minReportError
      fprintf(surfFid,'%s,%s,%s,%s,%f,%f,%f\n',d1{1},data.out.surf.frame_1{outIdx},d2{1},data.out.surf.frame_2{outIdx},abs(data.out.surfaceError(outIdx)),data.out.surf.gps_time_1(outIdx),data.out.surf.gps_time_2(outIdx));
    end
    
  end
end
fclose(surfFid);

fprintf('   Writing BOTTOM output log ...\n');
bedFid = fopen(fullfile(reportBasePath,'bottomErrorLog.txt'),'w+');
fprintf(bedFid,'SEASON1,FRAME1,SEASON2,FRAME2,ABSERROR,GPS1,GPS2\n');
% CREATE BOTTOM OVER LIMIT REPORT
for outIdx = 1:length(data.out.bot.gps_time_1)
  if abs(data.out.bottomError(outIdx)) > minReportError
    % QUERY FOR SEASON
    query1 = sprintf('SELECT season_name from rds_seasons WHERE season_id=(SELECT season_id from rds_segments where segment_id=(SELECT segment_id from rds_frames WHERE frame_name=''%s''));',data.out.bot.frame_1{outIdx});
    [s,d1] = ops_query(query1);
    if s == 2
      d1{1} = 'UNKNOWN';
    end
    query2 = sprintf('SELECT season_name from rds_seasons WHERE season_id=(SELECT season_id from rds_segments where segment_id=(SELECT segment_id from rds_frames WHERE frame_name=''%s''));',data.out.bot.frame_2{outIdx});
    [s,d2] = ops_query(query2);
    if s == 2
      d2{1} = 'UNKNOWN';
    end
    
    % WRITE ERRORS TO FILE IF THEY ARE > MINREPORTERROR
    if abs(data.out.bottomError(outIdx)) > minReportError
      fprintf(bedFid,'%s,%s,%s,%s,%f,%f,%f\n',d1{1},data.out.bot.frame_1{outIdx},d2{1},data.out.bot.frame_2{outIdx},abs(data.out.bottomError(outIdx)),data.out.bot.gps_time_1(outIdx),data.out.bot.gps_time_2(outIdx));
    end 
  end
end
fclose(bedFid);