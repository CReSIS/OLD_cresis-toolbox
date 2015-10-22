function opsLayerData = layerDataToOps(lyr,settings)
% opsLayerData = layerDataToOps(layerDataFn,settings)
%
% 1. Converts CReSIS layerData to the OPS format.
% 2. Interpolates layerData onto a fixed scale based on point paths in the OPS database.
% 3. Removes any large gaps in the layerData (see also data_gaps_check.m)
% 4. Removes any duplicate points in the layerData (all output are set to type=auto/2)
%
% Input:
%   data: Combined layerData structure from opsMergeLayerData.
%   settings: structure with the following fields
%     .layerFilter (optional) = REGULAR EXPRESSION OF LAYER NAMES TO INSERT
%       for more information on layerFilter see also runOpsBulkInsert.m
%     .location = string ('arctic','antarctic')
%     .seasonName = string
%     .sysName = string ('rds','accum',...)
%
% Output:
%   opsLayerData = structure with fields:
%       properties.point_path_id = integer array 
%       properties.twtt = double array
%       properties.type = integer arry (1 or 2) 1:manual, 2:auto
%       properties.quality = integer array (1, 2 or 3) 1:good, 2:moderate, 3:derived
%       properties.lyr_name = string ('surface','bottom', ...)
%
% Author: Kyle W. Purdon
%
% see also opsGetPath opsCreateLayerPoints

% SET DAFAULTS AND CHECK INPUT
if ~exist('settings','var')
  settings = struct();
end
if ~isfield(settings,'layerFilter') || isempty(settings.layerFilter)
  settings.layerFilter = inline('~isempty(regexp(x,''.*''))');
end

% % LOAD THE LAYERDATA FILE
% lyr = load(layerDataFn,'GPS_time','Latitude','Longitude','Elevation','layerData');
% 
% % LOAD ECHOGRAM DATA IF LAYERDATA DOES NOT EXIST IN FILE
% if ~isfield(lyr, 'layerData')
%   lyr = load(layerData_fn,'GPS_time','Latitude','Longitude','Elevation','Surface','Truncate_Bins','Elevation_Correction','Time');
%   lyr = uncompress_echogram(lyr);
%   lyr.layerData{1}.value{1}.data = NaN*zeros(size(lyr.Surface));
%   lyr.layerData{1}.value{2}.data = lyr.Surface;
%   lyr.layerData{1}.quality = ones(size(lyr.Surface));
% end
% 
% % DONT SUPPORT NEW LAYERDATA FOR NOW (CRESIS > OPS ONLY)
% % CHECK WHICH LAYERDATA FORMAT THIS IS (CReSIS OR OPS)
% % newLd = false;
% if ~isfield(lyr.layerData{1},'value')
%   error('NEW LAYERDATA FORMAT NOT SUPPORTED YET')
%   % newLd = true;
% end

% LOAD AND COMBINE LAYERDATA
% lyr = opsMergeLayerData(layerDataFn,layerDataFn_Pre,layerDataFn_Post);

% GET OPS PATH INFORMATION
opsCmd;
pathParam.properties.location = settings.location;
pathParam.properties.season = settings.seasonName;
pathParam.properties.start_gps_time = min(lyr.GPS_time);
pathParam.properties.stop_gps_time = max(lyr.GPS_time);
pathParam.properties.nativeGeom = true;
[~,pathData] = opsGetPath(settings.sysName,pathParam);

% BUILD UP A STRUCTURE COMBINE AUTO/MANUAL AND REMOVE DUPLICATES
lyrCombined = [];
opsLayerData = [];
for layerIdx = 1:length(lyr.layerData)
  
  % SET DEFAULT LAYER NAMES (FOR SURFACE AND BOTTOM)
  if ~isfield(lyr.layerData{layerIdx},'name')
    if layerIdx == 1
      lyr.layerData{layerIdx}.name = 'surface';
    elseif layerIdx == 2
      lyr.layerData{layerIdx}.name = 'bottom';
    end
  end
  
  % CHECK IF THIS LAYER SHOULD BE PROCESSED
  if ~settings.layerFilter(lyr.layerData{layerIdx}.name)
    continue;
  end
  
  % ADD MANUAL LAYER POINTS (type = 1)
  layerSource = lyr.layerData{layerIdx}.value{1}.data;
  goodIdxs = find(~isnan(layerSource) & isfinite(layerSource));
  lyrCombined.gps_time = lyr.GPS_time(goodIdxs);
  lyrCombined.lat = lyr.Latitude(goodIdxs);
  lyrCombined.lon = lyr.Longitude(goodIdxs);
  lyrCombined.elev = lyr.Elevation(goodIdxs);
  lyrCombined.twtt = double(layerSource(goodIdxs));
  lyrCombined.quality = lyr.layerData{layerIdx}.quality(goodIdxs);
  
  % ADD AUTOMATIC LAYER POINTS (type = 2)
  layerSource = lyr.layerData{layerIdx}.value{2}.data;
  goodIdxs = find(~isnan(layerSource) & isfinite(layerSource));
  lyrCombined.gps_time = cat(2,lyrCombined.gps_time,lyr.GPS_time(goodIdxs));
  lyrCombined.lat = cat(2,lyrCombined.lat,lyr.Latitude(goodIdxs));
  lyrCombined.lon = cat(2,lyrCombined.lon,lyr.Longitude(goodIdxs));
  lyrCombined.elev = cat(2,lyrCombined.elev,lyr.Elevation(goodIdxs));
  lyrCombined.twtt = cat(2,lyrCombined.twtt,double(layerSource(goodIdxs)));
  lyrCombined.quality = cat(2,lyrCombined.quality,lyr.layerData{layerIdx}.quality(goodIdxs));
  
  % SORT COMBINED LAYER DATA BY GPS_TIME
  if ~issorted(lyrCombined.gps_time)
      [lyrCombined.gps_time,sortIdxs] = sort(lyrCombined.gps_time);
      lyrCombined.lat = lyrCombined.lat(sortIdxs);
      lyrCombined.lon = lyrCombined.lon(sortIdxs);
      lyrCombined.elev = lyrCombined.elev(sortIdxs); 
      lyrCombined.twtt = lyrCombined.twtt(sortIdxs);
      lyrCombined.quality = lyrCombined.quality(sortIdxs);
  end
  
  % CHECK FOR EMPTY LAYER COMBINED (WARN AND STOP IF SURFACE)
  if isempty(lyrCombined.twtt)
      if strcmp(lyr.layerData{layerIdx}.name,'surface')
          warning('LayerData file with no surface. No data will be loaded for this frame, please report this.');
          opsLayerData = [];
          return;
      else
          lyrCombined = [];
          continue;
      end
  end
  
  % FIND DUPLICATES AND REMOVE
  [~,notDupIdxs] = unique(lyrCombined.gps_time);
  newGpsTime = nan(size(lyrCombined.gps_time));
  newGpsTime(notDupIdxs) = lyrCombined.gps_time(notDupIdxs);
  lyrCombined.gps_time = newGpsTime;
  clear newGpsTime;
  
  % REMOVE ALL OF THE DUPLICATE DATA FROM THE CURRENT COMBINED LAYERDATA
  keepIdxs = ~isnan(lyrCombined.gps_time);
  lyrCombined.gps_time = lyrCombined.gps_time(keepIdxs);
  lyrCombined.lat = lyrCombined.lat(keepIdxs);
  lyrCombined.lon = lyrCombined.lon(keepIdxs);
  lyrCombined.elev = lyrCombined.elev(keepIdxs);
  lyrCombined.twtt = lyrCombined.twtt(keepIdxs);
  lyrCombined.quality = lyrCombined.quality(keepIdxs);
  lyrCombined.quality(isnan(lyrCombined.quality)) = 1; % CORRECT FOR NAN QUALITY
  lyrCombined.lyr_name = lower(lyr.layerData{layerIdx}.name); % ADD A NAME FIELD
  
  % SUSBSET PATH DATA TO EXTENT OF COMBINED LAYERDATA
  keepPathIdxs = intersect(find(pathData.properties.gps_time <= max(lyrCombined.gps_time)),find(pathData.properties.gps_time >= min(lyrCombined.gps_time)));
  keepPathIdxs = pathData.properties.gps_time <= max(lyrCombined.gps_time) & pathData.properties.gps_time >= min(lyrCombined.gps_time);
  
  % FIND GAPS IN DATA
  pathAlongTrack = geodetic_to_along_track(pathData.properties.Y(keepPathIdxs),pathData.properties.X(keepPathIdxs),pathData.properties.elev(keepPathIdxs));
  master_along_track_at_slave_times = interp1(pathData.properties.gps_time(keepPathIdxs), pathAlongTrack, lyrCombined.gps_time);
  dataGapIdxs = data_gaps_check_mex(pathAlongTrack,master_along_track_at_slave_times,50,20);
  
  % INTERPOLATE COMBINED LAYERDATA ONTO OPS PATH, STORE IN THE OUTPUT
  opsLayerData(end+1).properties.point_path_id = pathData.properties.id(keepPathIdxs);
  opsLayerData(end).properties.twtt = interp1(lyrCombined.gps_time,lyrCombined.twtt,pathData.properties.gps_time(keepPathIdxs),'pchip');
  opsLayerData(end).properties.type = ones(size(pathData.properties.gps_time(keepPathIdxs)))*2;
  opsLayerData(end).properties.quality = interp1(lyrCombined.gps_time,lyrCombined.quality,pathData.properties.gps_time(keepPathIdxs),'nearest');
  opsLayerData(end).properties.lyr_name = lyrCombined.lyr_name;
  
  % REMOVE GAPS IN DATA
  opsLayerData(end).properties.point_path_id = double(opsLayerData(end).properties.point_path_id(~dataGapIdxs));
  opsLayerData(end).properties.twtt = opsLayerData(end).properties.twtt(~dataGapIdxs);
  opsLayerData(end).properties.type = double(opsLayerData(end).properties.type(~dataGapIdxs));
  opsLayerData(end).properties.quality = double(opsLayerData(end).properties.quality(~dataGapIdxs));
  
%   % REMOVE NANs IN DATA AND KEEP ONLY VALUES THAT HAVE A SURFACE
%   noSurfIdxs = setdiff(opsLayerData(end).properties.point_path_id,opsLayerData(1).properties.point_path_id);
%   if ~isempty(noSurfIdxs)
%       warning('%d Bottom Layer Points Without a Surface Were Removed Before Insert',length(noSurfIdxs))
%   end
%   keepIdxs = ~isnan(opsLayerData(end).properties.twtt) & ~ismember(opsLayerData(end).properties.point_path_id,noSurfIdxs);
%   opsLayerData(end).properties.point_path_id = opsLayerData(end).properties.point_path_id(keepIdxs);
%   opsLayerData(end).properties.twtt = opsLayerData(end).properties.twtt(keepIdxs);
%   opsLayerData(end).properties.type = opsLayerData(end).properties.type(keepIdxs);
%   opsLayerData(end).properties.quality = opsLayerData(end).properties.quality(keepIdxs);
  
  lyrCombined = []; % RESET COMBINED LAYERDATA STRUCTURE
  
end
end