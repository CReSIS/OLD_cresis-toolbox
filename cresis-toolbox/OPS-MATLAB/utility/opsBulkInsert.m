function opsBulkInsert(settings)
%
% opsBulkInsert(settings)
%
% Loads data in bulk from the CReSIS filesystem to the OPS.
%
% Author: Kyle W. Purdon
%
% see also runOpsBulkInsert

%% SETUP DEFAULTS AND OTHER INPUTS

global gRadar;
opsCmd;

% SET THE LOG BASE DIRECTORY
if settings.logsOn
  logsBaseFn = gRadar.tmp_path;
  fprintf('Logs will be saved to: %s\n',logsBaseFn);
end

% LAYER FILTER DEFAULTS
if ~isfield(settings,'layerFilter') || isempty(settings.layerFilter)
  settings.layerFilter = inline('~isempty(regexp(x,''.*''))');
end

% PATH SPACING DEFAULT
if ~isfield(settings,'pathSpacing') || isempty(settings.pathSpacing)
  settings.pathSpacing = 15;
end

% GET A BOOLEAN VALUE FOR EACH RUN OPTION
insertPathCmd = any(settings.runType == [1,4,5]);
insertLayerCmd = any(settings.runType == [2,4,5]);
insertAtmCmd = any(settings.runType == [3,5]);

%% LOAD THE PARAM SPREADSHEET AND SAVE LOCAL VARIABLES FROM THE SETTINGS

params = read_param_xls(ct_filename_param(settings.paramFn));
settings.radarName = params(1).radar_name;
settings.seasonName = params(1).season_name;

%% CHECK FILE INPUTS

fprintf('Checking file inputs ...\n');
for param_idx = 1:length(params)
  param = params(param_idx);
  if param.cmd.generic == 1
    if ~isempty(regexpi(param.cmd.notes,'do not process'))
      warning('You have enabled a segment with ''do not process'' in the cmd.notes, dbcont to continue');
      keyboard
    end
    recordsFn = ct_filename_support(param,'','records');
    framesFn = ct_filename_support(param,'','frames');
    layerDir = ct_filename_out(param,settings.layerDataPath,param.day_seg);
    if ~exist(recordsFn,'file') && insertPathCmd
      error('  %s: missing %s\n', param.day_seg, recordsFn);
    elseif ~exist(framesFn,'file') && insertPathCmd
      error('  %s: missing %s\n', param.day_seg, framesFn);
    elseif ~exist(layerDir,'dir') && insertLayerCmd
      error('  %s: missing %s\n', param.day_seg, layerDir);
    else
      fprintf('  %s checked\n', param.day_seg);
    end
  end
end

%% USER CONFIRMATION

switch settings.runType
  case 1
    runStr = 'path';
  case 2
    runStr = 'layer';
  case 3
    runStr = 'atm';
  case 4
    runStr = 'path,layer';
  case 5
    runStr = 'path,layer,atm';
end

confirmParams = {sprintf('SERVER: \t %s',gOps.serverUrl),'',sprintf('SEASON NAME: \t %s',settings.seasonName),'',...
  sprintf('RADAR NAME: \t %s',settings.radarName),'',sprintf('RUN TYPE: \t %s',runStr),'',...
  sprintf('SYSTEM NAME: \t %s',settings.sysName),'',sprintf('PATH SPACING: \t %0.2f meters',settings.pathSpacing),'',...
  sprintf('LOCATION: \t %s',settings.location),''};

confirmButton = questdlg(confirmParams,'Confirm Settings','YES:LOAD','NO:QUIT','NO:QUIT');

switch confirmButton
  case 'NO:QUIT'
    error('PROCESS STOPPED BY USER.');
end

%% PATH INSERTION

if insertPathCmd
  
  failedSegments = {}; % STORE SEGMENTS THAT FAIL
  
  % START LOGGING
  if settings.logsOn
    pathLogFn = fullfile(logsBaseFn,strcat('OPS_PATHINSERT_',datestr(now,'yyyy.mm.dd'),'_',datestr(now,'HH.MM.SS'),'.txt'));
    diary(pathLogFn);
  end
  
  % FOR EACH SEGMENT PROCESS THE INPUT AND PUSH THE DATA TO THE SERVER
  for paramIdx = 1:length(params)
    
    try
      
      % CONFIRM THAT GENERIC IS NOT FLAGGED
      param = params(paramIdx);
      if param.cmd.generic ~= 1
        continue;
      end
      
      fprintf('Loading path for segment %s ... ',param.day_seg);
      
      start = tic; % SET UP TIMING
      
      % GET THE RECORDS AND FRAMES FILENAMES
      recordsFn = ct_filename_support(param,'','records');
      framesFn = ct_filename_support(param,'','frames');
      
      load(framesFn); % LOAD FRAMES FILE
      records = load(recordsFn); % LOAD RECORDS FILE
      if isfield(records,'records')
        records = records.records; % SUPPORT OLD RECORDS FORMAT
      end
      
      % GET THE START GPS TIME OF EACH FRAME
      frameStartGpsTime = zeros(1,length(frames.frame_idxs));
      for frmIdx = 1:length(frames.frame_idxs)
        if frmIdx == length(frames.frame_idxs)
          frameStartGpsTime(1,frmIdx) = records.gps_time(frames.frame_idxs(frmIdx)); % ADD THE END SEGMENT
        else
          frameStartGpsTime(1,frmIdx) = records.gps_time(frames.frame_idxs(frmIdx)); % ADD 1:N-1 SEGMENTS
        end
      end
      
      % INTERPOLATE RECORDS GPS TIME ONTO THE GIVEN SPACING (DEFAULT = 15m)
      alongTrack = geodetic_to_along_track(records.lat,records.lon,records.elev);
      newAlongTrack = 0:settings.pathSpacing:alongTrack(end);
      positive_idxs = [1, 1+find(diff(alongTrack) > 0)];
      outGpsTime = interp1(alongTrack(positive_idxs),records.gps_time(positive_idxs),newAlongTrack,'pchip');
      
      % INTERPOLATE RECORDS VALUES ONTO NEW GPS TIME
      outLon = interp1(records.gps_time,records.lon,outGpsTime,'pchip');
      outLat = interp1(records.gps_time,records.lat,outGpsTime,'pchip');
      outElev = interp1(records.gps_time,records.elev,outGpsTime,'pchip');
      outRoll = interp1(records.gps_time,records.roll,outGpsTime,'pchip');
      outPitch = interp1(records.gps_time,records.pitch,outGpsTime,'pchip');
      outHeading = interp1(records.gps_time,records.heading,outGpsTime,'pchip');
      
      % ERROR CHECK OUTPUT DATA
      if any(find(outHeading>(pi*2))) || any(find(outHeading<(-pi*2)))
        warning('OUTPUT HEADING OUT OF BOUND 2pi <> -2pi');
        keyboard;
      end
      if any(find(outPitch>(pi*2))) || any(find(outPitch<(-pi*2)))
        warning('OUTPUT PITCH OUT OF BOUND 2pi <> -2pi');
        keyboard;
      end
      if any(find(outRoll>(pi*2))) || any(find(outRoll<(-pi*2)))
        warning('OUTPUT ROLL OUT OF BOUND 2pi <> -2pi');
        keyboard;
      end
      if any(find(outLon>180)) || any(find(outLon<-180))
        warning('OUTPUT LONGITUDE OUT OF BOUND 180 <> -180');
        keyboard;
      end
      if any(find(outLat>90)) || any(find(outLat<-90))
        warning('OUTPUT LATITUDE OUT OF BOUND 90 <> -90');
        keyboard;
      end
      
      % BUILD STRUCTURE FOR opsCreatePath()
      outData.geometry.coordinates = [outLon' outLat'];
      outData.properties.location = settings.location;
      outData.properties.season = settings.seasonName;
      outData.properties.season_group = settings.seasonGroup;
      outData.properties.radar = settings.radarName;
      outData.properties.segment = param.day_seg;
      outData.properties.gps_time = outGpsTime;
      outData.properties.elev = outElev;
      outData.properties.roll = outRoll;
      outData.properties.pitch = outPitch;
      outData.properties.heading = outHeading;
      outData.properties.frame_count = length(frames.frame_idxs);
      outData.properties.frame_start_gps_time = frameStartGpsTime;
      
      clear records frames
      
      mstop = toc(start); % RECORD MATLAB COMPUTATION TIME
      
      % PUSH DATA TO THE SERVER
      [status,message] = opsCreatePath(settings.sysName,outData);
      
      pstop = toc(start)-mstop; % RECORD PYTHON COMPUTATION TIME
      
      if status ~= 1
        fprintf('\n');
        warning(message);
        failedSegments{end+1} = param.day_seg; % STORE THE SEGMENT IF THE SERVER PUSH FAILED
      end
      
      % REPORT TIMING AND STATUS
      fprintf('\n\t-> Time: Matlab %2.2fs Python %2.2fs\n',mstop,pstop);
      fprintf('\t-> Status: %s\n',message);
      
    catch ME
      
      diary OFF
      ME.getReport()
      failedSegments{end+1} = param.day_seg; % STORE THE SEGMENT IF ANYTHING FAILED
      
    end
  end
  
  % REPORT FAILED SEGMENTS
  if ~isempty(failedSegments)
    fprintf('\n\nThere were issues loading paths for segments:\n');
    for failIdx = 1:length(failedSegments)
      fprintf('\t%s\n',failedSegments{failIdx});
    end
    fprintf('No layer points will be loaded for these segments.\n\n');
  end
  
  diary OFF % STOP LOGGING
  
end

%% LAYER INSERTION

if insertLayerCmd
  
  mstart = tic;
  
  failedFrames = {}; % STORE FRAMES THAT FAIL
  failedLayers = {}; % STORE LAYERS THAT FAIL
  frameName = '';
  segDayObj = [];
  opsLayerDataSub = [];
  
  if ~exist('failedSegments','var')
    failedSegments = {};
  end
  
  % START LOGGING
  if settings.logsOn
    pathLogFn = fullfile(logsBaseFn,strcat('OPS_LAYERINSERT_',datestr(now,'yyyy.mm.dd'),'_',datestr(now,'HH.MM.SS'),'.txt'));
    diary(pathLogFn);
  end
  
  % BUILD THE DAY / DAYSEG LIST FOR PROCESSING
  for paramIdx = 1:length(params)
    
    % CONFIRM THAT GENERIC IS NOT FLAGGED
    param = params(paramIdx);
    if param.cmd.generic ~= 1
      continue;
    end
    
    % DO NOT LOAD LAYERS FOR FAILTED SEGMENTS
    if any(strcmp(param.day_seg,failedSegments))
      fprintf('Skipping layers for segment %s, path loading failed for this segment.\n',param.day_seg);
      continue;
    end
    
    load(ct_filename_support(param,'','frames')); % LOAD FRAMES FILE
    
    % CROSS VERIFY FRAMES FILE WITH PARAM.CMD.FRMS LIST
    if isempty(param.cmd.frms)
      param.cmd.frms = 1:length(frames.frame_idxs);
    end
    [validFrms,keepIdxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
    if length(validFrms) ~= length(param.cmd.frms)
      badMask = ones(size(param.cmd.frms));
      badMask(keepIdxs) = 0;
      warning('Nonexistent frames in param.cmd.frms (e.g. frame "%g" is invalid). These will be removed.',param.cmd.frms(find(badMask,1)));
      param.cmd.frms = validFrms;
    end
    
    for frmId = 1:length(param.cmd.frms)
      
      % BUILD UP DAY/DAYSEG/DAYSEGFRM OBJECT
      segDayObj = cat(1,segDayObj,{param.day_seg(1:end-3) param.day_seg sprintf('%s_%03d',param.day_seg,frmId)});
      
    end
    
  end
  
  % GET UNIQUE LIST OF DAYS
  segDays = unique(segDayObj(:,1));
  
  % PROCESS EACH DAY
  for segDayIdx = 1:length(segDays)
    
    curDay = segDays{segDayIdx};
    fprintf('Processing data for day %s ...\n',curDay);
    
    % LOAD ALL THE LAYERDATA FOR THE CURRENT DAY
    layerBaseDir = ct_filename_out(param,settings.layerDataPath,'',true);
    dayLayerData = opsMergeLayerData(layerBaseDir,curDay);
    opsLayerData = layerDataToOps(dayLayerData,settings);
    
    % CHECK FOR EMPTY LAYERS
    emptyLayerIdxs = [];
    for layerIdx = 1:length(opsLayerData)
      if isempty(opsLayerData(layerIdx).properties.twtt)
        emptyLayerIdxs(end+1) = layerIdx;
      end
    end
    opsLayerData(emptyLayerIdxs) = []; % SET THE EMPTY STRUCTURES
    if isempty(opsLayerData)
      message = 'WARNING: NO LAYERS TO INSERT, LAYERS ARE EMPTY'; % SAVE MESSAGE IF ALL ARE EMPTY
      fprintf('\n\t\t\t-> Time: Matlab %2.2fs Python %2.2fs\n',0.00,0.00);
      fprintf('\t\t\t-> Status: %s\n',message);
      failedFrames{end+1} = frameName;
    end
    
    % GET A LIST OF SEGMENTS FOR THE CURRENT DAY
    segList = unique(segDayObj(find(strcmp(segDayObj(:,1),curDay)),2));
    
    % PROCESS EACH SEGMENT
    for segIdx = 1:length(segList)
      
      curSeg = segList{segIdx};
      fprintf('Loading layers for segment %s ...\n',curSeg);
      
      segInfoParam.properties.segment = curSeg;
      segInfoParam.properties.season = settings.seasonName;
      [~,segData] = opsGetSegmentInfo(settings.sysName,segInfoParam);
      
      mstop = toc(mstart); % RECORD MATLAB COMPUTATION TIME
      
      % PROCESS EACH FRAME
      for frmIdx = 1:length(segData.properties.frame)
        
        try
          curFrame = segData.properties.frame{frmIdx};
          if ~any(strcmp(segDayObj(:,3),curFrame))
            continue;
          end
          
          fprintf('\tLoading frame %s ...\n',curFrame);
          
          % PROCESS EACH LAYER
          for layerIdx = 1:length(opsLayerData)
            
            curLayer = opsLayerData(layerIdx).properties.lyr_name;
            
            try
              
              ptic = tic;
              fprintf('\t\tLoading layer %s ...',curLayer);
              
              % GET POINT PATH IDS FOR FRAME TIME RANGE
              pathParam.properties.location = settings.location;
              pathParam.properties.season = settings.seasonName;
              pathParam.properties.start_gps_time = segData.properties.start_gps_time(frmIdx);
              pathParam.properties.stop_gps_time = segData.properties.stop_gps_time(frmIdx);
              [~,pathData] = opsGetPath(settings.sysName,pathParam);
              [~,keepIdxs] = ismember(pathData.properties.id,opsLayerData(layerIdx).properties.point_path_id);
              keepIdxs = keepIdxs(find(keepIdxs~=0));
              
              % SUBSET OPSLAYERDATA TO FRAME
              opsLayerDataSub.properties.point_path_id = opsLayerData(layerIdx).properties.point_path_id(keepIdxs);
              opsLayerDataSub.properties.twtt = opsLayerData(layerIdx).properties.twtt(keepIdxs); 
              opsLayerDataSub.properties.type = opsLayerData(layerIdx).properties.type(keepIdxs); 
              opsLayerDataSub.properties.quality = opsLayerData(layerIdx).properties.quality(keepIdxs); 
              opsLayerDataSub.properties.lyr_name = opsLayerData(layerIdx).properties.lyr_name; 

              % PUSH DATA TO THE SERVER
              [status,message] = opsCreateLayerPoints(settings.sysName,opsLayerDataSub);
              pstop = toc(ptic);
              
              opsLayerDataSub = []; % RESET SUBSET STRUCTURE
              
              if status ~= 1
                
                % REPORT A FAILED LAYER
                failedLayers{end+1} = {curFrame,opsLayerData(layerIdx).properties.lyr_name};
                fprintf('\n');
                warning(message);
                
              else
                
                fprintf('\n\t\t\t-> Time: Matlab %2.2fs Python %2.2fs\n',mstop,pstop);
                fprintf('\t\t\t-> Status: %s\n',message);
                
              end
              
            catch ME
              
              % REPORT A FAILED FRAME
              fprintf('\n');
              warning(sprintf('%s at line %d in file %s.',ME.message,ME.stack(1).line,ME.stack(1).name));
              failedLayers{end+1} = {curFrame,opsLayerData(layerIdx).properties.lyr_name}; % STORE THE LAYER NAME IF ANYTHING FAILED
              continue
              
            end
          end
          
        catch ME
          
          fprintf('\n');
          warning(sprintf('%s at line %d in file %s.',ME.message,ME.stack(1).line,ME.stack(1).name));
          failedFrames{end+1} = curFrame; % STORE THE FRAME NAME IF ANYTHING FAILED
          continue
          
        end
      end
      mstart = tic; % RESTART MATLAB TIMING (END OF SEGMENT LOAD)
    end
  end
  
  % REPORT FAILED FRAMES
  if ~isempty(failedFrames)
    fprintf('\n\nThere were issues loading data for frames:\n');
    for failIdx = 1:length(failedFrames)
      fprintf('\t%s\n',failedFrames{failIdx});
    end
  end
  
  % REPORT FAILED LAYERS
  if ~isempty(failedLayers)
    fprintf('\n\nThere were issues loading data for layers:\n');
    for failIdx = 1:length(failedLayers)
      fprintf('\t%s %s\n',failedLayers{failIdx}{1},failedLayers{failIdx}{2});
    end
  end
  
  diary OFF
end

%% ATM INSERTION

if insertAtmCmd
  
  atmFailed=false;
  segDayObj = [];
  
  if ~exist('failedSegments','var')
    failedSegments = {};
  end
  
  % START LOGGING
  if settings.logsOn
    pathLogFn = fullfile(logsBaseFn,strcat('OPS_ATMINSERT_',datestr(now,'yyyy.mm.dd'),'_',datestr(now,'HH.MM.SS'),'.txt'));
    diary(pathLogFn);
  end
  
  start = tic; % SET UP TIMING
  
  % BUILD THE DAY / DAYSEG LIST FOR PROCESSING
  for paramIdx = 1:length(params)
    
    % CONFIRM THAT GENERIC IS NOT FLAGGED
    param = params(paramIdx);
    if param.cmd.generic ~= 1
      continue;
    end
    
    % DO NOT LOAD LAYERS FOR FAILTED SEGMENTS
    if any(strcmp(param.day_seg,failedSegments))
      fprintf('Skipping layers for segment %s, path loading failed for this segment.\n',param.day_seg);
      continue;
    end
    
    % BUILD UP DAY/DAYSEG/DAYSEGFRM OBJECT
    segDayObj = cat(1,segDayObj,{param.day_seg(1:end-3) param.day_seg});
    
  end
  
  % CREATE PARAM FOR ATM LAYER
  opsAtmLayerParam.properties.lyr_name = 'atm';
  opsAtmLayerParam.properties.lyr_group_name = 'lidar';
  opsAtmLayerParam.properties.lyr_description = 'atm l2 lidar surface';
  opsAtmLayerParam.properties.public = true;
  
  % PUSH DATA TO THE SERVER
  try
    [~,~] = opsCreateLayer(settings.sysName,opsAtmLayerParam);
  catch ME
    fprintf('\t-> Layer ''atm'' exists, layer points will be added to layer.\n');
  end
  
  % GET UNIQUE LIST OF DAYS
  segDays = unique(segDayObj(:,1));
  
  % PROCESS EACH DAY
  for segDayIdx = 1:length(segDays)
    
    try
      curDay = segDays{segDayIdx};
      fprintf('Loading ATM data for day %s ...\n',curDay);
      
      % GET ATM L2 FILENAMES
      atmFns = get_filenames_atm(settings.location,curDay,settings.data_support_path);
      if isempty(atmFns)
        warning('No atm data. Skipping segment %s',param.day_seg);
        continue;
      end
      
      % CONVERT ATM L2 TO OPS FORMAT
      opsAtmData = atmToOps(atmFns,settings);
      
      mstop = toc(start); % RECORD MATLAB COMPUTATION TIME
      ptic = tic;
      
      % PUSH DATA TO THE SERVER
      [status,message] = opsCreateLayerPoints(settings.sysName,opsAtmData);
      
      pstop = toc(ptic);
      
      if status ~= 1
        
        atmFailed = true;
        fprintf('\n');
        warning(message);
        
      else
        
        fprintf('\t-> Time: Matlab %2.2fs Python %2.2fs\n',mstop,pstop);
        fprintf('\t-> Status: %s\n',message);

      end
      
    catch ME
      fprintf('\n');
      warning(sprintf('%s at line %d in file %s.',ME.message,ME.stack(1).line,ME.stack(1).name));
      atmFailed = true;
      continue
    end
    
  end
  
  % REPORT FAILED ATM
  if atmFailed
    fprintf('\n\nThere were issues loading atm data.\n');
  end
  
  diary OFF
  
end