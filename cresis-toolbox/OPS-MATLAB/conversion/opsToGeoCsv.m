% CONVERT OPS2 CSV FORMAT TO GEOGRAPHIC SEARCH GUI FORMAT

% OPS2 FORMAT
% LAT,LON,ELEVATION,ROLL,PITCH,HEADING,UTCSOD,UTCDATE,SURFACE,BOTTOM,THICKNESS,SURFACE_TYPE,BOTTOM_TYPE,SURFACE_QUALITY,BOTTOM_QUALITY,SEASON,FRAME
% -42.43523945,68.75737351,3303.50817,0.05867,0.03308,1.59267,66668.000,20080627,726.997,9312.405,8585.409,2,2,1,1,2008_Greenland_TO,20080627_06_005
%f%f%f%f%f%f%f%s%f%f%f%d%d%d%d%s%s

% GEOGRAPHIC SEARCH GUI FORMAT
% LAT,LON,UTCTime,THICK,ELEVATION,FRAME,SURFACE,BOTTOM,QUALITY,SEASON
% 80.052114,-58.091100,55501.1396,1100.72,1478.7405,2010032401040,442.35,1543.07,1,2010_Greenland_DC8,
%f%f%f%f%f%s%f%f%d%s

%% =========================================================================
% USER INPUT

opsCsvFn = 'C:\Users\kpurdon.HOME\Downloads\OPS_CReSIS_L2_CSV_GOOD_jbK1TzFfax.csv';
outFn = 'C:\Users\kpurdon.HOME\Downloads\OPS_CReSIS_L2_CSV_GOOD_jbK1TzFfax_CROSSOVER2.csv';
keepEvery = 5;

%% =========================================================================
% AUTOMATED SECTION

% OPEN OPS2 CSV FILE
fid = fopen(opsCsvFn,'r');
if fid == -1
  error('Failed to open csv file.');
end

% READ OPS2 CSV FILE
opsCsvData = textscan(fid,'%f%f%f%f%f%f%f%s%f%f%f%d%d%d%d%s%s','delimiter',',','headerlines',1);

% CONVERT UTCDATE,UTCSOD TO GPSTIME
outTime = zeros(length(opsCsvData{7}),1);
for timeIdx = 1:length(opsCsvData{7})
  outTime(timeIdx) = datenum(str2double(opsCsvData{8}{timeIdx}(1:4)),str2double(opsCsvData{8}{timeIdx}(5:6)),str2double(opsCsvData{8}{timeIdx}(7:8)),0,0,opsCsvData{7}(timeIdx));
end

% ORDER DATA
[~,orderIdxs] = sort(outTime);
for dataIdx = 1:length(opsCsvData)
  opsCsvData{dataIdx} = opsCsvData{dataIdx}(orderIdxs);
end

% CLOSE OPS2 CSV FILE
cid = fclose(fid);
if cid == -1
  fprintf('Failed to close csv file.\n')
end

% OPEN OUTPUT CSV FILE
outFid = fopen(outFn,'w+');
if outFid == -1
  error('Failed to open new csv file.');
end

% WRITE OUTPUT CSV FILE
for dataIdx = 1:keepEvery:length(opsCsvData{1})
  
  fprintf(outFid,'%f,%f,%f,%f,%f,%s,%f,%f,%d,%s\n',...
    opsCsvData{1}(dataIdx),...
    opsCsvData{2}(dataIdx),...
    opsCsvData{7}(dataIdx),...
    opsCsvData{11}(dataIdx),...
    opsCsvData{3}(dataIdx),...
    opsCsvData{17}{dataIdx},...
    opsCsvData{9}(dataIdx),...
    opsCsvData{10}(dataIdx),...
    opsCsvData{15}(dataIdx),...
    opsCsvData{16}{dataIdx});
  
end

% CLOSE OUTPUT CSV FILE
outCid = fclose(outFid);
if outCid == -1
  fprintf('Failed to close new csv file.\n')
end