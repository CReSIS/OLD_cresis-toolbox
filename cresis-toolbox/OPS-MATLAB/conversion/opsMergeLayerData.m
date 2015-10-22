function lyrOut = opsMergeLayerData(layerDataDayBaseFn,dayStr)
% lyrOut = opsMergeLayerData(layerDataDayBaseFn)
%
% Merges layerData files (for an entire day) into a single concatenated output structure.
%
% Input:
%   layerDataDayBaseFn: Base path to a day's directory of layerData files.
%
% Output:
%   lyrOut: same structure as CReSIS layerdata, with all concatenated data.
%
% Author: Kyle W. Purdon
%
% see also layerDataToOps opsBulkInsert

lyrOut = struct('GPS_time',[],'Latitude',[],'Longitude',[],'Elevation',[],'layerData',...
  {{struct('value',{{struct('data',[]) struct('data',[])}},'quality',[]) ...
  struct('value',{{struct('data',[]) struct('data',[])}},'quality',[])}});

layerDataFns = get_filenames(layerDataDayBaseFn,strcat('Data_',dayStr),'','.mat','recursive');

for fileIdx = 1:length(layerDataFns)
       
  curLayerData = load(layerDataFns{fileIdx});
  if ~isfield(curLayerData, 'layerData')
    % LOAD ECHOGRAM DATA IF LAYERDATA DOES NOT EXIST IN FILE
    curLayerData = uncompress_echogram(curLayerData);
    lyrOut.layerData{1}.value{1}.data = cat(2,lyrOut.layerData{1}.value{1}.data,NaN*zeros(size(curLayerData.Surface)));
    lyrOut.layerData{1}.value{2}.data = cat(2,lyrOut.layerData{1}.value{2}.data,curLayerData.Surface);
    lyrOut.layerData{1}.quality = cat(2,lyrOut.layerData{1}.quality,ones(size(curLayerData.Surface)));
  else
    lyrOut.layerData{1}.value{1}.data = cat(2,lyrOut.layerData{1}.value{1}.data,curLayerData.layerData{1}.value{1}.data);
    lyrOut.layerData{1}.value{2}.data = cat(2,lyrOut.layerData{1}.value{2}.data,curLayerData.layerData{1}.value{2}.data);
    lyrOut.layerData{1}.quality = cat(2,lyrOut.layerData{1}.quality,curLayerData.layerData{1}.quality);
    lyrOut.layerData{2}.value{1}.data = cat(2,lyrOut.layerData{2}.value{1}.data,curLayerData.layerData{2}.value{1}.data);
    lyrOut.layerData{2}.value{2}.data = cat(2,lyrOut.layerData{2}.value{2}.data,curLayerData.layerData{2}.value{2}.data);
    lyrOut.layerData{2}.quality = cat(2,lyrOut.layerData{2}.quality,curLayerData.layerData{2}.quality);
  end
  
  lyrOut.GPS_time = cat(2,lyrOut.GPS_time,curLayerData.GPS_time);
  lyrOut.Latitude = cat(2,lyrOut.Latitude,curLayerData.Latitude);
  lyrOut.Longitude = cat(2,lyrOut.Longitude,curLayerData.Longitude);
  lyrOut.Elevation = cat(2,lyrOut.Elevation,curLayerData.Elevation);
  
end
end
