function layerData_fns = getLayerDataPaths(param_fn,season_dir)
% Get a cell of layerData filenames for an entire seasons based on the
% values in a param sheet. Only processed layerData will be included in the
% output cell array.
%
% Input:
%   param_fn: Absoulute path and filename to a CReSIS param sheet.
%   season_dir: Directory of season data on scratch2.
%
% Output:
%   layerData_fns: A cell array of absoulte paths with filename to all
%   "valid" layerData (.m) files for the input param season.
%
% "Valid" as defined as NOT containing a "do not process" flag in the param
% spreadsheet.
%
% Author: Kyle W. Purdon
%
% Example:
%   param_fn = 'H:\scripts\params\mcords_param_2012_Greenland_P3.xls';
%   season_dir = 'Z:\mdce\mcords2\2012_Greenland_P3\';
%   layerData_fns = getLayerDataPaths(param_fn,season_dir);
%

% Load the param sheet
params = read_param_xls(param_fn);

% Create holder for segment list (filtered w/ dnp)
good_segs = {};

% Get a list of segments and dnp filters
for seg_idx = 1:length(params)
  % Save the segment name
  seg = params(seg_idx).day_seg;
  
  % Check if "do not process" is noted.
  notes_string = params(seg_idx).cmd.notes;
  dnp_flag = ~isempty(strfind(lower(notes_string),'do not process'));
  
  % Save seg if dnp is not (1)
  if ~dnp_flag
    good_segs{end+1,1} = seg;
  end
end

% Get all of the layerData filenames for the season
layerData_fns = get_filenames(fullfile(season_dir,'CSARP_layerData'),'Data_','','.mat','recursive');

% Get a list of segments from the layerData_fns (Segments will repeat)
layerData_fns_segs = cell(length(layerData_fns),1);

for frame_idx = 1:length(layerData_fns)
  layerData_fns_segs{frame_idx} = layerData_fns{frame_idx}(length(layerData_fns{frame_idx})-35:length(layerData_fns{frame_idx})-25);
end

% Get the layerData_fns_segs indeces that match a seg from good_segs. Then
% keep the corresponding layerData_fns as the final product.
keep_paths = {};
for seg_idx = 1:length(good_segs)
  match_idxs = strcmp(good_segs(seg_idx),layerData_fns_segs);
  if any(match_idxs)
    keep_paths{end+1,1} = layerData_fns(match_idxs);
  end
end

layerData_fns = cat(1,keep_paths{:});
end