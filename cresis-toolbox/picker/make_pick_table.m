function make_pick_table(out_path,param)
%
% Creates a table for use in the picking spreadsheet.
%
% make_pick_table(out_path,param)
%
% out_path: Absoulute path to an output folder.
%
% param: Structure containing the following fields.
%   .radar_name
%   .season_name
%
% EXAMPLE
% -------------------------------------------------
% param.radar_name ='mcords2';
% param.season_name='2012_Greenland_P3';
% out_path = 'C:\Users\kpurdon\Documents\';
% make_pick_table(out_path,param);
% -------------------------------------------------
%
% Author: Kyle W. Purdon

% Check param structure for correct fields.
if sum(isfield(param,{'radar_name','season_name'})) ~= 2
  error('Invalid input structure. param must contain both fields radar_name and season_name. Please try again.');
end

% Create data directory
if ispc
  data_dir = fullfile('Z:\mdce',param.radar_name,param.season_name,'CSARP_layerData\');
else
  data_dir = fullfile('/cresis/scratch2/mdce',param.radar_name,param.season_name,'CSARP_layerData/');
end

% Get list of segments
data_fns = get_filenames(data_dir,'Data','','.mat','recursive');
segs = cell(length(data_fns),1);
for fns_idx = 1:length(data_fns)
  segs{fns_idx} = data_fns{fns_idx}(end-18:end-8);
end
segs = sort(segs)';
[seg_list] = unique(segs)';

% Get the number of frames per segment
frms_per_seg = zeros(length(seg_list),1);
for seg_idx = 1:length(seg_list)
  match_segs = strcmp(seg_list{seg_idx},segs);
  frms_per_seg(seg_idx) = sum(match_segs);
end
  
% Set up variables for output table
output_fn = fullfile(out_path,strcat(param.season_name,'_PickingTable.csv'));
fid = fopen(output_fn,'w+');
% Write header
fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s\n','Complete','Segment','#Frames','Picker','PickStatus','Notes','QC','QCNotes','QCStatus');

% Write data to file
for f_idx = 1:length(seg_list)
  fprintf(fid,'%d,%s,%d,%s,%d,%s,%s,%s,%d\n',0,seg_list{f_idx},frms_per_seg(f_idx),'',0,'','','',0);
end

fclose(fid);

end