function run_concatenate_thickness_files(csv_base_dir,kml_base_dir,param)
% run_concatenate_thickness_files(csv_base_dir,kml_base_dir,param)
%
% Concatenates segments into whole season,
% and then creates browse files in CSV and KML format. Note that this
% recreates all of the files in the directory (i.e. it always does the
% whole mission).
%
% Example
%  csv_base_dir = '/cresis/scratch2/mdce/mcords/2009_Antarctica_DC8/CSARP_post/csv';
%  kml_base_dir = '/cresis/scratch2/mdce/mcords/2009_Antarctica_DC8/CSARP_post/kml';
%  param.season_name = '2009_Antarctica_DC8'
%  run_concatenate_thickness_files(csv_base_dir,kml_base_dir,param);
%
% Author: John Paden, Logan Smith, Shashanka Jagarlapudi
%
% See also create_posting.m, concatenate_thickness_files.m

season_fn = fullfile(csv_base_dir, [param.season_name '.csv']);
concatenate_thickness_files(csv_base_dir,'Data_*.csv',season_fn,',');

% Create KML browse file for whole season
kml_out_fn = fullfile(kml_base_dir, ['Browse_' param.season_name '.kml']);
kml_write_cresis(season_fn, kml_out_fn, param.season_name,'mission',[50 0]);

% Create CSV browse file for whole season
[season_path season_name season_ext] = fileparts(season_fn);
browse_season_fn = fullfile(season_path,['Browse_' season_name season_ext]);
fid = fopen(season_fn);
fid_browse = fopen(browse_season_fn,'w');
buffer = fgets(fid);
fwrite(fid_browse,buffer,'char');
num_skip = 0;
while ~feof(fid)
  buffer = fgets(fid);
  if num_skip == 0
    fwrite(fid_browse,buffer,'char');
    num_skip = 50;
  else
    num_skip = num_skip - 1;
  end
end
fclose(fid);
fclose(fid_browse);
  
return;
