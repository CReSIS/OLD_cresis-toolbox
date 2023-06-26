function dbf_to_txt(input_dir_fn,output_dir)
% Converts Shapefile attributes (DBF) to TXT
%
% Format: dbf_to_txt(input_dir_fn,output_dir)
%
% input_dir_fn = absolute path with filename and extension for the DBF.
% output_dir = absoulte path with filename and extention (.txt only) for the output.
%
% INPUT DBF MUST BE FROM "Interpolate CReSIS Data" (ArcGIS Tool)
% REQUIRED FORMAT:
% 'A_SURF','A_BED','THICK','DATATYPE','SEASON','LAT','LON','UTCTime','ELEVATION','YYYYMMDD','SEGMENT','FRAME','SURFACE','BOTTOM','QUALITY'\n
%
% See also dbfread.m

% Create an Output Dir & Filename
sep_idxs = strfind(input_dir_fn,filesep);
filename = input_dir_fn(sep_idxs(end)+1:length(input_dir_fn)-4);
output_dir_fn = fullfile(output_dir);

% Read the attribute data from the DBF
[dbf_data,dbf_header] = dbfread(input_dir_fn);

% Check for the correct number of attributes.
if length(dbf_header) ~= 15
  error('Incorrect fields in DBF file, Check SHP attributes');
end

% CReSIS Shapefile DBF FORMAT
% 'A_SURF','A_BED','THICK','DATATYPE','SEASON','LAT','LON','UTCTime','ELEVATION','YYYYMMDD','SEGMENT','FRAME','SURFACE','BOTTOM','QUALITY'

output_fid = fopen(output_dir_fn,'w+');
% Print Header to file
fprintf(output_fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
  'A_SURF','A_BED','THICK','DATATYPE','SEASON','LAT','LON','UTCTime','ELEVATION','YYYYMMDD','SEGMENT','FRAME','SURFACE','BOTTOM','QUALITY');
% Print Data to file  
% Pay attention to the sequence of fields
for print_idx = 1:length(dbf_data)
  fprintf(output_fid,'%6.2f,%6.2f,%6.2f,%s,%s,%2.6f,%2.6f,%5.4f,%4.4f,%d,%d,%d,%6.2f,%6.2f,%01d\n',...
    dbf_data{print_idx,13},dbf_data{print_idx,14},dbf_data{print_idx,4},dbf_data{print_idx,15},dbf_data{print_idx,12},...
    dbf_data{print_idx,1},dbf_data{print_idx,2},dbf_data{print_idx,3},dbf_data{print_idx,5},dbf_data{print_idx,6},...
    dbf_data{print_idx,7},dbf_data{print_idx,8},dbf_data{print_idx,9},dbf_data{print_idx,10},dbf_data{print_idx,11});
end
close_fid = fclose(output_fid);
end