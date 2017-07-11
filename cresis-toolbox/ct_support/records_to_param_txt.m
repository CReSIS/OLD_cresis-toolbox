% script records_to_param_csv
%
% Function which extracts the param settings from the records file.
% Useful for recovering parameter settings.
%
% Author: John Paden

%% User Settings
record_dir = '/scratch/csarp_support/records/rds/2016_Antarctica_DC8/';
output_base_dir = '/scratch/tmp/';

%% Read in records files
fns = get_filenames(record_dir, 'records_','','.mat');
params = [];
for fn_idx = 1:length(fns)
  records_fn = fns{fn_idx};
  fprintf('Loading %s\n', records_fn);
  
  load(records_fn,'param_records');
  
  if fn_idx == 1
    params = param_records;
  else
    params(fn_idx).day_seg = param_records.day_seg;
    params(fn_idx).cmd = param_records.cmd;
    params(fn_idx).vectors = param_records.vectors;
    params(fn_idx).records = param_records.records;
    params(fn_idx).get_heights = param_records.get_heights;
    params(fn_idx).csarp = param_records.csarp;
    params(fn_idx).combine = param_records.combine;
    params(fn_idx).radar = param_records.radar;
  end
  
end

%% Save params output file
[output_dir,radar_type,radar_name] = ct_output_dir(params(1).radar_name);
params_fn = fullfile(output_base_dir,sprintf('%s_params_%s.mat',output_dir,params(1).season_name));
fprintf('Output %s\n', params_fn);
save(params_fn, 'params');

%% Save command TXT file
params_cmd_fn = fullfile(output_base_dir,sprintf('%s_params_%s_cmd.txt',output_dir,params(1).season_name));
fprintf('Output %s\n', params_cmd_fn);
fid = fopen(params_cmd_fn,'wb');

fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
  'Date','1','frms','create_records','create_frames','get_heights','csarp','combine_wf_chan','generic','mission_names','notes');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
  'YYYYMMDD','Segment','r','b','b','b','b','b','r','t','t');
for param_idx = 1:length(params)
  param = params(param_idx);
  
  % Print out data row
  fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
    param.day_seg(1:8),param.day_seg(10:11),'','','','','','','',param.cmd.mission_names,param.cmd.notes);
  
end

fclose(fid);

%% Save vectors TXT file
params_vectors_fn = fullfile(output_base_dir,sprintf('%s_params_%s_vectors.txt',output_dir,params(1).season_name));
fprintf('Output %s\n', params_vectors_fn);
fid = fopen(params_vectors_fn,'wb');

fprintf(fid,'\nVectors WorkSheet\n');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'Date', ...
  '1', 'file.start_idx', 'file.stop_idx', 'file.basedir', 'file.adc_folder_name', 'file.prefix', 'file.midfix','file.regexp','gps.time_offset');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
  'YYYYMMDD','Segment','r','r','t','t','t','t','r','r');
for param_idx = 1:length(params)
  param = params(param_idx);
  
  % Print out data row
  fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
    param.day_seg(1:8),param.day_seg(10:11),'','','','','','','',param.cmd.mission_names,param.cmd.notes);
  fprintf(fid,'%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%g\n', param.day_seg(1:8),param.day_seg(10:11), ...
    param.vectors.file.start_idx,param.vectors.file.stop_idx, ...
    param.vectors.file.base_dir,param.vectors.file.adc_folder_name, ...
    param.vectors.file.file_prefix,'', ...
    '',param.vectors.gps.time_offset);
end

fclose(fid);
