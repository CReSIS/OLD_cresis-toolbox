% Function for copying backups in the field from the scratch drive to the
% NAS (or some other archive).

param = [];
param.season_name = '2016_Antarctica_DC8';
input_dirs = {};
input_dirs{end+1} = {gRadar.support_path,'records','kuband'};
input_dirs{end+1} = {gRadar.support_path,'frames','kuband'};
input_dirs{end+1} = {gRadar.support_path,'records','rds'};
input_dirs{end+1} = {gRadar.support_path,'frames','rds'};
input_dirs{end+1} = {gRadar.support_path,'records','snow'};
input_dirs{end+1} = {gRadar.support_path,'frames','snow'};
input_dirs{end+1} = {gRadar.support_path,'gps'};
input_dirs{end+1} = {gRadar.out_path,'kuband'};
input_dirs{end+1} = {gRadar.out_path,'rds'};
input_dirs{end+1} = {gRadar.out_path,'snow'};

backup_dir = '/net/field1/landing/scratch';

copyfile_filematch = '20161015';

% ========================================================================
% Automated Section
% ========================================================================

for input_idx = 1:length(input_dirs)
  input_dir = input_dirs{input_idx};
  
  in_fns = get_filenames(fullfile(input_dir{:}),'',copyfile_filematch,'',struct('recursive',true,'type','f'));

  for in_fn_idx = 1:length(in_fns)
    in_fn = in_fns{in_fn_idx};
    out_fn = fullfile(backup_dir,in_fn(length(input_dir{1})+1:end));
    out_fn_dir = fileparts(out_fn);
    if ~exist(out_fn_dir,'dir'); mkdir(out_fn_dir); end;
    copy_file = false;
    if ~exist(out_fn)
      copy_file = true;
    else
      in_finfo = dir(in_fn);
      out_finfo = dir(out_fn);
      if in_finfo.datenum ~= out_finfo.datenum ...
          || in_finfo.bytes ~= out_finfo.bytes
        copy_file = true;
      end
    end
    if copy_file
      fprintf('Copying %s\n  %s\n', in_fn, out_fn);
      copyfile(in_fn,out_fn)
    end
  end
end


