function convert_mat_version(fn,file_ver)
% convert_mat_version(fn,file_ver)
%
% Convert mat files to specified version. Default file_ver is '-v7.3'.
% fn can be a directory in which case all .mat files in that directory will
% be converted.
%
% Example:
% convert_mat_version('/cresis/snfs1/data/MCoRDS/BAS/20170122/');

if nargin < 2
  file_ver = '-v7.3';
end

if exist(fn,'dir')
  fns = get_filenames(fn,'','','.mat',struct('recursive',true));
  for fn_idx = 1:length(fns)
    fn = fns{fn_idx};
    convert_mat_version(fn,file_ver)
  end
  return;
end

tmp = load(fn);
fn_tmp = fullfile(tempdir,fn);
fprintf('%s\n  => %s\n', fn, fn_tmp);
movefile(fn,fn_tmp);
fprintf('  %s: %s\n', file_ver, fn);
save(fn,file_ver,'-struct','tmp');
