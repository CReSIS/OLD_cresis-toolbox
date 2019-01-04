function convert_mat_version(fn)

if exist('fn','dir')
  fns = get_filenames(fn,'','','.mat')
  for fn_idx = 1:length(fns)
    fid = fns{fn_idx};
    convert_mat_version(fns)
  end
  return;
end

tmp = load(fn);
movefile(fn,tmp_dir);
save(fn,'-v7.3','-struct','tmp');


