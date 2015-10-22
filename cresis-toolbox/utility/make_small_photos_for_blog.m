fn_dir = 'C:\Users\dangermo\Documents\Travel\Antarctica_2013\BLOG\combined';

[tmp fn_dir_name] = fileparts(fn_dir);
out_fn_dir = fullfile(fn_dir,sprintf('%s_lowres', fn_dir_name));
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end

fns = get_filenames(fn_dir,'','','*.jpg');

max_file_size = 250e3;

for idx = 1:length(fns)
  fn = fns{idx};
  [tmp fn_name fn_ext] = fileparts(fn);
  out_fn = fullfile(out_fn_dir, [fn_name fn_ext]);
  fprintf('%s to %s\n', fn, out_fn);
  A = imread(fn);
  A = imresize(A,0.5);
  imagesc(A);
  axis equal;
  pause
  done = false;
  quality = 80;
  while ~done
    imwrite(A,out_fn,'Quality',quality);
    finfo = dir(out_fn);
    if finfo.bytes > max_file_size && quality > 0
      quality = quality - 1;
    else
      done = true;
    end
  end
end

return;
