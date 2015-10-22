% script convert_matlab_files_to_v6
%
% Converts a whole directory of .MAT files to version 6 (no compression)
%
% Author: John Paden

fn_dir = '/cresis/scratch1/paden/CSARP_post';

fns = get_filenames(fn_dir,'','','.mat',struct('recursive',true));

convert_matlab_files_to_v6_tstart = tic;

for fn_idx = 1:length(fns)
  if mod(fn_idx-1,500) == 0
    fprintf('File %d of %d (%.1f sec)\n', fn_idx, length(fns), ...
      toc(convert_matlab_files_to_v6_tstart));
  end
  fn = fns{fn_idx};
  tmp = load(fn);
  save(fn,'-v6','-struct','tmp')
end


return;
