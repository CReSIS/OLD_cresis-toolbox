function ct_file_lock(fns,lock_state)
% ct_file_lock(fns,lock_state)
%
% Locks or unlocks a file or a cell array of files
%
% fns: cell array of paths to files or a string containing path to file.
%   Each file should contain a 'file_version' field.
%
% Author: John Paden
%
% See also: ct_filename_data, ct_filename_out, ct_filename_support,
%  ct_filename_tmp, ct_filename_gis

if ischar(fns)
  fn = fns;
  if ~exist(fn,'file')
    warning('File does not exist: %s', fn);
  else
    tmp = load(fn,'file_version');
    if isfield(tmp,'file_version')
      file_version = tmp.file_version(~isletter(tmp.file_version));
      if lock_state
        file_version = [file_version, 'L'];
      end
      save(fn,'-append','file_version');
    else
      file_version = '1';
      if lock_state
        file_version = [file_version, 'L'];
      end
      save(fn,'-append','file_version');
    end
  end
else
  % fns should be a cell array of strings
  for idx=1:length(fns)
    ct_file_lock(fns{idx},lock_state);
  end
end
