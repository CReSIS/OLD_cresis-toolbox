function ct_file_lock(fns,lock_state,file_version)
% ct_file_lock(fns,lock_state,file_version)
%
% Locks or unlocks a file or a cell array of files
%
% fns: cell array of paths to files or a string containing path to file.
%   Each file should contain a 'file_version' field.
% lock_state: one of three states:
%   0: file is not locked (no letter in file_version)
%   1: file is locked ('L' in file_version)
%   2: file should be deleted ('D' in file_version)
% file_version: in case the file does not have a file_version field, use
% this value, default values is '1'
%
% Author: John Paden
%
% See also: ct_filename_data, ct_filename_out, ct_filename_support,
%  ct_filename_tmp, ct_filename_gis

if ischar(fns)
  fn = fns;
  if exist(fn,'dir')
    % Lock/unlock all files in the directory
    fns = get_filenames(fn,'','','', ...
      struct('recursive',0,'regexp','.mat$|.nc$'));
    ct_file_lock(fns,lock_state);
  elseif ~exist(fn,'file')
    warning('File does not exist: %s', fn);
  else
    tmp = load(fn,'file_version');
    if isfield(tmp,'file_version')
      file_version = tmp.file_version(isstrprop(tmp.file_version,'digit'));
    elseif ~exist('file_version','var')
      % If file_version was not passed in, assume "1"
      file_version = '1';
    end
    if lock_state==1
      file_version = [file_version, 'L'];
    elseif lock_state==2
      file_version = [file_version, 'D'];
    end
    file_version = sprintf('%8s',file_version);
    save(fn,'-append','file_version');
  end
else
  % fns should be a cell array of strings
  for idx=1:length(fns)
    ct_file_lock(fns{idx},lock_state);
  end
end
