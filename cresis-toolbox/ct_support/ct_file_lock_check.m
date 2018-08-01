function file_locked = ct_file_lock_check(fns,check_mode)
% file_locked = ct_file_lock_check(fns,check_mode)
%
% Throws an error if filename specified by fn is locked.
%
% fns: cell array of paths to files or a string containing path to file.
%   Each file should contain a 'file_version' field. If this field contains
%   an 'L', then an error is thrown.
% check_mode: 0, returns list of lock states, 1, throws an error as soon as
%   a locked file is found
%
% Author: John Paden
%
% See also: ct_filename_data, ct_filename_out, ct_filename_support,
%  ct_filename_tmp, ct_filename_gis

if ~exist('check_mode','var')
  check_mode = 1;
end

if ischar(fns)
  fn = fns;
  file_locked = false;
  if exist(fn,'file')
    [fn_dir,fn_file,fn_ext] = fileparts(fn);
    if strcmpi(fn_ext,'.nc')
      tmp = [];
      try
        file_version = ncread(fn,'file_version');
        tmp.file_version = file_version;
      end
    else
      warning off;
      tmp = load(fn,'file_version');
      warning on;
    end
    if isfield(tmp,'file_version')
      if any(tmp.file_version=='L')
        if check_mode
          error('File %s is locked.', fn);
        else
          file_locked = true;
        end
      end
    end
  end
else
  % fns should be a cell array of strings
  for idx=1:length(fns)
    file_locked(idx) = ct_file_lock_check(fns{idx}, check_mode);
  end
end
