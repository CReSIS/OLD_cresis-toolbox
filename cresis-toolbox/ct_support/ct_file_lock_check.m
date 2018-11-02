function file_locked = ct_file_lock_check(fns,check_mode)
% file_locked = ct_file_lock_check(fns,check_mode)
%
% Throws an error if filename specified by fn is locked.
%
% fns: cell array of paths to files or a string containing path to file.
%   Each file should contain a 'file_version' field. If this field contains
%   an 'L', then an error is thrown.
% check_mode: Integer specifying one of these modes:
%   0: returns list of lock states
%   1: throws an error as soon as a locked file is found unless
%      gRadar.ct_file_lock_check == false.
%   2: If gRadar.ct_file_lock_check == true and locked file, then ask the
%      user if the file lock should be removed or throw and error.
%   3: Same as 2 except also set the file_version to 'D' for deletable if
%      the file is not locked or if the user asks to remove the lock.
%      * Used by cluster monitoring programs to mark output files that
%      already exist from a previous run.
%   4: Checks for existence of file. If it exists and does not have 'D'
%      in file_version then the function returns true.
%      * Used by cluster monitoring programs to check success criteria.
%
% Author: John Paden
%
% See also: ct_filename_data, ct_filename_out, ct_filename_support,
%  ct_filename_tmp, ct_filename_gis

if ~exist('check_mode','var')
  check_mode = 1;
end

global gRadar;
if isfield(gRadar,'ct_file_lock_check') && ~gRadar.ct_file_lock_check
  no_stdio = true;
  if check_mode ~= 3 && check_mode ~= 4
    check_mode = 0;
  end
else
  no_stdio = false;
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
      try
        tmp = load(fn,'file_version');
      catch
        % Corrupt file so consider it to not be locked
        return;
      end
      warning on;
    end
    if isfield(tmp,'file_version')
      file_version = tmp.file_version;
    else
      file_version = '1';
    end
    switch (check_mode)
      case 0
        % Do nothing except return the file_locked state
        if any(file_version=='L')
          file_locked = true;
        end
        
      case 1
        if any(file_version=='L')
          error('File %s is locked.', fn);
        end
        
      case {2,3}
        if any(file_version=='L')
          % Check with user
          if no_stdio
            uinput = 1;
          else
            fprintf('<strong>File is locked: %s</strong>\nChoose one of these options:\n  1: Remove lock on this file\n  2: Disable gRadar.ct_file_lock_check (which disabled file lock checking globally)\n  3: Stop execution\n', fn);
            uinput = [];
            while isempty(uinput) || ~isnumeric(uinput)
              uinput = input('? ');
            end
          end
          if uinput==1
            if ~no_stdio
              fprintf('  Removing lock\n');
            end
            if check_mode == 2
              % Remove lock
              ct_file_lock(fn,0);
            else
              % Set file to delete mode
              ct_file_lock(fn,2);
            end
          elseif uinput==2
            fprintf('  Removing lock and setting gRadar.ct_file_lock_check=0\n');
            gRadar.ct_file_lock_check = 0;
            if check_mode == 2
              % Remove lock
              ct_file_lock(fn,0);
            else
              % Set file to delete mode
              ct_file_lock(fn,2);
            end
          else
            error('File %s is locked.', fn);
          end
        end

      case 4
        if ~any(file_version=='D')
          % File exists and is not in deleted state, return success
          file_locked = true;
        end
        
    end
  end
else
  % fns should be a cell array of strings
  for idx=1:length(fns)
    file_locked(idx) = ct_file_lock_check(fns{idx}, check_mode);
  end
end
