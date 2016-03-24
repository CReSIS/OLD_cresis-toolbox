function torque_compile(fun,hidden_depend_funs,force_compile)
% torque_compile(fun,hidden_depend_funs,force_compile)
%
% Compiles functions to be used with the torque_* functions. It tries
% to detect if any of the dependent functions to "fun" have changed
% and will recompile if they have.
%
% Makes use of "fdep.m" which is an external party function with separate
% license rules.
%
% MODE1:
% Inputs:
%  fun: specific function to be called by torque_create_task
%    default is to check every function, specifying a function here allows
%    the date check part of torque_compile to go faster
%  hidden_depend_funs: cell vector of hidden dependent functions cells
%    default is to use gRadar.sched.hidden_depend_funs
%    Each cell in the cell vector has two fields:
%    {1} = filename
%    {2} = data check level, 0: do not check, 1: check always, 2: only
%          check when fun not specified
%  force_compile: force a compile even if dependent functions have not changed
%    default is true
%
% Author: John Paden
%
% See also: torque_batch_list torque_cleanup torque_compile ...
%   torque_create_task torque_hold torque_job_list torque_job_status ...
%   torque_new_batch torque_print torque_rerun

% Find dependencies that are directly called

if ~exist('fun','var')
  fun = [];
end

if ~exist('hidden_depend_funs','var') || isempty(hidden_depend_funs)
  global gRadar;
  hidden_depend_funs = gRadar.sched.hidden_depend_funs;
end

if ~exist('force_compile','var')
  force_compile = 1;
end

matlab_ver = ver('matlab');
% 8.3 was first version of builtin, but it has bugs. So using 8.6 as
% threshold especially since fdep no longer works in that version
use_builtin_fdep = str2double(matlab_ver.Version) >= 8.6;

worker_fn = fullfile(getenv('MATLAB_TORQUE_PATH'),'worker_task.m');
if ~force_compile
  % If any of the functions in depend_fun are newer, then recompile
  
  test_date = dir(fullfile(getenv('MATLAB_TORQUE_PATH'),'run_worker_task.sh'));
  
  if length(test_date) == 0
    % File does not exist, have to compile to make the file
    force_compile = true;
  else
    fun_info = dir(worker_fn);
    if fun_info.datenum > test_date.datenum
      force_compile = true;
    end
    
    if ~isempty(fun)
      try
        if use_builtin_fdep
          flist = matlab.codetools.requiredFilesAndProducts(fun);
        else
          flist = fdep(fun,'-q');
          flist = flist.fun;
        end
        for fun_idx = 1:length(flist)
          fun_info = dir(flist{fun_idx});
          if fun_info.datenum > test_date.datenum
            force_compile = true;
            break;
          end
        end
      catch
        force_compile = true;
      end
    end
    
    if ~force_compile
      for fun_idx = 1:length(hidden_depend_funs)
        fn = hidden_depend_funs{fun_idx}{1};
        level = hidden_depend_funs{fun_idx}{2};
        if level == 1 || level == 2 && isempty(fun)
          try
            if use_builtin_fdep
              flist = matlab.codetools.requiredFilesAndProducts(fn);
            else
              flist = fdep(fn,'-q');
              flist = flist.fun;
            end
            for fun_idx = 1:length(flist)
              fun_info = dir(flist{fun_idx});
              if fun_info.datenum > test_date.datenum
                force_compile = true;
                break;
              end
            end
          catch
            force_compile = true;
            break;
          end
        end
      end
    end
    
  end
end

if force_compile
  worker_fn_dir = fileparts(worker_fn);
  cmd = sprintf('mcc -m -d %s -R ''-singleCompThread,-nodisplay'' %s', worker_fn_dir, worker_fn);
  
  for dep_idx = 1:length(hidden_depend_funs)
    cmd = [cmd ' ' hidden_depend_funs{dep_idx}{1}];
  end
  
  fprintf('Start Compiling %s\n\n', datestr(now));
  fprintf('  %s\n', cmd);
  system(cmd);
  fprintf('\nDone Compiling %s\n', datestr(now));
end

return;
