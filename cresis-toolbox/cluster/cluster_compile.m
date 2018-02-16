function cluster_compile(fun,hidden_depend_funs,force_compile,ctrl)
% cluster_compile(fun,hidden_depend_funs,force_compile,ctrl)
%
% Compiles functions to be used with the cluster_* functions. It tries
% to detect if any of the dependent functions to "fun" have changed
% and will recompile if they have.
%
% Makes use of "fdep.m" which is an external party function with separate
% license rules.
%
% MODE1:
% Inputs:
%  fun: specific function to be called by cluster_create_task
%    default is to check every function, specifying a function here allows
%    the date check part of cluster_compile to go faster
%  hidden_depend_funs: cell vector of hidden dependent functions cells
%    default is to use gRadar.cluster.hidden_depend_funs
%    Each cell in the cell vector has two fields:
%    {1}: filename
%    {2}: data check level, 0: do not check, 1: check always, 2: only
%          check when fun not specified
%  force_compile: force a compile even if dependent functions have not changed
%    default is true if ctrl is undefined or if ctrl.cluster.type is 'slurm' or
%    'torque'
% ctrl: ctrl structure returned from cluster_new_batch, default is empty
%  .cluster: cluster parameters
%    .type: (only used to determine default state of force_compile)
%    .mcc: string containing 'system' or 'eval' and specifies which of
%      these two functions will be used to execute the mcc command
%      system is generally preferred because it only requires the matlab
%      compiler license while the command line is running and then releases
%      it. Calling eval does not release the license until the matlab
%      session ends.
%
% Author: John Paden
%
% See also: cluster_chain_stage, cluster_cleanup, cluster_compile
%   cluster_exec_job, cluster_get_batch, cluster_get_batch_list, 
%   cluster_hold, cluster_job, cluster_new_batch, cluster_new_task,
%   cluster_print, cluster_run, cluster_submit_batch, cluster_submit_task,
%   cluster_update_batch, cluster_update_task

% Find dependencies that are directly called

if ~exist('ctrl','var')
  ctrl = [];
end

if ~isfield(ctrl,'cluster')
  ctrl.cluster = [];
end

if ~isfield(ctrl.cluster,'mcc') 
  ctrl.cluster.mcc = 'system';
end

if ~exist('fun','var')
  fun = {};
end

% Support the legacy format of a single string containing a function
if ischar(fun)
  fun = {fun};
end

if ~exist('hidden_depend_funs','var') || isempty(hidden_depend_funs)
  global gRadar;
  hidden_depend_funs = gRadar.cluster.hidden_depend_funs;
end

if isfield(ctrl.cluster,'type') && all(~strcmpi(ctrl.cluster.type,{'torque','slurm'}))
  return;
end  

if ~exist('force_compile','var')
  force_compile = 1;
end

matlab_ver = ver('matlab');
% 8.3 was first version of builtin, but it has bugs. So using 8.6 as
% threshold especially since fdep no longer works in that version
use_builtin_fdep = str2double(matlab_ver.Version) >= 8.6;

cluster_job_fn = fullfile(getenv('MATLAB_CLUSTER_PATH'),'cluster_job.m');
if ~force_compile
  % If any of the functions in depend_fun are newer, then recompile
  
  test_date = dir(fullfile(getenv('MATLAB_CLUSTER_PATH'),'run_cluster_job.sh'));
  
  if length(test_date) == 0
    % File does not exist, have to compile to make the file
    force_compile = true;
  else
    fun_info = dir(cluster_job_fn);
    if fun_info.datenum > test_date.datenum
      force_compile = true;
    end
    
    if ~isempty(fun)
      for input_idx = 1:length(fun)
        try
          if use_builtin_fdep
            flist = matlab.codetools.requiredFilesAndProducts(fun{input_idx});
          else
            warning off
            flist = fdep(fun{input_idx},'-q');
            warning on
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
              warning off
              flist = fdep(fn,'-q');
              warning on
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
  cluster_job_fn_dir = fileparts(cluster_job_fn);
  cmd = sprintf('mcc -m -d %s -R ''-singleCompThread,-nodisplay'' %s', cluster_job_fn_dir, cluster_job_fn);
  
  for dep_idx = 1:length(hidden_depend_funs)
    cmd = [cmd ' ' hidden_depend_funs{dep_idx}{1}];
  end
  
  % Check to make sure the working directory is not a package or class
  % directory. This messes up the compiler if it uses functions from those
  % packages or classes.
  working_dir = pwd;
  [~,working_dir_name] = fileparts(working_dir);
  if ~isempty(working_dir_name) && (working_dir_name(1) == '+' || working_dir_name(1) == '@')
    warning('Before dbcont, change directory out of any package or class directories, because this can cause mcc to fail. (E.g. "cd ..")');
    keyboard
  end
  
  fprintf('Start Compiling %s\n\n', datestr(now));
  fprintf('  %s\n', cmd);
  if strcmpi(ctrl.cluster.mcc,'system')
    system(cmd);
  else
    eval(cmd);
  end
  fprintf('\nDone Compiling %s\n', datestr(now));
end

return;
