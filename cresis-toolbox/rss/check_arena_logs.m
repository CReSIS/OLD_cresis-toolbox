function logs = check_arena_logs(param,param_override)
% logs = check_arena_logs(param,param_override)
%
% Checks the contents of ARENA log files
%
% See run_check_arena_logs.m for an example of how to run
%
% Author: John Paden
%
% See also: run_check_arena_logs.m, check_arena_logs.m, read_arena_logs.m


%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Determine which files to load

[~,config_fn_name] = fileparts(param.records.config_fn);
config_timestamp = config_fn_name(1:15);
metadata_dir = param.day_seg(1:8);

log_dir = fullfile(param.data_support_path,param.season_name,metadata_dir);

fns = get_filenames(log_dir,config_timestamp,'','.txt');

logs = read_arena_logs(fns);

%% Print out notes

%% Check for issues
check_arena_logs_sub(logs);


end

function check_arena_logs_sub(logs,parent_str)

if ~exist('parent_str','var')
  parent_str = '';
end

log_fieldnames = fieldnames(logs);

voltage_threshold = 0.05;
temperature_threshold = 65;
alignment_threshold = 8;
link_error_threshold = 1;

for field_idx = 1:length(log_fieldnames)
  fieldname = log_fieldnames{field_idx};
  if 0
    % Debug
    fprintf('%s\n',fieldname);
  end
  if isstruct(logs.(fieldname))
    if isempty(parent_str)
      new_parent_str = fieldname;
    else
      new_parent_str = [parent_str '.' fieldname];
    end
    check_arena_logs_sub(logs.(fieldname), new_parent_str);
  else
    % Voltage field
    if ~isempty(regexp(fieldname,'V')) && isnumeric(logs.(fieldname))
      mean_val = mean(logs.(fieldname));
      std_val = std(logs.(fieldname));
      if std_val/mean_val > voltage_threshold
        fprintf(2,'%s.%s: mean(%.2g) std(%.2g)\n', parent_str, fieldname, mean_val, std_val);
      end
    end
    
    % Temperature field
    if ~isempty(regexpi(fieldname,'temp')) && isnumeric(logs.(fieldname))
      mean_val = mean(logs.(fieldname));
      max_val = max(logs.(fieldname));
      std_val = std(logs.(fieldname));
      if max_val > temperature_threshold
        fprintf(2,'%s.%s: mean(%.2g) max(%.2g) std(%.2g)\n', parent_str, fieldname, ...
          mean_val, max_val, std_val);
      end
    end
    
    % Alignment field
    if ~isempty(regexpi(fieldname,'alignment'))
      alignment = [];
      for idx=1:length(logs.(fieldname))
        tmp = textscan(logs.(fieldname){idx},'%d:%d:%d','CollectOutput',true);
        if numel(tmp{1}) < 3
          continue
        end
        if length(alignment) < tmp{1}(1)
          alignment{tmp{1}(1)}(1:2,1) = tmp{1}(2:3);
        else
          alignment{tmp{1}(1)}(1:2,end+1) = tmp{1}(2:3);
        end
      end
      for idx = 1:length(alignment)
        min_min = min(alignment{idx}(1,:));
        min_max = max(alignment{idx}(1,:));
        max_min = min(alignment{idx}(2,:));
        max_max = max(alignment{idx}(2,:));
        if abs(min_max-min_min) >= alignment_threshold ...
            || abs(max_max-max_min) >= alignment_threshold
          fprintf(2,'%s.%s: %d [%d to %d] : [%d to %d]\n', parent_str, fieldname, idx, ...
            min_min, min_max, ...
            max_min, max_max);
        end
      end
    end
    
    % Link error count field
    if ~isempty(regexpi(fieldname,'LinkErrorCount'))
      logs.(fieldname) = cellfun(@(x) x(3:end), logs.(fieldname), 'UniformOutput', false);
      logs.(fieldname) = hex2dec(logs.(fieldname));
      min_val = min(logs.(fieldname));
      max_val = max(logs.(fieldname));
      if max_val-min_val >= link_error_threshold
        fprintf(2,'%s.%s: min(%d) max(%d)\n', parent_str, fieldname, min_val, max_val);
      end
    end
     
    % Notes field
    if ~isempty(regexpi(fieldname,'notes'))
      fprintf('Notes:\n');
      for idx=1:length(logs.(fieldname))
        fprintf('%s: %s\n\n', datestr(epoch_to_datenum(logs.nmeaTime(idx))), logs.(fieldname){idx});
      end
    end   
    
  end
end

end
