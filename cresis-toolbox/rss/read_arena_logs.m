function logs = read_arena_logs(fns)
% logs = read_arena_logs(fns)
%
% Reads ARENA log files
%
% fns = get_filenames('/data/logs/','','','.txt');
% logs = read_arena_logs(fns);
%
% Author: John Paden

logs = struct(); logs = logs([]);

for fn_idx = 1:length(fns)
  
  fn = fns{fn_idx};
  
  [~,fn_name] = fileparts(fn);
  
  % Grab date string from filename
  date_str = fn_name(1:15);
  log_idx = 1;
  while log_idx <= length(logs)
    if strcmp(logs(log_idx).date_str,date_str)
      % Found existing log structure with matching date string
      break;
    end
    log_idx = log_idx + 1;
  end
  logs(log_idx).date_str = date_str;
  
  % Create a structure field evaluation string from the filename
  %   e.g. .ARENA.CTU.ctu.gps from 20180905_163632_ARENA-CTU-ctu-gps
  [token,remain] = strtok(fn_name(17:end),'-.');
  token(~isstrprop(token,'alphanum')) = '_';
  eval_str = sprintf('.%s',token);
  while ~isempty(remain)
    [token,remain] = strtok(remain,'-');
    eval_str = [eval_str sprintf('.%s',token)];
  end
  
  % Open the file
  [fid,msg] = fopen(fn,'rb');
  if fid == -1
    error(msg);
  end
  
  % Read each line of the file in
  while ~feof(fid)
    line_str = fgets(fid);
    line_str(line_str==10|line_str==13) = [];
    if isempty(line_str)
      continue;
    end
    [token,remain] = strtok(line_str,':');
    if length(remain) < 2
      continue;
    end
    remain = remain(2:end);
    token(~isstrprop(token,'alphanum')) = '_';
    if isstrprop(token(1),'digit')
      token = ['N' token];
    end
    token(~isstrprop(token(1),'alphanum')) = '_';
    
    val = str2double(remain);

    try
      cmd_str = sprintf('idx = 1+length(logs(log_idx)%s.%s);', eval_str,token);
      eval(cmd_str)
    catch
      idx = 1;
    end
    
    if isnan(val)
      try
        if length(remain) > 2 && strcmpi(remain(1:2),'0x')
          val = hex2dec(remain(3:end));
          cmd_str = sprintf('logs(log_idx)%s.%s(idx) = val;', eval_str,token);
        else
          error;
        end
      catch
        if length(remain)>2 && any(strcmpi(remain(end-1:end),{':C',':V'}))
          [val,~] = strtok(remain,':');
          val = str2double(val);
          cmd_str = sprintf('logs(log_idx)%s.%s(idx) = val;', eval_str,token);
        else
          cmd_str = sprintf('logs(log_idx)%s.%s{idx} = ''%s'';', eval_str,token,remain);
        end
      end
    else
      if isempty(regexpi(token,'time')) && ~strcmp(token,'ppsCntr')
        cmd_str = sprintf('logs(log_idx)%s.%s(idx) = val;', eval_str,token);
      else
        if strcmpi(token,'relTimeCntr')
          val = val / 10e6;
        end
        try
          cmd_str = sprintf('logs(log_idx)%s.%s(idx) = val;', eval_str,token);
          cmd_str = [cmd_str sprintf('logs(log_idx)%s.%s_str{idx} = ''%s'';', eval_str,token,datestr(epoch_to_datenum(val)))];
        catch
          cmd_str = sprintf('logs(log_idx)%s.%s(idx) = val;', eval_str,token);
        end
      end
    end
    if ~isempty(cmd_str)
      try
        eval(cmd_str);
      catch ME
        ME.getReport
        keyboard
      end
    end
    
    %     fprintf('%s\n',line_str);
  end
  
  % Close the file
  fclose(fid);
end
