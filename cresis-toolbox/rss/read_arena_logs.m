function logs = read_arena_logs(fns)
% logs = read_arena_logs(fns)
%
% Reads ARENA log files
%
% fns = get_filenames('/data/20180927/logs/','','','.txt');
% logs = read_arena_logs(fns);
%
% log_idx = find(strcmp('20180905_164210',{logs.date_str}));
% for idx=1:length(logs(log_idx).notes.nmeaTime_str)
%   fprintf('%s: %s\n', logs(log_idx).notes.nmeaTime_str{idx}, logs(log_idx).notes.notes{idx});
% end
%
% Author: John Paden
%
% See also: run_check_arena_logs.m, check_arena_logs.m, read_arena_logs.m

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
  
  if 1
    %% Method 1
    fseek(fid,0,-1);
    C=textscan(fid,'%s%s%s%s%*[^\n]','delimiter',':');
    [log_fieldnames,unique_idxs,map_idxs] = unique(C{1});
    field_nums = str2double(C{2});
    str_mask = isnan(str2double(C{2}));
    empty_state3 = cellfun(@(x) ~isempty(x), C{3});
    empty_state4 = cellfun(@(x) ~isempty(x), C{4});
    
    new_log = [];
    for field_idx = 1:length(unique_idxs)
      field_name = log_fieldnames{field_idx};
      if empty_state3(unique_idxs(field_idx))
        if C{3}{unique_idxs(field_idx)}(1)=='C'
          if isempty(regexpi(field_name,'Temp'))
            field_name = [field_name 'Temp'];
          end
        elseif C{3}{unique_idxs(field_idx)}(1)=='V'
          if isempty(regexpi(field_name,'V'))
            field_name = [field_name 'V'];
          end
        else
          C{2}(map_idxs==field_idx) = cellfun(@(x,y) [x ':' y], C{2}(map_idxs==field_idx), C{3}(map_idxs==field_idx),'UniformOutput',false);
          str_mask(unique_idxs(field_idx)) = true;
        end
      end
      if empty_state3(unique_idxs(field_idx)) && C{3}{unique_idxs(field_idx)}(1)~='C' && C{3}{unique_idxs(field_idx)}(1)~='V'
        C{2}(map_idxs==field_idx) = cellfun(@(x,y) [x ':' y], C{2}(map_idxs==field_idx), C{4}(map_idxs==field_idx),'UniformOutput',false);
        str_mask(unique_idxs(field_idx)) = true;
      end
      field_name(~isstrprop(field_name,'alphanum')) = '_';
      if isstrprop(field_name(1),'digit')
        field_name = ['N' field_name];
      end
      if str_mask(unique_idxs(field_idx))
        new_log.(field_name) = C{2}(map_idxs==field_idx);
      else
        new_log.(field_name) = field_nums(map_idxs==field_idx);
      end
    end
    cmd_str = sprintf('try; logs(log_idx)%s = merge_structs(logs(log_idx)%s,new_log); catch; logs(log_idx)%s = new_log; end;',eval_str,eval_str,eval_str);
    eval(cmd_str);
    if 0
      % Debug
      fprintf('logs(log_idx)%s',eval_str)
      eval(sprintf('logs(log_idx)%s',eval_str))
      fprintf('\n');
      keyboard
    end
    
  else
    %% Method 2
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
  end
  
  % Close the file
  fclose(fid);
end
