function [param,data] = emquest_txt_reader(fn)

fid = fopen(fn,'r');

% Read param header
line = 0;
[param,line] = emquest_txt_reader_param(fid, line);

% Read data
format_line = fgets(fid); line = line + 1;
if ~strcmp(sprintf('Format\r\n'),format_line)
  warning('%d: Should be "Format", was %s', line, format_line);
end

while ~feof(fid)
  fline = fgets(fid); line = line + 1;
  wspace = isstrprop(fline, 'wspace');
  indentation = find(~wspace,1);
  fline = fline(indentation:end-2);
  
  [polarization,remain] = strtok(fline,' ');
  [data_type,remain] = strtok(remain(2:end),sprintf('\r\n'));
  data_type(data_type==' ') = '_';
  
  if ~any(strcmpi(data_type,{'Log_Magnitude','Phase'}))
    return;
  end
  fprintf('Reading %s\n', data_type);
  
  fline = fgets(fid); line = line + 1;
  wspace = isstrprop(fline, 'wspace');
  indentation = find(~wspace,1);
  fline = fline(indentation:end-2);
  
  [axis,remain] = strtok(fline,sprintf(' \r\n'));
  [units,remain] = strtok(remain,sprintf(' \r\n'));
  
  if strcmp(units,'(MHz)')
    axis_scale = 1e6;
  else
    keyboard
  end
  
  data.(polarization).(data_type).freq = []; freq_idx = 0;
  data.(polarization).(data_type).angle = [];
  data.(polarization).(data_type).val = zeros(2000,1000); % Preallocate?
  
  axis1_indentation = indentation;
  
  while ~feof(fid)
    cur_pos = ftell(fid);
    fline = fgets(fid);
    wspace = isstrprop(fline, 'wspace');
    indentation = find(~wspace,1);
    if indentation ~= axis1_indentation
      fseek(fid,cur_pos,-1);
      break;
    end
    line = line + 1;
    fline = fline(indentation:end-2);
    
    freq_idx = freq_idx + 1;
    data.(polarization).(data_type).freq(freq_idx) = str2double(fline) * axis_scale;
    
    %fprintf('Reading %.1f MHz\n',
    %data.(polarization).(data_type).freq(freq_idx)/1e6); % DEBUG
    
    fline = fgets(fid); line = line + 1;
    wspace = isstrprop(fline, 'wspace');
    indentation = find(~wspace,1);
    fline = fline(indentation:end-2);
    
    [axis2,remain] = strtok(fline,sprintf(','));
    [units2,remain] = strtok(remain(2:end),sprintf('\r\n'));
    
    axis2_indentation = indentation;
    angle_idx = 0;
    while ~feof(fid)
      
      cur_pos = ftell(fid);
      fline = fgets(fid);
      wspace = isstrprop(fline, 'wspace');
      indentation = find(~wspace,1);
      if indentation ~= axis2_indentation
        fseek(fid,cur_pos,-1);
        break;
      end
      line = line + 1;
      fline = fline(indentation:end-2);
      
      [pat_angle,remain] = strtok(fline,sprintf(','));
      [val,remain] = strtok(remain(2:end),sprintf('\r\n'));
      
      angle_idx = angle_idx + 1;
      data.(polarization).(data_type).angle(angle_idx) = str2double(pat_angle);
      data.(polarization).(data_type).val(freq_idx, angle_idx) = str2double(val);
    end
    
  end
  data.(polarization).(data_type).val = data.(polarization).(data_type).val(1:freq_idx,1:angle_idx);
  
end

fclose(fid);

end

function [cell_out,line] = emquest_txt_reader_param(fid,line)

cell_out = {};
cell_idx = 1;
cell_out{cell_idx} = struct();

fline = fgets(fid); line = line + 1;
wspace = isstrprop(fline, 'wspace');
indentation = find(~wspace,1);
fline = fline(indentation:end-2);
if length(fline) == 80 && all(fline == '-')
  %fprintf('%d: Found starting dash\n', line); % DEBUG
else
  error('%d: Format error', line);
end

found_end_dash = false;
while ~found_end_dash
  fline = fgets(fid); line = line + 1;
  wspace = isstrprop(fline, 'wspace');
  indentation = find(~wspace,1);
  fline = fline(indentation:end-2);
  if length(fline) == 80 && all(fline == '-')
    %fprintf('%d: Found ending dash\n', line); % DEBUG
    % Check to see if next line is at the same indentation
    cur_pos = ftell(fid);
    fline = fgets(fid);
    fseek(fid,cur_pos,-1);
    wspace = isstrprop(fline, 'wspace');
    next_indentation = find(~wspace,1);
    if next_indentation == indentation+5
      % This is a struct array
      cell_idx = cell_idx + 1;
      cell_out{cell_idx} = struct();
    else
      found_end_dash = true;
    end
  else
    end_name_idx = find(fline == '=',1) - 1;
    name = fline(1:end_name_idx);
    val = fline(end_name_idx+2:end);
    %fprintf('%d: %s = %s\n', line, name, val); % DEBUG
    if isempty(val)
      % Check next line to see if it is dash line
      cur_pos = ftell(fid);
      fline = fgets(fid);
      fseek(fid,cur_pos,-1);
      wspace = isstrprop(fline, 'wspace');
      next_indentation = find(~wspace,1);
      fline = fline(next_indentation:end-2);
      if indentation == next_indentation && fline(1) == '-'
        % If it is a dash line at the same indentation then it is a struct
        % and not an empty val
        [val,line] = emquest_txt_reader_param(fid, line);
      end
    end
    name = name(name ~= '_');
    cell_out{cell_idx}.(name) = val;
  end
end

end




