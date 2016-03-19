function [params] = read_param_xls_generic(param_fn, generic_ws, params)
% [params] = read_param_xls_generic(param_fn, generic_ws, params)
%
% Support function for read_param_xls, generic worksheet reader.
%
% param_fn: MS Excel .xls parameter spreadsheet to load in
% generic_ws: worksheet name to load in
% params: add new fields to this structure
% 
% Author: Theresa Stumpf, John Paden
%
% See also: read_param_xls

% ======================================================================%
% CREATING THE PARAM STRUCTURE ARRAY FROM PARAM_STARTER.XLS
% ======================================================================%
warning('off','MATLAB:xlsread:Mode');

% =======================================================================
% Create Generic Parameters
% =======================================================================
sheet_name = generic_ws;
fprintf('Reading sheet %s of xls file: %s\n', sheet_name, param_fn);
warning off MATLAB:xlsfinfo:ActiveX
[status, sheets] = xlsfinfo(param_fn);
warning on MATLAB:xlsfinfo:ActiveX
sheet_idx = strmatch(generic_ws,sheets,'exact');
if isempty(sheet_idx)
  fprintf('  Sheet not found\n');
  return
end

[num txt] = xlsread(param_fn,sheet_name,'','basic');

num_header_rows = find(strcmp( 'Date' , txt(1:end,1) )) + 1;
if isempty(num_header_rows)
  error('Could not find required "Date" field in first column.');
end

rows = max(size(num,1), size(txt,1)) - num_header_rows;
field_names = txt(num_header_rows-1,:);
field_types = txt(num_header_rows,:);

missing_types = find(cellfun('isempty',field_types) & ~cellfun('isempty',field_names));
if ~isempty(missing_types)
  for col = missing_types
    fprintf('  Field %s (column %s/%d) is missing type definition in row 2\n',field_names{col}, char(65+mod(col-1,26)), col);
  end
  error('Type definitions are missing');
end

for idx = 1:rows
  row = idx + num_header_rows;
  day_seg = sprintf('%08.0f_%02.0f',num(row,1),num(row,2));
  if strcmpi(generic_ws,'cmd')
    params(idx).day_seg = day_seg;
  elseif (idx > length(params) || ~strcmpi(params(idx).day_seg,day_seg))
    error('The date segment order of sheets cmd and %s do not match at row %d in the excel spreadsheet. Each sheet must have the same date and segment list.',sheet_name,row);
  end
  for col = 3:length(field_names)
    try
    if ~isempty(field_names{col})
      period_idxs = find(field_names{col} == '.');
      
      if field_types{col} ~= 'a'
        field_type = field_types{col};
        is_array = false;
      elseif field_types{col}(1) == 'a'
        field_type = field_types{col}(2);
        is_array = true;
      else
        error('Unsupported field type (%s)', field_types{col});
      end
      if field_type == 'b'
        val = read_param_xls_boolean(row,col,num,txt);
      elseif field_type == 't'
        val = read_param_xls_text(row,col,num,txt);
      else
        val = read_param_xls_general(row,col,num,txt);
      end
      if ~is_array
        if length(period_idxs) == 0
          params(idx).(generic_ws).(field_names{col}) = val;
        elseif length(period_idxs) == 1
          params(idx).(generic_ws).(field_names{col}(1:period_idxs-1)).(field_names{col}(period_idxs+1:end)) = val;
        else
          error('read_param_xls_generic:toomanyperiods', ...
            ' Fieldnames %s with more than 1 period not supported', field_names{col})
        end
      else
        paren_idx = find(field_types{col}(4:end) == '(');
        if isempty(paren_idx)
          % Size of the struct array field (the size field should go first)
          array_field_name = field_types{col}(4:end);
          params(idx).(generic_ws).(array_field_name) = rmfield(struct('size',cell(1,val)),'size');
        else
          % Fields of the struct array
          array_field_name = field_types{col}(4+(0:paren_idx-2));
          array_field_idx = str2double(field_types{col}(4+paren_idx:end-1));
          
          if ~isfield(params(idx).(generic_ws),array_field_name) ...
              || length(params(idx).(generic_ws).(array_field_name)) >= array_field_idx
            % Only include the field if we don't exceed the declared
            % size of the array
            if length(period_idxs) == 0
              params(idx).(generic_ws).(array_field_name)(array_field_idx).(field_names{col}) = val;
            elseif length(period_idxs) == 1
              params(idx).(generic_ws).(array_field_name)(array_field_idx).(field_names{col}(1:period_idxs-1)).(field_names{col}(period_idxs+1:end)) = val;
            else
              error('read_param_xls_generic:toomanyperiods', ...
                ' Fieldnames %s with more than 1 period not supported', field_names{col})
            end
          else
            if ~isempty(val)
              fprintf('  Warning in row %d, column %s/%d (field %s)\n', row, char(65+mod(col-1,26)), col, field_names{col});
              fprintf('    Struct array %s should have %d elements, but has more (ignoring)\n', array_field_name, length(params(idx).(generic_ws).(array_field_name)));
            end
          end
        end
      end
    end
    catch ME
      if ~strcmp(ME.identifier,'read_param_xls_generic:toomanyperiods')
        fprintf('  Error in row %d, column %s/%d (field %s)\n', row, char(65+mod(col-1,26)), col, field_names{col});
      end
      rethrow(ME)
    end
  end
  
end

return

