function read_param_xls_print(param_fn, generic_ws, params)
% read_param_xls_print(param_fn, generic_ws, params)
%
% Support function for read_param_xls, generic worksheet reader.
%
% param_fn: MS Excel .xls parameter spreadsheet to load in
% generic_ws: worksheet name to load in
% params: add new fields to this structure
%
% Author: John Paden
%
% See also: read_param_xls

%% Setup
% ======================================================================%
warning('off','MATLAB:xlsread:Mode');

% Get the headers from the parameter spreadsheet
[headers] = read_param_xls_headers(param_fn, generic_ws);

%% Print headers
for idx = 1:length(params)
  % Print field names
  fprintf('Day\t1');
  for idx = 1:length(headers)
    fprintf('\t%s', headers(idx).field_names);
  end
  fprintf('\n');

  % Print field types
  fprintf('YYYYMMDD\tSegment');
  for idx = 1:length(headers)
    fprintf('\t%s', headers(idx).field_types);
  end
  fprintf('\n');
end

%% Print values
for param_idx = 1:length(params)
  % Print field values
  fprintf('%s\t%d', params(param_idx).day_seg(1:8), str2double(params(param_idx).day_seg(10:11)));
  for h_idx = 1:length(headers)
    if headers(h_idx).array_type == '0'
      eval_str = sprintf('params(param_idx).%s.%s',headers(h_idx).generic_ws,headers(h_idx).field_names);
      value = eval(eval_str);
      if isempty(value)
        fprintf('\t');
      elseif strcmpi(headers(h_idx).field_types,'t')
        fprintf('\t%s',value);
      elseif strcmpi(headers(h_idx).field_types,'r')
        fprintf('\t%s',mat2str_generic(value));
      elseif strcmpi(headers(h_idx).field_types,'b')
        if (value)
          fprintf('\t1');
        else
          fprintf('\t');
        end
      end
      
    elseif headers(h_idx).array_type == 'a'
      eval_str = sprintf('numel(params(param_idx).%s.%s)',headers(h_idx).generic_ws,headers(h_idx).array_field_name);
      value = eval(eval_str);
      if headers(h_idx).array_field_idx <= value
        eval_str = sprintf('params(param_idx).%s.%s(%d).%s',headers(h_idx).generic_ws,headers(h_idx).array_field_name,headers(h_idx).array_field_idx,headers(h_idx).field_names);
        value = eval(eval_str);
        if isempty(value)
          fprintf('\t');
        elseif strcmpi(headers(h_idx).field_types(2),'t')
          fprintf('\t%s',value);
        elseif strcmpi(headers(h_idx).field_types(2),'r')
          fprintf('\t%s',mat2str_generic(value));
        elseif strcmpi(headers(h_idx).field_types(2),'b')
          if (value)
            fprintf('\t1');
          else
            fprintf('\t');
          end
        end
      end
      
    elseif headers(h_idx).array_type == 'as'
      eval_str = sprintf('numel(params(param_idx).%s.%s)',headers(h_idx).generic_ws,headers(h_idx).array_field_name);
      value = eval(eval_str);
      fprintf('\t%d',value);

    elseif headers(h_idx).array_type == 'c'
      
    elseif headers(h_idx).array_type == 'cs'
      eval_str = sprintf('numel(params(param_idx).%s.%s)',headers(h_idx).generic_ws,headers(h_idx).array_field_name);
      value = eval(eval_str);
      fprintf('\t%d',value);

    end
    
  end
  fprintf('\n');
end
