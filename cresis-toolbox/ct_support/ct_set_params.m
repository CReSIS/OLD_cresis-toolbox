function params = ct_set_params(params,field,value,filter_field,filter_regexp)
% params = ct_set_params(params,field,value,filter_field,filter_regexp)
%
% CReSIS toolbox utility function for setting parameter spreadsheet
% parameters in bulk. It allows a single field to be set for a range of
% segments rather than having to do it one at a time.
%
% params: struct array of parameters. Usually from read_param_xls/parameter
%   spreadsheet.
% field: a string containing the field to be set (e.g. 'cmd.generic')
% value: the value to assign to the field
% filter_field: Optional. Specifies a field to filter on (only struct
%   elements with a matching field will then have there field updated)
% filter_regexp: A regular expression that specifies the filtering to do
%   on filter_field.
%
% params: the updated struct array
%
% Example:
%  params = ct_set_params(params,'cmd.generic',0);
%  params = ct_set_params(params,'cmd.generic',1,'day_seg','20140401_03|20140307_11');
%  params = ct_set_params(params,'cmd.generic',0,'cmd.notes','Do not process');
%  params = ct_set_params(params,'cmd.generic',1,'cmd.notes','^((?!Do not process).)*$');
%
% Authors: John Paden
%
% See also: read_param_xls, ct_set_params, ct_filename_param

if exist('filter_field','var') && ~isempty(filter_field)
  use_filter = true;
else
  use_filter = false;
end

for param_idx = 1:length(params)
  
  param = params(param_idx);
  
  if use_filter
    if 1
      cmd = sprintf('filter_value = params(param_idx).%s;', filter_field);
      eval(cmd);
    else
      str = filter_field;
      [token,str] = strtok(str,'.');
      field_names = '';
      while ~isempty(token)
        field_names = sprintf('%s.(''%s'')',field_names,token);
        [token,str] = strtok(str,'.');
      end
      
      cmd = sprintf('filter_value = params(param_idx)%s;', field_names);
      eval(cmd);
    end
  end
  
  if ~use_filter || ~isempty(regexpi(filter_value,filter_regexp))
    
    if 1
      cmd = sprintf('params(param_idx).%s = value;', field);
      eval(cmd);
    else
      str = field;
      [token,str] = strtok(str,'.');
      field_names = '';
      while ~isempty(token)
        field_names = sprintf('%s.(''%s'')',field_names,token);
        [token,str] = strtok(str,'.');
      end
      
      cmd = sprintf('params(param_idx)%s = value;', field_names);
      eval(cmd);
    end
  end
end
