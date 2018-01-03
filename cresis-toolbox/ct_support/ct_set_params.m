function params = ct_set_params(params,field,value,day_seg)
% params = ct_set_params(params,field,value,day_seg)
%
% CReSIS toolbox utility function for setting parameter spreadsheet
% parameters in bulk. It allows a single field to be set for a range of
% segments rather than having to do it one at a time.
%
% params: struct array of parameters. Usually from read_param_xls/parameter
%   spreadsheet.
% field: a string containing the field to be set (e.g. 'cmd.generic')
% value: the value to assign to the field
% day_seg: Optional. A regular expression that specifies which segments to
%   effect. The default is '.*' which is all segments.
%
% params: the updated struct array
%
% Example:
%  params = ct_set_params(params,'cmd.generic',0);
%  params = ct_set_params(params,'cmd.generic',1,'20140401_03|20140307_11');
%
% Authors: John Paden
%
% See also: read_param_xls, ct_set_params, ct_filename_param

if ~exist('day_seg','var') || isempty(day_seg)
  day_seg = '.*';
end

for param_idx = 1:length(params)
  
  param = params(param_idx);
  
  if ~isempty(regexpi(param.day_seg,day_seg))
    
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
