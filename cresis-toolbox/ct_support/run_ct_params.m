% script run_ct_params
%
% A script for testing and using the class ct_param
% NOT FINISHED
% 
% Author: John Paden
%
% See also: ct_param.m


global obj;
if ~isempty(obj)
  delete(obj);
end
fn = ct_filename_param('accum_param_2018_Antarctica_TObas.xls')
obj = ct_params(fn);

