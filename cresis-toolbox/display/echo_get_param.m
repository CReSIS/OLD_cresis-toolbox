function param = echo_get_param(mdata)
% param = echo_get_param(mdata)
%
% Returns the processing parameters structure used to create the data
% mdata. mdata is loaded from qlook.m or array.m echogram outputs.
%
% mdata: structure loaded from echogram file (e.g. CSARP_qlook or
% CSARP_standard)
%
% fn = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_standard/20140512_01/Data_20140512_01_018.mat';
% mdata = load(fn);
% param = echo_get_param(mdata);
%
% Author: John Paden
%

if isfield(mdata,'param_array')
  param = mdata.param_array;
elseif isfield(mdata,'param_qlook')
  param = mdata.param_qlook;
else
  error('There is no param_array or param_qlook field in mdata.');
end
