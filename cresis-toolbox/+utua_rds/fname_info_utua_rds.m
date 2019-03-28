function fname = fname_info_utua_rds(fn)
% fname = fname_info_utua_rds(fn)
%
% Parses UTUA filename
%
% fn = UTUA RDS file name
%
% fname = structure containing each of the fields of the MCORDS
%   filename
%
%   For example: 20180523-233532.tdms or 20180523-233532.mat
%  .datenum = time stamp converted to Matlab date number. Common
%     usage: "[year month day hour min sec] = datevec(fname.datenum)"
%
% Example
%  fn = 'E:\ct_tmp\headers\rds\2018_Alaska_SO\20180523\20180523-233532.mat';
%  fname = utua_rds.fname_info_utua_rds(fn)
%
% Author: John Paden
%
% See also datenum

[path name ext] = fileparts(fn);

[day_str time_str] = strtok(name,'-');

fname.datenum = datenum([day_str time_str(2:end)],'yyyymmddHHMMSS');
