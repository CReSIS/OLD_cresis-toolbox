function fname = fname_info_mcrds(fn)
% fname = fname_info_mcrds(fn)
%
% Detects the type of MCRDS filename and then parses it.
%
% fn = MCRDS file name
%
% fname = structure containing each of the fields of the MCORDS
%   filename
%
%   For example: data.20080730212231.0500.raw
%  .name = data
%  .file_idx = 500
%  .datenum = time stamp converted to Matlab date number. Common
%     usage: "[year month day hour min sec] = datevec(fname.datenum)"
%
% Example
%  fn = '/cresis/data1/MCRDS/2008_Greenland/data.20080730212231.0500.raw';
%  fname = fname_info_mcrds(fn)
%
% Author: John Paden
%
% See also datenum

[path name ext] = fileparts(fn);
fn = [name ext];

[fname.name fn] = strtok(fn,'.');

[date_str fn] = strtok(fn,'.');
fname.datenum = datenum(date_str,'yyyymmddHHMMSS');

[fname.file_idx fn] = strtok(fn(2:end),'.');
fname.file_idx = str2double(fname.file_idx);

return;
