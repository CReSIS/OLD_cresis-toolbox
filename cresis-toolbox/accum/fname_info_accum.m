function fname = fname_info_accum(fn)
% fname = fname_info_accum(fn)
%
% Detects the type of accum filename and then parses it.
%
% fn = accum file name
%
% fname = structure containing each of the fields of the accum
%   filename
%
%   For example:
%     accum_20110329_10330963_0072.dat
%  .name = accum
%  .board = 0
%  .group = 1
%  .file_idx = 0
%  .datenum = time stamp converted to Matlab date number. Common
%     usage: "[year month day hour min sec] = datevec(fname.datenum)"
%
%     data.03292011.0105.dat
%
% Example
%  fn = '/mnt/tmp/20110329-accum/data.03292011.0105.dat';
%  fname = fname_info_accum(fn)
%
%  fn = '/mnt/tmp/20110329-accum/accum_20110329_10330963_0079.dat';
%  fname = fname_info_accum(fn)
%
% Author: John Paden
%
% See also datenum

[path name ext] = fileparts(fn);
fn = [name ext];

if name(1:4) == 'data'
  [fname.name fn] = strtok(fn,'.');
  [date_str fn] = strtok(fn,'.');
  year = str2double(date_str(5:8));
  month = str2double(date_str(1:2));
  day = str2double(date_str(3:4));
  fname.datenum = datenum(year,month,day,0,0,0);
  [fname.file_idx] = str2double(strtok(fn,'.'));
  fname.board = 0;
  fname.group = 0;
elseif name(1:5) == 'accum'
  [fname.name fn] = strtok(fn,'_');
  [day_str fn] = strtok(fn,'_');
  year = str2double(day_str(1:4));
  month = str2double(day_str(5:6));
  day = str2double(day_str(7:8));
  [time_str fn] = strtok(fn,'_');
  hour = str2double(time_str(1:2));
  minute = str2double(time_str(3:4));
  sec = str2double(time_str(5:end))/100;
  fname.datenum = datenum(year,month,day,hour,minute,sec);
  [fname.file_idx] = str2double(strtok(fn(2:end),'.'));
  fname.board = 0;
  fname.group = 0;
else
  error('Name type not supported');
end

return;
