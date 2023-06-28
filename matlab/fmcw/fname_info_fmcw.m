function fname = fname_info_fmcw(fn)
% fname = fname_info_fmcw(fn)
%
% Detects the type of FMCW Radar filename and then parses it.
%
% fn = FMCW file name
%
% fname = structure containing each of the fields of the FMCW
%   filename
%
%   For example: data04.11022009.0050.dat
%  .name = 'data'
%  .board = 0
%  .group = '04'
%  .file_idx = 50
%  .ver = 1
%  .rx = 1
%  .datevec = time stamp converted to Matlab date number. Common
%     usage: "[year month day hour min sec] = datevec(fname.datenum)"
%
% See also datenum

fname.dir_info = dir(fn);

[path name ext] = fileparts(fn);
fn = [name ext];

if strcmpi(ext,'.bin')
  [fname.name fn] = strtok(fn,'_');
  [fname.group fn] = strtok(fn,'_');
  fname.group = str2double(fname.group);
  [day_str fn] = strtok(fn,'_');
  year = str2double(day_str(1:4));
  month = str2double(day_str(5:6));
  day = str2double(day_str(7:8));
  [time_str fn] = strtok(fn,'_');
  hour = str2double(time_str(1:2));
  minute = str2double(time_str(3:4));
  sec = str2double(time_str(5:end));
  fname.datenum = datenum(year,month,day,hour,minute,sec);
  [fname.file_idx] = str2double(strtok(fn(2:end),'.'));
  fname.board = 0;
elseif name(1:4) == 'data'
  [fname.name fn] = strtok(fn,'0123456789');
  [fname.group fn] = strtok(fn,'.');
  fname.group = str2double(fname.group);

  [time_stamp fn] = strtok(fn,'.');

  fname.datenum = datenum(str2double(time_stamp(5:8)), ...
    str2double(time_stamp(1:2)), str2double(time_stamp(3:4)), 0, 0, 0);

  fname.file_idx = strtok(fn,'.');
  fname.file_idx = str2double(fname.file_idx);

  fname.board = 0;
elseif strcmp(name(1:6),'kuband') || strcmpi(name(1:4),'snow')
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
elseif strcmpi(name(1:4),'fred')
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
