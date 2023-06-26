function fname = fname_info_mcords(fn)
% fname = fname_info_mcords(fn)
%
% Detects the type of MCORDS filename and then parses it.
%
% fn = MCORDS file name
%
% fname = structure containing each of the fields of the MCORDS
%   filename
%
%   For example: mcords.rec003.0001.r4-2.20100512132256.dat
%  .name = mcords
%  .group = 'rec003'
%  .file_idx = 1
%  .radar_num = 4
%  .rx = 2
%  .datevec = time stamp converted to Matlab date number. Common
%     usage: "[year month day hour min sec] = datevec(fname.datenum)"
%
%   For example: r1-1.20091102092921.0000.dat
%  .name = ''
%  .group = ''
%  .file_idx = 0
%  .radar_num = 1
%  .rx = 1
%  .datevec = time stamp converted to Matlab date number. Common
%     usage: "[year month day hour min sec] = datevec(fname.datenum)"
%
% See also datenum

[path name ext] = fileparts(fn);
fn = [name ext];

rev_idx = regexp(fn,'r\d+-');

if rev_idx == 1
  % File naming convention version 1
  fname.rev = 1;
  
  fname.name = '';
  fname.group = '';
  [fname.radar_num fn_rem] = strtok(fn(2:end),'-');
  fname.radar_num = str2double(fname.radar_num);
  [fname.rx fn_rem] = strtok(fn_rem(2:end),'.');
  fname.rx = str2double(fname.rx);
  
  [time_stamp fn_rem] = strtok(fn_rem,'.');

  fname.datenum = datenum(str2double(time_stamp(1:4)), ...
    str2double(time_stamp(5:6)), str2double(time_stamp(7:8)), ...
    str2double(time_stamp(9:10)), str2double(time_stamp(11:12)), ...
    str2double(time_stamp(13:14)));
  
  fname.file_idx = strtok(fn_rem,'.');
  fname.file_idx = str2double(fname.file_idx);
  
elseif rev_idx == 20
  % File naming convention version 2
  fname.rev = 2;
  
  [fname.name fn_rem] = strtok(fn,'.');
  [fname.group fn_rem] = strtok(fn_rem,'.');
  [fname.file_idx fn_rem] = strtok(fn_rem,'.');
  fname.file_idx = str2double(fname.file_idx);
  [fname.radar_num fn_rem] = strtok(fn_rem(2:end),'-');
  fname.radar_num = str2double(fname.radar_num(2:end));
  [fname.rx fn_rem] = strtok(fn_rem(2:end),'.');
  fname.rx = str2double(fname.rx);
  [time_stamp] = strtok(fn_rem,'.');

  fname.datenum = datenum(str2double(time_stamp(1:4)), ...
    str2double(time_stamp(5:6)), str2double(time_stamp(7:8)), ...
    str2double(time_stamp(9:10)), str2double(time_stamp(11:12)), ...
    str2double(time_stamp(13:14)));

elseif rev_idx == 15
  % File naming convention version 3
  fname.rev = 3;
  
  [fname.name fn_rem] = strtok(fn,'.');
  [fname.group fn_rem] = strtok(fn_rem,'.');
  [fname.radar_num fn_rem] = strtok(fn_rem(2:end),'-');
  fname.radar_num = str2double(fname.radar_num(2:end));
  [fname.rx fn_rem] = strtok(fn_rem(2:end),'.');
  fname.rx = str2double(fname.rx);
  [time_stamp fn_rem] = strtok(fn_rem,'.');
  [fname.file_idx] = strtok(fn_rem,'.');
  fname.file_idx = str2double(fname.file_idx);

  fname.datenum = datenum(str2double(time_stamp(1:4)), ...
    str2double(time_stamp(5:6)), str2double(time_stamp(7:8)), ...
    str2double(time_stamp(9:10)), str2double(time_stamp(11:12)), ...
    str2double(time_stamp(13:14)));
else
  error('Unsupported file version\n');
  
end

return;
