function fname = fname_info_accum2(fn)
% fname = fname_info_accum2(fn)
%
% Detects the type of accumulation filename and then parses it.
%
% fn = accum2 file name
%
% fname = structure containing each of the fields of the accumulation
%   filename
%
%   For example: mcords_0_03012011_171839_01_0000.bin
%   For example: accum2_00_20120317_183859_0011.bin
%  .name = accum2
%  .group = 0
%  .file_idx = 11
%  .datenum = time stamp converted to Matlab date number. Common
%     usage: "[year month day hour min sec] = datevec(fname.datenum)"
%
% Example
%  fn = '/cresis/scratch1/paden/mcords2/board_0/mcords_0_03012011_171839_01_0000.bin';
%  fname = fname_info_mcords2(fn)
%
%  fn = '/cresis/scratch1/paden/mcords2/board_0/mcords_1_03012011_171946_01_0009.bin';
%  fname = fname_info_mcords2(fn)
%
% Author: John Paden, Logan Smith
%
% See also datenum

[path name ext] = fileparts(fn);
fn = [name ext];

[fname.name fn] = strtok(fn,'_');

[fname.group fn] = strtok(fn,'_');
fname.group = str2double(fname.group);

[day_str fn] = strtok(fn,'_');
[time_str fn] = strtok(fn,'_');
fname.datenum = datenum([day_str time_str],'yyyymmddHHMMSS');

[fname.file_idx fn] = strtok(fn(2:end),'.');
fname.file_idx = str2double(fname.file_idx);

return;
