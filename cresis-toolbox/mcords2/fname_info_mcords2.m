function fname = fname_info_mcords2(fn)
% fname = fname_info_mcords2(fn)
%
% Detects the type of MCORDS filename and then parses it.
%
% fn = MCORDS2 file name
%
% fname = structure containing each of the fields of the MCORDS
%   filename
%
%   For example: mcords_0_03012011_171839_01_0000.bin
%  .name = mcords
%  .board = 0
%  .group = 1
%  .file_idx = 0
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
% Author: John Paden
%
% See also datenum

[path name ext] = fileparts(fn);
fn = [name ext];

[fname.name fn] = strtok(fn,'_');

[fname.board fn] = strtok(fn,'_');
fname.board = str2double(fname.board);

[day_str fn] = strtok(fn,'_');
[time_str fn] = strtok(fn,'_');
fname.datenum = datenum([day_str time_str],'yyyymmddHHMMSS');

[fname.group fn] = strtok(fn,'_');
fname.group = str2double(fname.group);

[fname.file_idx fn] = strtok(fn(2:end),'.');
fname.file_idx = str2double(fname.file_idx);
