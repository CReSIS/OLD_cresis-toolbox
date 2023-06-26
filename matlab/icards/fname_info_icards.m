function fname = fname_info_icards(fn)
% fname = fname_info_icards(fn)
%
% Detects the type of ICARDS filename and then parses it.
%
% fn = ICARDS file name
%
% fname = structure containing each of the fields of the ICARDS
%   filename
%
%   For example: may19_01.009
%  .name = may19_01
%  .file_idx = 009
%
% Example
%  fn = 'S:\data\ICARDS\2001\may19\may19_01.009';
%  fname = fname_info_icards(fn)
%
% Author: John Paden
%
% See also datenum

[path name ext] = fileparts(fn);
fn = [name ext];

[fname.name fn] = strtok(fn,'.');

[fname.file_idx fn] = strtok(fn(2:end),'.');
fname.file_idx = str2double(fname.file_idx);

return;
