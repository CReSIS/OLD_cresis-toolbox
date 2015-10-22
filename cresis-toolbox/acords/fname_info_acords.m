function fname = fname_info_acords(fn)
% fname = fname_info_acords(fn)
%
% Detects the type of acords filename and then parses it.
%
% fn = acords file name
%
% fname = structure containing each of the fields of the accum
%   filename
%
%   For example:
%     basename.file_idx (may09_03.6 or may09_03b.13)
%  .basename = can be anything without dots
%  .file_idx = 0 to 99
%
% Example
%  fn = '/cresis/snfs1/data/ACORDS/airborne2003/may09/may09_03b.35';
%  fname = fname_info_acords(fn)
%
% Author: John Paden
%
% See also datenum

[path name ext] = fileparts(fn);
fn = [name ext];

fname.basename = name;
fname.file_idx= str2double(ext(2:end));

return;
