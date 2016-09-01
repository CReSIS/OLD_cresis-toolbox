function fname = fname_info_acords(fn,param)
% fname = fname_info_acords(fn)
%
% Detects the type of acords filename and then parses it.
%
% fn = acords file name
% param.hnum = the index of the header within a raw file (since ACORDS supports
% mutiple different headers in a file)
% param.file_version
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

if ~isfield(param,'hnum')
  param.hnum = 1;
end
if ~isfield(param,'file_version')
  param.file_version = 406;
end

[hdr] = basic_load_acords(fn,struct('datatype',0,'file_version',param.file_version));

[path name ext] = fileparts(fn);
fn = [name ext];

fname.basename = name;
fname.file_idx= str2double(ext(2:end));

if strcmp(fn(1:3),'apr')
  month = '04';
  day = fn(4:5);
elseif strcmp(fn(1:3),'may')
  month = '05';
  day = fn(4:5);
elseif strcmp(fn(1:3),'nov')
  month = '11';
  day = fn(4:5);
  year = '2004';
end

[fname.name fn] = strtok(fn,'_');
if ~exist('year','var')
  year = ['20', fn(2:3)];
end


day_str = [year month day];
[Y,M,D,h,m,s] = datevec(epoch_to_datenum(hdr(param.hnum).time(1)));
time_str = [num2str(h) num2str(m) num2str(round(s))];
fname.datenum = datenum([day_str time_str],'yyyymmddHHMMSS');



return;
