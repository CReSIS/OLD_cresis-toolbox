function fname = fname_info_arena(fn)
% fname = fname_info_arena(fn)
%
% Detects the type of Arena filename and then parses it.
%
% fn = Arena file name
%
% fname = structure containing each of the fields of the MCORDS
%   filename
%
%   For example: 20161108_090811_system.xml
%  .type = 'system'
%  .datenum = time stamp converted to Matlab date number. Common
%     usage: "[year month day hour min sec] = datevec(fname.datenum)"
%
%   For example: 20161108_092527__1804.dat
%  .file_idx = 1804
%  .datenum = time stamp converted to Matlab date number. Common
%     usage: "[year month day hour min sec] = datevec(fname.datenum)"
%
% Example
%  fn = '/cresis/snfs1/data/HF_Sounder/2016_Greenland_TO/20161108A/20161108_092527_system.xml';
%  fname = fname_info_arena(fn)
%
%  fn = '/cresis/snfs1/data/HF_Sounder/2016_Greenland_TO/20161108A/20161108_092527__1791.dat';
%  fname = fname_info_arena(fn)
%
% Author: John Paden
%
% See also datenum

[~,fn,fn_ext] = fileparts(fn);

fname = [];
if strcmpi(fn_ext,'.xml')
  [day_str fn] = strtok(fn,'_');
  [time_str fn] = strtok(fn,'_');
  fname.datenum = datenum([day_str time_str],'yyyymmddHHMMSS');
  
  fname.type = strtok(fn(2:end),'.');
  
elseif strcmpi(fn_ext,'.dat')
  [day_str fn] = strtok(fn,'_');
  [time_str fn] = strtok(fn,'_');
  fname.datenum = datenum([day_str time_str],'yyyymmddHHMMSS');

  if fn(2) == '_'
    fn = fn(2:end);
  else
    [~,fn] = strtok(fn(2:end),'_');
  end

  [fname.file_idx fn] = strtok(fn(2:end),'.');
  fname.file_idx = str2double(fname.file_idx);
end

return;
