function [prj_struct data] = read_prj(fn,data)
% [prj_struct] = read_prj(fn)
%
% Reads in a .prj projection file (well known text or WKT format).
% Do not pass more than one argument to the function or it will fail.
%
% Input:
% fn = string with filename to load
%
% Output:
% prj_struct = structure containing the prj contents in a struct
%   .fields: cell array of all the fields associated with this structure
%   .s_STRUCT_NAME: sub-structure (which will contain its own fields
%      and sub-structures)
% 
% Example:
% fn = '/cresis/data3/GIS_data/greenland/kpurdon/MOA_BoundaryLines/ORIGINAL_Files/moa_coastline.prj';
% prj = read_prj(fn);
%
% Author: John Paden

if ~exist('data','var')
  % The first time the function is called we read in the whole file
  fid = fopen(fn,'r');
  data = fread(fid,inf,'uint8=>char').';
  fclose(fid);
end

% Recursively read in the fields
prj_struct.fields = {};
while ~isempty(data)
  [data,field,type] = read_prj_field(data);
  if type == 0 || type == 1
    % Just a normal field
    prj_struct.fields{end+1} = field;
  elseif type == 2
    % Recurse since this is a structure
    [prj_struct.(['s_' field]) data] = read_prj(fn,data);
  end
  if length(data) >= 1 && data(1) == ']'
    data = data(2:end);
    return;
  end
end

return;

function [data,field,type] = read_prj_field(data)

data_idx = 1;
while data_idx < length(data)
  if data(data_idx) == '"'
    % Beginning of string (read until '"')
    type = 0;
    data_idx = data_idx + 1;
    new_data_idx = data_idx-1 + find(data(data_idx:end) == '"',1);
    if isempty(new_data_idx)
      error('string with out ending "');
    end
    field = data(data_idx:new_data_idx-1);
    data = data(new_data_idx+1:end);
    return;
  elseif isstrprop(data(data_idx),'digit')
    % Beginning of number (read until ',', or ']')
    type = 1;
    new_data_idx = data_idx-1 + find(data(data_idx:end)==',' | data(data_idx:end)==']',1);
    if isempty(new_data_idx)
      error('numeric with out ending ,]');
    end
    field = str2double(data(data_idx:new_data_idx-1));
    data = data(new_data_idx:end);
    return;
  elseif isstrprop(data(data_idx),'alpha')
    % Beginning of struct (read until '[')
    type = 2;
    new_data_idx = data_idx-1 + find(data(data_idx:end)=='[',1);
    if isempty(new_data_idx)
      error('struct with out [');
    end
    field = data(data_idx:new_data_idx-1);
    data = data(new_data_idx+1:end);
    return;
  end
  data_idx = data_idx + 1;
end

% End of file
field = '';
type = 3;
data = [];

return;


