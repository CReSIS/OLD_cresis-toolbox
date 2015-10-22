function [params] = read_param_xls(param_fn, day_seg_filter, generic_ws)
% [params] = read_param_xls(param_fn, day_seg_filter, generic_ws)
%
% Reads params structure array in from .xls file.
%
% This params structure array is used by master/master_accum to control the
% data process.
%
% param_fn = filename to MS Excel params spreadsheet
% day_seg_filter = regular expression string, only the matching day-segment
%   parameters will be returned
% generic_ws = cell array of strings, loads generic worksheets specified in
%   the cell array in addition to the standard worksheets (e.g. analysis,
%   update_data_files). Legacy support allows generic_ws to be a string
%   and this function will detect this and place it into a cell array.
%
% Example:
%   param_fn = '/mnt/scratch2/csarp_support/documents/mcords_param_2011_Greenland_P3.xls';
%   [params] = read_param_xls(param_fn);
%
%   param_fn = '/mnt/scratch2/csarp_support/documents/mcords_param_2011_Greenland_P3.xls';
%   [params] = read_param_xls(param_fn,'20110329_02');
%
%   param_fn = '/mnt/scratch2/csarp_support/documents/accum_param_2011_Greenland_P3.xls';
%   [params] = read_param_xls(param_fn);
%
% Author: Brady Maasen, John Paden
%
% See also: master

cell_boolean = @read_param_xls_boolean;
cell_text = @read_param_xls_text;
cell_read = @read_param_xls_general;

% ======================================================================%
% CREATING THE PARAM STRUCTURE ARRAY FROM PARAM_STARTER.XLS
% ======================================================================%
warning('off','MATLAB:xlsread:Mode');

% =======================================================================
% Create Command Parameters
% =======================================================================
sheet_name = 'command';
fprintf('Reading sheet %s of xls file: %s\n', sheet_name, param_fn);

%% Load standard worksheets
[params] = read_param_xls_radar(param_fn);

%% Load the generic worksheets if specified
if exist('generic_ws','var')
  if ischar(generic_ws)
    % Legacy support to allow a string to be passed into generic_ws variable.
    generic_ws = {generic_ws};
  end
  for idx = 1:size(generic_ws,1)
    tmp = read_param_xls_generic(param_fn,generic_ws{idx,1},params);
    if size(generic_ws,2) > 1
      % Rename the worksheet variable
      [params.(generic_ws{idx,2})] = tmp.(generic_ws{idx,1});
    else
      params = tmp;
    end
  end
end

%% Just get the specific day_seg that was requested
if exist('day_seg_filter','var') && ~isempty(day_seg_filter)
  good_mask = logical(zeros(size(params)));
  for idx = 1:length(params)
    if ~isempty(regexp(params(idx).day_seg, day_seg_filter))
      good_mask(idx) = 1;
    end
  end
  params = params(good_mask);
  if isempty(params)
    warning('No segment day_seg matched regexp %s .', day_seg_filter);
  end
end

return

