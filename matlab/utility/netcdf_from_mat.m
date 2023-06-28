function netcdf_from_mat(cdf_fn,mat,varmap)
% netcdf_from_mat(cdf_fn,mat,varmap)
%
% Writes Matlab structure OR Matlab file to NetCDF4 file.
%
% cdf_fn = output netcdf filename
% mat = struct OR character string of matlab filename
% varmap = structure array where each element of the structure defines
%  a variable mapping from the mat file to the new cdf file. If a variable
%  is not specified in the varmap structure, then a default mapping occurs.
%  In other words, varmap is optional. However, if you do specify a mapping
%  for a variable, you need to specify at least these three fields
%   varmap(n).mat_name = 'GPS_time';
%   varmap(n).cdf_name = 'time';
%   varmap(n).dim_name = {'time'};
%  The following fields are optional. "eval" causes the passed in function
%  to be evaluated on the variable in question before it is written to the
%  netcdf file. This is useful for converting units.
%   varmap.eval = @(x) (x - utc_leap_second(x(1)));
%  The "attributes" field is a 2N length cell array where each pair
%  of values provides an attribute name and an attribute value.
%   varmap.attributes = {'units' 'seconds since 1970-01-01 00:00:00', ...
%     'calendar', 'standard', ...
%     'long_name', 'Time of day UTC', ...
%     'standard_name', 'time', ...
%     'axis', 'T'};
%
% For each field, collect all field names and dimensions of each field
% Create the minimum number of dimensions (removes duplicate dimensions
%   since it is recommended to keep number of dimensions to less than 100)
% Create the variables referencing these dimensions
%   Add matlab_class attribute for reading the variables out
%   Add matlab_size attribute for reading the variables out
% Put the data in the variables
%
% Not currently preserved:
% * P-element cell arrays of 1xN character strings are converted to MxP character
%   strings where M is the longest string and strings are filled with zero.
%   There are given a new class "cell_string".
% * struct and cell arrays are flattened into 1xN at each level in the
%   hierarchy
% * All 1xN strings are converted to fixed length with the longest string
%   in the file setting this length, all strings are then null/zero filled
%
% Author: John Paden
%
% See also: run_netcdf_from_mat.m, netcdf_type.m, netcdf_to_mat.m

if ~exist('varmap','var')
  varmap = [];
end

%% Load variables
if ischar(mat)
  % Passed in a filename
  mat = load(mat);
end
 
%% Get list of variables from the MAT file
vars = netcdf_from_mat_get_vars(mat,'');

%% Get a unique list of all the mapped variable's dimensions (mapped_var_dim_names)
mapped_var_dim_names = {};
for map_idx = 1:length(varmap)
  mapped_var_dim_names = cat(2,mapped_var_dim_names, ...
    reshape(varmap(map_idx).dim_name,[1 numel(varmap(map_idx).dim_name)]));
end
if ~isempty(mapped_var_dim_names)
  % If empty, unique returns erroneous 0-by-1 cell array
  mapped_var_dim_names = unique(mapped_var_dim_names);
end

%% Find the size of each dimension (mapped_var_dim_sizes)
mapped_var_dim_sizes = -1 * ones(size(mapped_var_dim_names));
if ~isempty(vars)
  [vars.map_idx] = deal(zeros(size(vars)));
end
for map_idx = 1:length(varmap)
  var_idx = strmatch(varmap(map_idx).mat_name,{vars.name},'exact');
  
  if isempty(var_idx)
    %  Variable does not exist
    continue;
  end
  
  % Record the mapping for later
  vars(var_idx).map_idx = map_idx;
  
  % Check to make sure the number of dimensions match
  num_dims_spec = length(varmap(map_idx).dim_name);
  num_dims = sum(vars(var_idx).size ~= 1);
  if num_dims == 0
    num_dims = 1;
  end
  if num_dims ~= num_dims_spec
    error('Number of dimensions for %s (%d) does not match what is in mat file (%d)', ...
      varmap(map_idx).mat_name, num_dims_spec, num_dims);
  end
  
  % Go through each dimension and make sure the size of each dimension
  % is consistent
  for dim_idx = 1:length(varmap(map_idx).dim_name)
    dim = strmatch(varmap(map_idx).dim_name{dim_idx}, mapped_var_dim_names, 'exact');
    if length(varmap(map_idx).dim_name) == 1
      dim_size = find(vars(var_idx).size ~= 1,1);
      if isempty(dim_size)
        dim_size = prod(vars(var_idx).size);
      else
        dim_size = vars(var_idx).size(dim_size);
      end
    else
      dim_size = vars(var_idx).size(dim_idx);
    end
    if mapped_var_dim_sizes(dim) == -1
      mapped_var_dim_sizes(dim) = dim_size;
    elseif mapped_var_dim_sizes(dim) ~= dim_size
      error('Dimension of %s (%d) does not match expected size (%d)', ...
        varmap(map_idx).mat_name, dim_size, mapped_var_dim_sizes(dim));
    end
  end
end

%% Remove mapped var dimensions that are never used
good_mapped_var_dims = mapped_var_dim_sizes ~= -1;
mapped_var_dim_sizes = mapped_var_dim_sizes(good_mapped_var_dims);
mapped_var_dim_names = mapped_var_dim_names(good_mapped_var_dims);

%% Sort through dimensions, create a list of unique dim sizes (dims)
dims = [];
max_char_length = NaN;
for var_idx=1:length(vars)
  if vars(var_idx).map_idx == 0
    % Only consider dimensions of variables that are not in map_var list
    if strcmp(vars(var_idx).class,'char') && length(vars(var_idx).size) == 2
      if isnan(max_char_length) || max_char_length < prod(vars(var_idx).size)
        max_char_length = prod(vars(var_idx).size);
      end
    else
      dims = cat(2,dims,vars(var_idx).size);
    end
  end
end
if ~isnan(max_char_length)
  max_char_length = max_char_length + 1; % Allow for zero termination
  dims = [dims max_char_length];
end
[dims dim_idxs] = unique(dims);

%% Create NetCDF file
cdf_fn_dir = fileparts(cdf_fn);
if ~exist(cdf_fn_dir,'dir')
  mkdir(cdf_fn_dir);
end

try
  if exist(cdf_fn,'file')
    delete(cdf_fn); % Don't try to merge with existing file, start over
  end
  ncid = netcdf.create(cdf_fn,bitor(netcdf.getConstant('CLOBBER'),netcdf.getConstant('NETCDF4')));
catch ME
  warning('Exception during file creation');
  ME
  keyboard
end
netcdf.setFill(ncid,'NC_NOFILL');

%% Write dimension information to netcdf file (dim_ids, map_var.dim_ids)
% -------------------------------------------------------------------------
% Automatically generated dimensions
dim_ids = [];
for dim_idx = 1:length(dims)
  dim_name = sprintf('auto_%d',dims(dim_idx));
  dim_ids(dim_idx) = netcdf.defDim(ncid, dim_name, dims(dim_idx));
end
% Mapped variable dimensions
mapped_var_dim_ids = [];
for dim_idx = 1:length(mapped_var_dim_sizes)
  mapped_var_dim_ids(dim_idx) = netcdf.defDim(ncid, mapped_var_dim_names{dim_idx}, ...
    mapped_var_dim_sizes(dim_idx));
end

%% Write variable definitions to netcdf file 
% -------------------------------------------------------------------------
for var_idx = 1:length(vars)
  
  %% Find dimension ids and variable name
  vars(var_idx).dim_ids = [];
  if vars(var_idx).map_idx == 0
    % Unmapped variable dimensions
    var_name = vars(var_idx).name;
    if strcmp(vars(var_idx).class,'char') && length(vars(var_idx).size) == 2
      % Special case for strings
      vars(var_idx).dim_ids(1) = dim_ids(find(dims == 1));
      vars(var_idx).dim_ids(2) = dim_ids(find(dims == max_char_length));
    else
      % All other variables besides strings
      for dim_idx = 1:length(vars(var_idx).size)
        vars(var_idx).dim_ids(dim_idx) = dim_ids(find(dims == vars(var_idx).size(dim_idx)));
      end
    end
    
  else
    % Mapped variable dimensions
    map_idx = vars(var_idx).map_idx; % Short variable name
    var_name = varmap(map_idx).cdf_name;
    for dim_idx = 1:length(varmap(map_idx).dim_name)
      vars(var_idx).dim_ids(dim_idx) ...
        = mapped_var_dim_ids(strmatch(varmap(map_idx).dim_name(dim_idx), mapped_var_dim_names, 'exact'));
    end
  end
  
  %% Create variables
  vars(var_idx).id = netcdf.defVar(ncid,var_name, ...
    netcdf_type(vars(var_idx).class),vars(var_idx).dim_ids);
  netcdf.defVarFill(ncid,vars(var_idx).id,true,0);
  
  %% Construct special variable attributes
  if vars(var_idx).map_idx ~= 0
    map_idx = vars(var_idx).map_idx;
    if isfield(varmap(map_idx),'attributes')
      for attrib_idx = 1:2:length(varmap(map_idx).attributes)
        netcdf.putAtt(ncid,vars(var_idx).id, ...
          varmap(map_idx).attributes{attrib_idx}, ...
          varmap(map_idx).attributes{attrib_idx+1});
      end
    end
  end
  
  %% Construct general variable attributes
  netcdf.putAtt(ncid,vars(var_idx).id, 'matlab_class', vars(var_idx).class);
  netcdf.putAtt(ncid,vars(var_idx).id, 'matlab_size', vars(var_idx).size);
end

%% Write variable data to the file
netcdf.endDef(ncid);
for var_idx = 1:length(vars)
  
  % For mapped variables with eval field assigned, run the function
  % specified in eval on the data variable.
  if vars(var_idx).map_idx ~= 0
    map_idx = vars(var_idx).map_idx;
    if isfield(varmap(map_idx),'eval') ...
      && ~isempty(varmap(map_idx).eval)
      mat.(vars(var_idx).name) = varmap(map_idx).eval(mat.(vars(var_idx).name));
    end
  end
  
  if strcmp(vars(var_idx).class,'cell_string')
    %% Repackage cell array of strings into 2-D string array
    M = char(zeros(vars(var_idx).size,'uint8'));
    eval(['cell_array = mat.' vars(var_idx).name ';']);
    for cell_idx = 1:numel(cell_array)
      string_length = length(cell_array{cell_idx});
      M(end-2) = char(0); % null terminate
      M(end-1) = char(string_length/2^8); % encode original length of string
      M(end) = char(mod(string_length,2^8));
      M(1:length(cell_array{cell_idx}),cell_idx) = cell_array{cell_idx};
    end
    netcdf.putVar(ncid,vars(var_idx).id,M);
  elseif strcmp(vars(var_idx).class,'char')
    %% Repackage character string into fixed length string
    eval(['M = mat.' vars(var_idx).name '(:);']);
    if numel(M) < max_char_length
      M(max_char_length) = 0;
    end
    netcdf.putVar(ncid,vars(var_idx).id,M);
  elseif strcmp(vars(var_idx).class,'logical')
    netcdf.putVar(ncid,vars(var_idx).id,uint8(eval(['mat.' vars(var_idx).name ';'])));
  elseif any(strcmp(vars(var_idx).class,{'function_handle','inline'}))
    netcdf.putVar(ncid,vars(var_idx).id,char(eval(['mat.' vars(var_idx).name ';'])));
  else
    eval(['M = mat.' vars(var_idx).name ';']);
    % NetCDF does not allow empty variables to be filled with []
    % and we don't need to fill them anyway...
    if ~isempty(M)
      netcdf.putVar(ncid,vars(var_idx).id,M);
    end
  end
end
%netcdf.reDef(ncid);

%% Clean up and exit
netcdf.close(ncid);

return

function vars = netcdf_from_mat_get_vars(mat,parent)
%% vars = netcdf_from_mat_get_vars(mat,parent)
%
% Get variable information from Matlab structure
%
% Recursive function:
% mat: structure or cell array that we are pulling variables out of
% parent: parent variable name path string, e.g. 'field1.field2(1).' when
%   field2 is a struct array and 'field1.field2' for when field2 is a cell
%   array.

debug_level = 0;
if debug_level > 0
  indention = ' ' * ones(1,numel(find(parent=='.')) * 2);
end

vars = [];
if isstruct(mat)
  mat_fieldnames = fieldnames(mat);
elseif iscell(mat)
  mat_fieldnames = 1:length(mat);
end
for field_idx = 1:length(mat_fieldnames)
  if iscell(mat_fieldnames)
    field_name = mat_fieldnames{field_idx};
    field = mat.(field_name);
  else
    field_name = sprintf('{%d}',field_idx);
    field = mat{field_idx};
  end
  
  if isstruct(field) && ~isempty(field)
    %% Structure Array: recurse on each element of the structure array
    if debug_level > 0
      fprintf('%s%s\n', indention, field_name)
    end
    for idx = 1:numel(field)
      new_vars = netcdf_from_mat_get_vars(field(idx),[parent field_name sprintf('(%d)',idx) '.']);
      vars = cat(2,vars,new_vars);
    end
  elseif iscell(field) && ~isempty(field)
    %% Check to see if this is the special cell_string case (i.e.
    % a cell array of character strings)
    not_char = false;
    for idx = 1:numel(field)
      if ~ischar(field{idx})
        not_char = true;
        break;
      end
    end
    %% Get variable information
    if not_char
      %% Generic cell array, we treat it like a single element struct
      if debug_level > 0
        fprintf('%s%s\n', indention, field_name)
      end
      new_vars = netcdf_from_mat_get_vars(field,[parent field_name]);
      vars = cat(2,vars,new_vars);
    else
      %% Special cell_string case
      % We are going to repackage into 2-D char array, so we need to
      % get the maximum string length.
      max_char_length = 0;
      for idx = 1:numel(field)
        if length(field{idx}) > max_char_length
          max_char_length = length(field{idx});
        end
      end
      vars(end+1).name = [parent field_name];
      vars(end).size = [max_char_length+3 numel(field)];
      vars(end).class = 'cell_string';
      if debug_level > 0
        fprintf('%s%s', indention, vars(end).name)
        fprintf(' %s', vars(end).class)
        fprintf(' %d', vars(end).size)
        fprintf('\n');
      end
    end
  else
    %% Basic field types
    vars(end+1).name = [parent field_name];
    vars(end).class = class(field);
    if isa(field,'inline') || isa(field,'function_handle')
      field = char(field);
    end
    vars(end).size = size(field);
    if debug_level > 0
      fprintf('%s%s', indention, vars(end).name)
      fprintf(' %s', vars(end).class)
      fprintf(' %d', vars(end).size)
      fprintf('\n');
    end
  end
end

return;
