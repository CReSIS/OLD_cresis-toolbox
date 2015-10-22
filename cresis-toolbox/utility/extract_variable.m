function [varargout] = extract_variable(fns, varargin)
% [varargout] = extract_variable(fns, varargin)
%
% This function extracts variables from a given set of files
%
% Example Usage:
%  fns = get_filenames('CFF or FK Path','','','_pos.mat')
%  [Latitude, Longitude] = extract_variables(fns,'lat','lon');
%  or
%  [Latitude, Longitude, Breaks] = extract_variables(fns,'lat','lon');
%
% Author: William Blake

%% Perform Argument checks
if nargin ~= nargout 
    if nargin ~= nargout+1
        error('Arguements:Wrong','Wrong number of arguments');
    end
end
if ~iscell(fns)
    error('Format:Wrong','Argument 1 should be cell of filenames')
end
for i0 = 1:length(varargin)
    if ~ischar(varargin{i0})
        error('Format:Wrong','Arguemnt 2+ should be a variable name as a string')
    end
end
    
%% Determine Concatenation dimension
cat_dim     = zeros(1,length(varargin));
if length(fns) > 1 % Figure out what dimension is to be concatenated
    for idx0 = 1:length(varargin)
        tmp1 = whos('-file',fns{1},varargin{idx0});
        tmp2 = whos('-file',fns{2},varargin{idx0});
        % Check for existence
        if isempty(tmp1)
            error('Variable:Missing','Variable %s not present in %s', varargin{idx0}, fns{1});
        elseif isempty(tmp2)
            error('Variable:Missing','Variable %s not present in %s', varargin{idx0}, fns{1});
        end
        % Concatenate row vectors or matrix is 1st dimension
        if tmp1.size(1) == tmp2.size(1)
            cat_dim(idx0) = 1;
        end
        % Concatenate Matrix is 2nd dimension
        if (tmp1.size(2) == tmp2.size(2)) && (tmp1.size(1) ~= tmp2.size(1))
            cat_dim(idx0) = 2;
        end
        % Concatenate column vectors
        if tmp1.size(2) == 1 && tmp2.size(2) == 1
            cat_dim(idx0) = 2;
        end
        if cat_dim(idx0) == 0
            error('Inconsistent:Dimensions','Inconsistent Dimensions in %s',varargin{idx0});
        end
    end
else
    cat_dim = ones(1,length(varargin));
end

%% Determine Data Vector lengths
vars_len  = zeros(1,length(varargin));
break_idx = zeros(length(fns),length(varargin));
vars_siz_1  = zeros(1,length(varargin));
vars_siz_2  = zeros(1,length(varargin));
siz_1       = zeros(length(fns),length(varargin));
siz_2       = zeros(length(fns),length(varargin));
var_class   = cell(1,length(varargin));
var_complex = cell(1,length(varargin));

for i0 = 1:length(fns)
    for i1 = 1:length(varargin)
        tmp = whos('-file',fns{i0},varargin{i1});
        % Check for existence
        if isempty(tmp)
            error('Variable:Missing','Variable %s not present in %s', varargin{i1}, fns{i0});
        end
        % Set class and complex boolean
        var_class{i1}   = tmp.class;
        var_complex{i1} = tmp.complex;
        siz_1(i0,i1)    = tmp.size(1);
        siz_2(i0,i1)    = tmp.size(2);
        if cat_dim(i1) == 1
            % Check dimension for consistency
            if siz_1(max(1,i0-1),i1) ~= siz_1(i0,i1)
                error('Dimension:Mismatch','Dim 1 of %s does not match in file %.0f (%.0f) and %.0f (%.0f)', ...
                    varargin{i1},i0-1,siz_1(i0-1,i1),i0,siz_1(i0,i1));
            end
            vars_siz_1(i1) = siz_1(i0,i1);
            vars_siz_2(i1) = vars_siz_2(i1)+siz_2(i0,i1);
        elseif cat_dim(i1) == 2
            % Check dimension for consistency
            if siz_2(max(1,i0-1),i1) ~= siz_2(i0,i1)
                error('Dimension:Mismatch','Dim 2 of %s does not match in file %.0f (%.0f) and %.0f (%.0f)', ...
                    varargin{i1},i0-1,siz_2(i0-1,i1),i0,siz_2(i0,i1));
            end
            vars_siz_1(i1) = vars_siz_1(i1)+siz_1(i0,i1);
            vars_siz_2(i1) = siz_2(i0,i1);
        end
        break_idx(i0,i1) = break_idx(i0,i1)+tmp.size(2);
    end
end

%% Preallocate output data variables
varargout = cell(size(varargin));
for i0 = 1:length(varargin)
    if var_complex{i0} % If commplex cast as complex
        varargout{i0} = complex(ones(vars_siz_1(i0),vars_siz_2(i0),var_class{i0}));
    else
        varargout{i0} = ones(vars_siz_1(i0),vars_siz_2(i0),var_class{i0});
    end
end

%% Load in data
idx_start = ones(1,length(varargin));
for i0 = 1:length(fns)
    for i1 = 1:length(varargin)
        tmp_info = whos('-file',fns{i0},varargin{i1});
        tmp      = load(fns{i0},varargin{i1});
        if cat_dim(i1) == 1;
            idxs                  = idx_start(i1):idx_start(i1)+tmp_info.size(2)-1;
            s_data                = ['tmp.',varargin{i1}];
            varargout{i1}(:,idxs) = eval(s_data);
            idx_start(i1)         = idx_start(i1)+tmp_info.size(2);
        elseif cat_dim(i1) == 2;
            idxs                  = idx_start(i1):idx_start(i1)+tmp_info.size(1)-1;
            s_data                = ['tmp.',varargin{i1}];
            varargout{i1}(idxs,:) = eval(s_data);
            idx_start(i1)         = idx_start(i1)+tmp_info.size(1);
        end
    end
end

if nargin ~= nargout+1
    varargout{length(varargin)+1} = cumsum(break_idx(:,1));
end

return;
