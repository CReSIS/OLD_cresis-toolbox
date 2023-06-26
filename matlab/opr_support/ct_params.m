classdef ct_params < handle
% ct_params < handle
% 
% A class that manages cresis toolbox parameter Excel spreadsheets.
% NOT FINISHED
% 
% Author: John Paden
%
% See also: 
  
  properties (Constant)
    current_version = 5.0;
  end
  
  properties
    params
  end
  
  methods    
    %% constructor   
    function obj = ct_params(fn, day_seg_filter, generic_ws)
      
      if nargin == 0
        obj.params = [];
        
      elseif ischar(fn)
        %% Load standard worksheets
        warning('off','MATLAB:xlsread:Mode');
        [params] = read_param_xls_radar(param_fn);
        
        if isempty(params) || isempty(params(1).day_seg)
          warning('Parameter spreadsheet file is empty');
        end
        
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
        warning('on','MATLAB:xlsread:Mode');
        
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
        
        
        obj.read(fn, day_seg_filter, generic_ws);
        
      elseif isstruct(fn)
        
      else
        error('Input is not a char or struct.');
        
      end
    end
  
    function read(obj,fn)
      obj.read_cmd(fn);
      obj.read_records(fn);
    end
  
    function read_cmd(obj,fn)
    end
  
    function read_records(obj,fn)
    end
  
    function print(obj,fn)
    end
    
    function print_cmd(obj,fn)
    end
    
    function print_records(obj,fn)
    end
  
    function write(obj,fn)
    end
    
    function write_cmd(obj,fn)
    end
    
    function write_records(obj,fn)
    end
  end
  
  methods(Static)
    
  end
  
end

