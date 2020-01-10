classdef surfdata < handle
% tomo.surfdata(fn) < handle
% 
% A class that manages surf data in a frame. 
% 
% fn: optional string containing a surfdata filename to load, default
%   is to create an empty class
% 
% Example of a filepath:
% X:\ct_data\rds\2014_Greenland_P3\CSARP_surfData\20140401_03\Data_20140401_03_037.mat
% 
% Example of how to generate a filepath:  
%   param.radar_name   = 'rds';
%   param.season_name  = '2014_Greenland_P3';
%   param.day_seg = '20140401_03'
%   surfdata_ref = 'surfData';
%   frm = 37;
%   fn = fullfile(ct_filename_out(param,surfdata_ref,''),...
%     sprintf('Data_%s_%03d.mat',param.day_seg,frm));      
% 
% Author: Shane Chu, John Paden
%
% See also: run_surfdata.m, surfdata.m
  
  properties (Constant)
    current_version = 2.0;
  end
  
  properties
    % surf: A structure array; each element holds a single surface
    %  .name: string containing the name of the surface (e.g. 'top',
    %    'bottom')
    %  .x: Nsv by Nx array containing the x-coordinates of the surface points
    %  .y: Nsv by Nx array containing the y-coordinates of the surface points
    %  .visible: Logical scalar indicating whether or not the surface should
    %    be plotted.
    %  .top: Index into surf structure which indicates the "top"
    %    surface to use when this surface is selected.
    %  .active: Index into surf structure which indicates the "active"
    %    surface to use when this surface is selected.
    %  .mask: Index into surf structure which indicates the "ice mask"
    %    surface to use when this surface is selected.
    %  .control: Index into surf structure which indicates the "control"
    %    surface to use when this surface is selected.
    %  .quality: Index into surf structure which indicates the "quality"
    %    surface to use when this surface is selected.
    %  .plot_name_values: Plot (name,value) pairs stored in a cell array
    %    which will be passed in when ever this surface is plotted.
    surf
    
    % radar_name: string containing the radar name
    radar_name
    % season_name: string containing the season name
    season_name
    % day_seg: string containing the segment
    day_seg
    % frm: integer scalar containing the frame number
    frm

    % gps_time: GPS time of each column in surf.[xy]
    gps_time
    
    % theta: Nsv by 1 matrix containing DOA for each row in surf.[xy]
    % (zero theta is -z, theta increases toward positive y)
    theta

    % time: fast (range) time which allows surf.y values (which are in
    %   range bins) to be converted to two way travel time
    time
    
    % FCS: flight (SAR) coordinate system for each column in surf.[xy]
    %  .origin: origin, 3 by Nx
    %  .x: unit x-axis vector, 3 by Nx (points along-track)
    %  .y: unit y-axis vector, 3 by Nx (points left, but equal to cross(x,z))
    %  .z: unit z-axis vector, 3 by Nx (points up, but orthogonal to x)
    FCS
  end
  
  methods    
    %% constructor   
    function obj = surfdata(fn,param,frm,gps_time,theta,FCS)
      % no arguments => create an empty surf struct array
      if nargin == 0
        obj.surf = [];
        obj.radar_name = [];
        obj.season_name = [];
        obj.day_seg = [];
        obj.frm = [];
        obj.gps_time = [];
        obj.theta = [];
        obj.time = [];
        obj.FCS = [];
        
        % one argument => load file and create struct array accordingly
      else
        % Check to see if this is a direction of arrival method or a beam
        % forming method file.
        if ~exist('param','var') || ~isfield(param,'doa_method_flag')
          doa_method_flag = false;
        else
          doa_method_flag = param.doa_method_flag;
        end
        try
          obj.surf = [];
          obj.load_surfdata(fn,doa_method_flag);
          
        catch ME
          error('Failed to load filename (%s):\n %s', fn, ME.getReport());
        end
      end
    end
    
    %% load_surf
    function load_surfdata(obj, fn,doa_method_flag)
      % obj.load_surfdata(fn)
      %
      % Loads a new surf struct into the object. If fn is a filename,
      % it should be a surfData file. Note that this function will overwrite
      % the existing surf struct array.
      %
      % Input:
      %  fn: The file path OR a new struct array with all the fields in a
      %   surfData file.
      %  doa_method_flag: logical that specifies is this file stores DOA
      %   data or beam forming data. Default is false.
      
      if ~exist('doa_method_flag','var')
        doa_method_flag = false;
      end
      
      if ischar(fn)
        if obj.version(fn) < obj.current_version
          error('File "%s" is old. Run tomo.surfdata.update on the file.', fn);
        end
      end
      
      old_surf = obj.surf;
      try
        if ischar(fn)
          tmp = load(fn);
        else
          tmp = fn;
        end
        obj.surf = tmp.surf;
        
        for surf_idx = 1:length(obj.surf)
          obj.valid_surf(obj.surf(surf_idx),doa_method_flag);
        end
        
        obj.set_metadata(tmp,doa_method_flag);
        
      catch ME
        obj.surf = old_surf;
        if ischar(fn)
          error('Failed to load %s:\n %s', fn, ME.getReport());
        else
          error('Failed to load struct:\n %s', ME.getReport());
        end
      end
    end
    
    function set_metadata(obj, md, doa_method_flag)
      % obj.set_metadata(md)
      % 
      % Function for validating all the metadata fields
      %
      % md: struct containing gps_time, theta, FCS, radar_name,
      %   season_name, day_seg, and frm.
      % doa_method_flag: logical that specifies is this file stores DOA
      %  data or beam forming data. Default is false.
      
      obj.valid_metadata(md,doa_method_flag);
      
      obj.radar_name = md.radar_name;
      obj.season_name = md.season_name;
      obj.day_seg = md.day_seg;
      obj.frm = md.frm;
      obj.gps_time = md.gps_time;
      obj.theta = md.theta;
      obj.time = md.time;
      obj.FCS = md.FCS;
    end
    
    function valid_metadata(obj, md,doa_method_flag)
      % obj.valid_metadata(md)
      % 
      % Function for validating all the metadata fields
      %
      % md: struct containing gps_time, theta, FCS, radar_name,
      %   season_name, day_seg, and frm.
      % doa_method_flag: logical that specifies is this file stores DOA
      %  data or beam forming data. Default is false.
      
      if ~ischar(md.radar_name)
        error('radar_name must be a string');
      end
      
      if ~ischar(md.season_name)
        error('season_name must be a string');
      end
      
      if ~ischar(md.day_seg)
        error('day_seg must be a string');
      end
      
      if ~isnumeric(md.frm)
        error('frm must be a positive integer');
      end
      
      Nx = size(md.gps_time,2);

      if size(md.gps_time,1) ~= 1
        error('gps_time must be a row vector.');
      end
      
      if ~doa_method_flag
        if ~isempty(obj.surf) && size(obj.surf(1).x,2) ~= Nx
          error('gps_time must have same number of columns as surf.x');
        end
        
        Nsv = size(md.theta,1);
        
        if size(md.theta,2) ~= 1
          error('theta must be a column vector.');
        end
        
        if ~isempty(obj.surf) && size(obj.surf(1).x,1) ~= Nsv
          error('theta must have same number of rows as surf.x');
        end
      end
      
      if size(md.time,2) ~= 1
        error('time must be a column vector of numbers');
      end
      
      if ~isnumeric(md.time)
        error('time must be a column vector of numbers');
      end
      
      if size(md.FCS.origin,2) ~= Nx
        error('FCS.origin must have same number of columns as gps_time.');
      end
      
      if size(md.FCS.x,2) ~= Nx
        error('FCS.x must have same number of columns as gps_time.');
      end
      
      if size(md.FCS.y,2) ~= Nx
        error('FCS.y must have same number of columns as gps_time.');
      end
      
      if size(md.FCS.z,2) ~= Nx
        error('FCS.z must have same number of columns as gps_time.');
      end
      
      if size(md.FCS.origin,1) ~= 3
        error('FCS.origin must have 3 rows.');
      end
      
      if size(md.FCS.x,1) ~= 3
        error('FCS.x must have 3 rows.');
      end
      
      if size(md.FCS.y,1) ~= 3
        error('FCS.y must have 3 rows.');
      end
      
      if size(md.FCS.z,1) ~= 3
        error('FCS.z must have 3 rows.');
      end
      
    end

    %% insert_surf
    function insert_surf(obj, surf_struct,doa_method_flag)
      % obj.insert_surf(surf_struct)
      %
      % Input: 
      %   surf_struct: A surf structure. Normally, the indexing fields
      %   such as top, active, mask, gt, quality should be cleared.
      %  doa_method_flag: logical that specifies is this file stores DOA
      %   data or beam forming data. Default is false.
      %
      % Result: 
      %   Adds a surf structure, surf_A, into the "surf" array of the object
      %
      % See also: surfdata.clear_references
      
      if ~exist('doa_method_flag','var')
        doa_method_flag = false;
      end
      % check if it is a valid surf structure
      obj.valid_surf(surf_struct,doa_method_flag);
      
      % check if name of the surface already existed in the class
      if isempty(obj.surf)
        obj.surf = [obj.surf surf_struct];
      elseif any(strcmpi(surf_struct.name,{obj.surf.name}))
        error('This surface is already inserted.');  
      else
        obj.surf = [obj.surf surf_struct];
      end
      
    end
    %% get_surf
    function surf = get_surf(obj, surf_name)
      % surf = obj.get_surf(surf_name)
      %
      % Input: 
      %   surf_name: Must be one of the following:
      %     1. a string containing a valid surface name
      %     2. a cell array of strings containing valid surface names
      %     3. a numeric array of surface indices 
      %
      % Return: 
      %   A surf structure array matching surf_name
      
      if isa(surf_name, 'char')
        index = find(strcmpi(surf_name, {obj.surf.name}));
        if index
          surf = obj.surf(index);
        else
          error('Surface "%s" does not exist.', surf_name);
        end
        
      elseif isa(surf_name, 'cell')
        surf = [];
        match_idxs = zeros(size(surf_name));
        for idx = 1:length(surf_name)
          new_match_idx = find(strcmpi(surf_name{idx}, {obj.surf.name}));
          if isempty(new_match_idx)
            error('Surface "%s" does not exist.', surf_name{idx});
          end
          match_idxs(idx) = new_match_idx;
        end
        surf = obj.surf(match_idxs);
        
      elseif isnumeric(surf_name)
        try
          surf = obj.surf(surf_name);
        catch ME
          error('Invalid indices provided:\n\n%s', ME.getReport);
        end
        
      else
        error('Input must be either a string containing a valid surface name, a cell array of strings containing valid surface names, or a numeric array of surface indices.');
      end    
    end
    
    %% set
    function set(obj, surf_name, varargin)
      % obj.set(surf_name, NAME_VALUE_PAIRS)
      %
      % Input: 
      % surf_name: surface name or cell array of surface names
      % NAME_VALUE_PAIRS: Pairs of input arguments where the first
      %   argument in the pair is the property to change and the second
      %   argument in the pair is the value.
      %
      % Return: 
      %   1. Nothing on success
      %   2. Throws an error if surface_old_name is not found.

      surf_idxs = obj.get_index(surf_name,true);
      
      for arg_idx = 3:2:nargin
        field_name = lower(varargin{arg_idx-2});
        switch field_name
          case 'name'
            if length(surf_idxs) ~= 1
              error('The name property can only be set for one surface at a time.');
            end
            new_name = varargin{arg_idx-1};
            if any(strcmpi(new_name, {obj.surf([1:surf_idxs-1, surf_idxs+1:end]).name}))
              error('The new name matches an existing name for another surface.');
            end
            obj.surf(surf_idxs).name = new_name;

          case {'top','active','mask','gt','quality'}
            new_surf_idx = obj.get_index(varargin{arg_idx-1});
            if length(new_surf_idx) ~= 1
              error('Field %s must be assigned to match a single existing surface.', varargin{arg_idx});
            end
            for surf_idx = surf_idxs
              obj.surf(surf_idx).(field_name) = new_surf_idx;
            end
          
          otherwise
            error('Invalid field %s', varargin{arg_idx-2});
            
        end
        
      end
      
    end
    
    %% set_surf
    function obj = set_surf(obj, surf_struct)
      % obj.set_surf(surf_struct)
      %
      % Input: 
      %   surf_struct: A surf structure. 
      %
      % Result: 
      %   1. Replace the surf structure if the name of 
      %   surf_struct mathes with another surf struct 
      %   in the surf struct array in the object.
      %   2. Throws an error if the surf struct is not valid.
      
      
      % check if it is a valid surf structure
      obj.valid_surf(surf_struct);
            
      % Sets an existing surface to the new field values
      match_idx = obj.get_index(surf_struct.name,true);
      obj.surf(match_idx) = surf_struct;     
    end
    
    %% remove_surf
    function remove_surf(obj, surface_name)
      % obj.remove_surf(surface_name_string)
      %
      % Input:  
      %   surface_name_string: A string that matches 
      %   the name of the surface in the surf struct 
      %   array of the object.
      % 
      % Result: 
      %   1. Removes surf struct in the surf struct array 
      %   of the object.
      %   2. Throws an error if the surface_name is not 
      %   in any of the surf struct in the surf struct array 
      %   of the object.
      
      match_idx = obj.get_index(surface_name,true);

      for surf_idx = 1:length(obj.surf)

        if obj.surf(surf_idx).top == match_idx
          obj.surf(surf_idx).top = [];
        elseif obj.surf(surf_idx).top > match_idx
            obj.surf(surf_idx).top = obj.surf(surf_idx).top -1;            
        end

        if obj.surf(surf_idx).active == match_idx
          obj.surf(surf_idx).active = [];
        elseif obj.surf(surf_idx).active > match_idx
          obj.surf(surf_idx).active = obj.surf(surf_idx).active -1;  
        end

        if obj.surf(surf_idx).mask == match_idx
          obj.surf(surf_idx).mask = [];
        elseif obj.surf(surf_idx).mask > match_idx
          obj.surf(surf_idx).mask = obj.surf(surf_idx).mask -1;  
        end

        if obj.surf(surf_idx).gt == match_idx
          obj.surf(surf_idx).gt = [];
        elseif obj.surf(surf_idx).gt > match_idx
          obj.surf(surf_idx).gt = obj.surf(surf_idx).gt -1;  
        end

        if obj.surf(surf_idx).quality == match_idx
          obj.surf(surf_idx).quality = [];
        elseif obj.surf(surf_idx).quality > match_idx
          obj.surf(surf_idx).quality = obj.surf(surf_idx).quality -1;  
        end
      end

      obj.surf = [obj.surf(1:match_idx-1) obj.surf(match_idx+1:end)];
    end
    
    %% adjust_surf
    function adjust_surf(obj, dbin, surface_name_string)
      % obj.adjust_surf(dbin)
      %
      % Input:  
      %   dbin: Number of bins to offset the surface
      %
      %   surface_name_string: A cell array of strings that match the name
      %   of the surface in the surf struct array of the object. If empty
      %   or undefined, all surfaces will be updated.
      % 
      % Result: 
      %   1. surf struct .y field will be adjusted by dbin

      if ~exist('surface_name_string','var')
        surface_name_string = [];
      end
      for surf_idx = 1:length(obj.surf)
        if isempty(surface_name_string) ...
            || ((ischar(surface_name_string) || iscell(surface_name_string)) ...
            && any(strcmpi(obj.surf(surf_idx).name,surface_name_string)))
          imagesc(obj.surf(surf_idx).y)
          colorbar;
          obj.surf(surf_idx).name
          obj.surf(surf_idx).y(1)
          keyboard
%           obj.surf(surf_idx).y = obj.surf(surf_idx).y + dbin;
        end
      end
      
    end
    
    %% save
    function [] = save_surfdata(obj, fn,doa_method_flag)
      % obj.save_surfdata(fn)
      %
      % Input: 
      %  fn: A string that will be the name of the saved surfData file.
      %  doa_method_flag: logical that specifies is this file stores DOA
      %   data or beam forming data. Default is false.
      %
      % Result:
      %   Saves the surf struct array in the object 
      %   as a .mat file. Creates directories as needed.
      
      if ~exist('doa_method_flag','var')
        doa_method_flag = false;
      end
      for surf_idx = 1:length(obj.surf)
        obj.valid_surf(obj.surf(surf_idx),doa_method_flag);
      end
      surf = obj.surf;
      version = obj.current_version;
      try
        surf = rmfield(surf,'h_plot');
      end
      valid_metadata(obj,obj,doa_method_flag);
      radar_name = obj.radar_name;
      season_name = obj.season_name;
      day_seg = obj.day_seg;
      frm = obj.frm;
      gps_time = obj.gps_time;
      theta = obj.theta;
      time = obj.time;
      FCS = obj.FCS;
      
      fn_dir = fileparts(fn);
      if ~exist(fn_dir,'dir')
        mkdir(fn_dir);
      end
      file_version = '1L';
      save(fn, 'surf', 'version', 'gps_time', 'theta', 'time', 'FCS', ...
        'radar_name', 'season_name', 'day_seg', 'frm', 'file_version', '-v7.3');
    end

    %% get_names
    function surf_names = get_names(obj)
      % surf_names = get_names(obj)
      %
      % Returns all the surface names in a cell array
      %
      % No inputs
      %
      % surf_names: cell array of surface names
      
      if ~isempty(obj.surf)
        surf_names = {obj.surf.name};
      else
        surf_names = {};
      end
    end
    
    %% get_index
    function surf_idx = get_index(obj, surf_name, error_on_fail)
      % obj.get_index(surf_name, error_on_fail)
      % 
      % Input: 
      % surf_name: A string that matches the name of a surface in the
      %   surf struct array in the object OR a cell array of similar
      %   strings OR an index array of indices into the surf structure
      % error_on_fail: Logical scalar indicating whether or not an error
      %   should be thrown if the surface is not found. Default is false.
      %
      % Return:
      % 1. The index of that surf struct in the surf 
      %   struct array in the object, if it is found.
      % 2. Throws an error if it is not found and error_on_fail is true
      %   or if an array of indices were input and contain an invalid
      %   index.
      
      if isnumeric(surf_name)
        % Just check that indices are valid
        valid_list = intersect(unique(surf_name),1:length(obj.surf));
        if length(valid_list) ~= length(unique(surf_name))
          error('Some of the indices in surf_name are invalid.');
        end
        surf_idx = surf_name;
        
      else
        if ~exist('error_on_fail','var') || isempty(error_on_fail)
          error_on_fail = false;
        end
        
        if ischar(surf_name)
          surf_name = {surf_name};
        elseif ~iscell(surf_name)
          error('Input type of surf_name is not valid.');
        end
        
        for idx = 1:length(surf_name)
          new_surf_idx = find(strcmp(surf_name{idx},{obj.surf.name}));
          if error_on_fail && isempty(new_surf_idx)
            error('Surface "%s" is not in the list of surfaces.', surf_name{idx});
          end
          surf_idx(idx) = new_surf_idx;
        end
      end
    end
   
    
    %% valid_surf
    function [] = valid_surf(obj, surf_struct,doa_method_flag)
      % obj.valid_surf(surf_struct)
      %
      % Input:
      %  surf_struct: A surf structure.
      %  doa_method_flag: logical that specifies is this file stores DOA
      %   data or beam forming data. Default is false. Field not used.
      %
      % Result: 
      % 1. Nothing, if the surf structure is valid.
      % 2. Throws an error if any of the condition 
      % is not satisfied.
      
      if ~exist('doa_method_flag','var')
        doa_method_flag = false;
      end
      
      % check if it's a structure
      if ~isstruct(surf_struct)
        error('The input type is not a structure');
      end      
      
      % check if all the fields exist
        if ~isfield(surf_struct, {'x', 'y', 'plot_name_values', ...
                          'name', 'top', 'active', ...
                          'mask', 'gt', ...
                          'quality','visible'})
           error('Struct missing field(s) that are necessary for a surface struct.');
        end
      
      % type check
        if ~isa(surf_struct.x, 'double')
          error('Invalid field type for the field x (Should be type double)');
        end
        
        if ~isa(surf_struct.y, 'double') && (~isa(surf_struct.y, 'logical'))
          error('Invalid field type for the field y (Should be type double or logical)');
        end

      if ~isa(surf_struct.plot_name_values, 'cell')
        error('Invalid field type for the field plot_name_values (Should be type cell).');
      end
      
      if ~isa(surf_struct.name, 'char')
        error('Invalid field type for the field name (Should be type char).');
      end
   
      if ~isa(surf_struct.top, 'double')
        error('Invalid field type for the field top (Should be type double).');
      end

      if ~isa(surf_struct.active, 'double')
        error('Invalid field type for the field active (Should be type double).');
      end
      
      if ~isa(surf_struct.mask, 'double')
        error('Invalid field type for the field mask (Should be type double).');
      end
      
      if ~isa(surf_struct.gt, 'double')
        error('Invalid field type for the field gt (Should be type double).');
      end
           
      if ~isa(surf_struct.quality, 'double')
        error('Invalid field type for the field quality (Should be type double).');
      end
      
      if ~isa(surf_struct.visible, 'logical')
        error('Invalid field type for the field visible (Should be type logical).');
      end
            
      % Check that size of x and y match each other
        if ~isequal(size(surf_struct.x), size(surf_struct.y))
          error('Size of the field x and field y do not match each other.');
        end
      % check if all the surface_struct in the object has the same dimension in x and y
      if ~isempty(obj.surf) && ...
          ~all(arrayfun(@(arr) isequal(size(surf_struct.x), size(arr.x)), obj.surf))
        % we only check x because x and y must have the same dimension
        % from the previous check.
        error('Field x and y sizes do not match this surfdata''s surf x and y sizes.');
      end
      % Check that surface indices point to valid surface indices:
      if ~isempty(surf_struct.top) ...
          && ~(surf_struct.top > 0 || surf_struct.top < length(obj.surf))
        error('top index out of range.');
      end
      
      if ~isempty(surf_struct.active) ...
          && ~(surf_struct.active > 0 || surf_struct.active < length(obj.surf))
        error('active index out of range.');
      end
      
      if ~isempty(surf_struct.mask) ...
          && ~(surf_struct.mask > 0 || surf_struct.mask < length(obj.surf))
        error('mask index out of range.');
      end
      
      if ~isempty(surf_struct.gt) ...
          && ~(surf_struct.gt > 0 || surf_struct.gt < length(obj.surf))
        error('gt index out of range.');
      end
      
      if ~isempty(surf_struct.quality) ...
          && ~(surf_struct.quality > 0 || surf_struct.quality < length(obj.surf))
        error('quality index out of range.');
      end
    end
    
    function [rmse,mean_diff,median_diff,min_diff,max_diff,surf_diff] ...
        = compare(obj, ref, sd_other, other, DOA_trim)
      % [rmse,mean_diff,median_diff,min_diff,max_diff,surf_diff] ...
      %   = tomo.compare(ref, sd_other, other)
      %
      % Compares one of the surfaces in this object to a surface in another
      % surfdata object.
      %
      % ref: reference layer (get_surf's surf_name argument)
      % sd_other: another surfdata object to compare to (can also be this
      %   surfdata object)
      % other: layer to compare (get_surf's surf_name argument)
      % DOA_trim: 
      %
      % rmse: root mean squared error of difference
      % mean_diff: mean of the absolute value of the difference
      % median_diff: median of the absolute value of the difference
      % min_diff: minimum of the absolute value of the difference
      % max_diff: maximum of the absolute value of the difference
      % surf_diff: ansolute vlaue difference matrix
      %
      % See also: tomo.run_compare_surfdata, tomo.compare_surfdata,
      %   tomo.surfdata
      
      ref = obj.get_surf(ref);
      other = sd_other.get_surf(other);
      
      if 1
      quality_mtx_ref = obj.surf(other.quality).y;
      quality_mtx_other = sd_other.surf(other.quality).y;
      
      bad_idx_ref = find(quality_mtx_ref==0);
      bad_idx_other = find(quality_mtx_other==0);
      
      ref.y(bad_idx_ref) = NaN;
      other.y(bad_idx_ref) = NaN;
      
      ref.y(bad_idx_other) = NaN;
      other.y(bad_idx_other) = NaN;
%       if any(isnan(ref.y(:))) || any(isnan(other.y(:)))
%         keyboard
%       end
      end
      
      surf_diff = abs(other.y(1+DOA_trim(1):end-DOA_trim(end)+1,:) ...
        - ref.y(1+DOA_trim(1):end-DOA_trim(end)+1,:));
      rmse        = sqrt(nanmean(abs(surf_diff(:)).^2));
      mean_diff   = nanmean(surf_diff(:));
      median_diff = nanmedian(surf_diff(:));
      min_diff    = nanmin(surf_diff(:));
      max_diff    = nanmax(surf_diff(:));
    end
  
  end
  
  methods(Static)
    function surf = clear_references(surf)
      for idx = 1:length(surf)
        surf(idx).top = [];
        surf(idx).active = [];
        surf(idx).mask = [];
        surf(idx).gt = [];
        surf(idx).quality = [];
      end
    end
        
    function surf = empty_surf(fn)
      % tomo.surfdata.empty_surf(fn)
      %
      % Returns an empty surf structure
      if exist('fn','var') && islogical(fn) && fn
        surf.theta = [];
      else
        surf.x = [];
        surf.y = [];
      end
      surf.plot_name_values = {'color','blue','marker','^'};
      surf.name = '';
      surf.top = [];
      surf.active = [];
      surf.mask = [];
      surf.gt= [];
      surf.quality = [];
      surf.visible = true;
      
    end
    
    function surfdata_ver = version(fn)
      % tomo.version(fn)
      %
      % Returns the version of the file contained in fn.
      
      tmp = load(fn);
      if isfield(tmp,'version')
        if tmp.version == 2.0
          surfdata_ver = tmp.version;
          return;
        else
          error('Invalid file version: %d.\n', tmp.version);
        end
      else
        surfdata_ver = 1.0;
        return;
      end
      
    end
    
    function surf = update_file(fn,fn_new,echogram_fn)
      % tomo.update_file(fn,fn_new,echogram_fn)
      %
      % Updates filename fn to latest version of surfdata and stores the
      % result in fn_new. The updated surf structure is returned.
      %
      % echogram_fn: the echoram filename is only required when converting
      %   from version 1.0 files. The echogram filename fields must match
      %   the surfData file or the updating will not work properly.
      
      fn_version = tomo.surfdata.version(fn);
      if fn_version == 1.0
        updated = load(fn);
        for surf_idx = 1:length(updated.surf)
          if strcmpi(updated.surf(surf_idx).name,'ice surface')
            updated.surf(surf_idx).name = 'top';
          elseif strcmpi(updated.surf(surf_idx).name,'surface gt')
            updated.surf(surf_idx).name = 'top gt';
          elseif strcmpi(updated.surf(surf_idx).name,'surface quality')
            updated.surf(surf_idx).name = 'top quality';
          end
          updated.surf(surf_idx).top = updated.surf(surf_idx).surf_layer;
          updated.surf(surf_idx).active = updated.surf(surf_idx).active_layer;
          updated.surf(surf_idx).mask = updated.surf(surf_idx).mask_layer;
          updated.surf(surf_idx).gt = updated.surf(surf_idx).control_layer;
          updated.surf(surf_idx).quality = updated.surf(surf_idx).quality_layer;
          if ischar(updated.surf(surf_idx).visible)
            if strcmpi(updated.surf(surf_idx).visible,'on')
              updated.surf(surf_idx).visible = true;
            else
              updated.surf(surf_idx).visible = false;
            end
          end
        end
        updated.surf = rmfield(updated.surf,{'surf_layer','active_layer','mask_layer','control_layer','quality_layer'});
        updated.version = 2.0;
        
        tmp = load(echogram_fn,'GPS_time','theta','Time','param_combine');
        updated.radar_name = tmp.param_combine.radar_name;
        updated.season_name = tmp.param_combine.season_name;
        updated.day_seg = tmp.param_combine.day_seg;
        updated.frm = tmp.param_combine.combine.frm;
        updated.gps_time = tmp.GPS_time;
        updated.theta = tmp.theta(:); % Make a column vector
        updated.time = tmp.Time(:); % Make a column vector
        updated.FCS.origin = tmp.param_combine.array_param.fcs{1}{1}.origin;
        updated.FCS.x = tmp.param_combine.array_param.fcs{1}{1}.x;
        updated.FCS.y = tmp.param_combine.array_param.fcs{1}{1}.y;
        updated.FCS.z = tmp.param_combine.array_param.fcs{1}{1}.z;
        
      elseif fn_version == 2.0
        updated = load(fn);
        
      else
        error('Invalid file version: %d.\n', fn_version);
        
      end
      
      if nargin >= 2
        obj = tomo.surfdata;
        obj.load_surfdata(updated);
        obj.save_surfdata(fn_new);
      end
      
      surf = updated.surf;
      
    end
    
    function [] = run_add_surf_from_dem()
      %% run add surf from dem
      % script run_add_surf_from_dem
      %
      % script for running tomo.surfdata.add_surf_from_dem
      %
      % Author:  Theresa Moore
      %
      % See also: tomo.surfdata
      
      
      % User Setup
      % =========================================================================
      param_override = [];
      params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'));
      params = ct_set_params(params,'cmd.generic',0);
      params = ct_set_params(params,'cmd.generic',1,'day_seg','20140506_01');
      params = ct_set_params(params,'cmd.frms',[2 3 4]);
%       params = ct_set_params(params,'add_surf_from_dem.ice_mask_fn',fullfile('greenland','IceMask','GimpIceMask_90m_v1.1.tif'));%'antarctica\DEM\BEDMAP2\original_data\bedmap2_tiff\bedmap2_icemask_grounded_and_shelves.tif';
      params = ct_set_params(params,'add_surf_from_dem.ice_mask_fn','canada/ice_mask/03_rgi50_ArcticCanadaNorth/03_rgi50_ArcticCanadaNorth.mat');%'antarctica\DEM\BEDMAP2\original_data\bedmap2_tiff\bedmap2_icemask_grounded_and_shelves.tif';
%       params = ct_set_params(params,'add_surf_from_dem.ice_mask_fn',ct_filename_gis(params,fullfile('greenland','IceMask','GimpIceMask_90m_v1.1.tif')));%'antarctica\DEM\BEDMAP2\original_data\bedmap2_tiff\bedmap2_icemask_grounded_and_shelves.tif';
%       param.add_surf_from_dem.ice_mask_fn = ct_filename_gis(param,fullfile('greenland','IceMask','GimpIceMask_90m_v1.1.tif'));%'antarctica\DEM\BEDMAP2\original_data\bedmap2_tiff\bedmap2_icemask_grounded_and_shelves.tif';
%         param.add_surf_from_dem.ice_mask_fn = ct_filename_gis([],'canada/ice_mask/03_rgi50_ArcticCanadaNorth/03_rgi50_ArcticCanadaNorth.mat');

      params = ct_set_params(params, 'add_surf_from_dem.dem_guard', 30e3);
%       params = ct_set_params(params, 'add_surf_from_dem.dem_per_slice_guard', 500e3);
      
      % Automated Section
      % =========================================================================
      % Input checking
      global gRadar;
      if exist('param_override','var')
        param_override = merge_structs(gRadar,param_override);
      else
        param_override = gRadar;
      end
      
      % Process each of the segments
      ctrl_chain = {};
      for param_idx = 1:length(params)
        param = params(param_idx);
        if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
          continue;
        end
        tomo.surfdata.add_surf_from_dem(param,param_override);
      end
    end
    
    function surf = add_surf_from_dem(param,param_override)
      %% add surf from dem
      % tomo.add_surf_from_dem(param,param_override)
      %
      % param:  struct with processing parameters or function handle to
      % script with processing parameters
      %
      % param_override: parameters in this struct override parameters in
      % param.  This struct must also contain the gRdar fields.
      % Typically global gRdadr; param_override = gRadar;
      %
      % Example:
      %   See run_add_surf_from_dem for how to run this function
      %   directly.  This function may be called from the run_master.m
      %   script using the param spreadsheet and enabling the cmd.generic
      %   column.
      %
      % Authors:  Theresa Moore, John Paden
      %
      % See also:  run_add_surf_from_dem
      %
      %
      %% General Setup
      % =================================================================
      physical_constants
      param = merge_structs(param,param_override);
      
      fprintf('==========================================================\n');
      fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
      fprintf('==========================================================\n');
      
      %% Input Checks
      % =================================================================
      
      % Load frames file
      load(ct_filename_support(param, '', 'frames'));
      
      % If no frames specified, then do all frames
      if isempty(param.cmd.frms)
        param.cmd.frms = 1:length(frames.frame_idxs);
      end
      
      % Remove frames that do not exist from param.cmd.frms list
      [valid_frms,keep_idxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
      if length(valid_frms) ~= length(param.cmd.frms)
        bad_mask = ones(size(param.cmd.frms));
        bad_mask(keep_idxs) = 0;
        warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
          param.cmd.frms(find(bad_mask,1)));
        param.cmd.frms = valid_frms;
      end
      
      
      if ~isfield(param,'add_surf_from_dem')
        param.add_surf_from_dem = [];
      end
      
      if ~isfield(param.add_surf_from_dem,'dem_res') || isempty(param.add_surf_from_dem.dem_res)
        param.add_surf_from_dem.dem_res = 10;
      end
      
      if ~isfield(param.add_surf_from_dem,'theta_vec') || isempty(param.add_surf_from_dem.theta_vec)
        param.add_surf_from_dem.theta_vec = -85:85;
      end
      
      if ~isfield(param.add_surf_from_dem,'dem_guard') || isempty(param.add_surf_from_dem.dem_guard)
        param.add_surf_from_dem.dem_guard = 12e3;
      end
      
      if ~isfield(param.add_surf_from_dem,'ice_mask_fn') || isempty(param.add_surf_from_dem.ice_mask_fn)
        param.add_surf_from_dem.ice_mask_fn = [];
      end
      
      if ~isfield(param.add_surf_from_dem,'sv_cal_fn') || isempty(param.add_surf_from_dem.sv_cal_fn)
        param.add_surf_from_dem.sv_cal_fn = [ct_filename_ct_tmp(param,'','add_surf_from_dem','theta_cal') '.mat'];
      end
      
      % .dem_per_slice_guard: additional region in meters around each slice to search for DEM points
      %   Setting too high slows the process down, setting too low will miss
      %   DEM points needed to properly represent the surface.
      if ~isfield(param.add_surf_from_dem,'dem_per_slice_guard') || isempty(param.add_surf_from_dem.dem_per_slice_guard)
        param.add_surf_from_dem.dem_per_slice_guard = 240;
      end
      
      if ~isfield(param.add_surf_from_dem,'surfdata_mode') || isempty(param.add_surf_from_dem.surfdata_mode)
        param.add_surf_from_dem.surfdata_mode.overwrite = true;
      end
      
       
      if ~isfield(param.add_surf_from_dem,'surf_out_path') || isempty(param.add_surf_from_dem.surf_out_path)
        param.add_surf_from_dem.surf_out_path = '';
      end
        
      % Set doa_method_flag to false (always)
      doa_method_flag = false;
      
    
      %% Setup processing
      % =====================================================================
      
      % Load in sar coordinates file
      % Currently supports default path only -- FIX THIS??
      sar_coord_fn = fullfile(ct_filename_out(param,'sar',''),'sar_coord.mat');
      sar_coord = load(sar_coord_fn);

      % Load records file
      records_fn = ct_filename_support(param,'','records');
      records = load(records_fn);
      
      
      % Loop over frames
      for frm_idx = 1:length(param.cmd.frms)
        frm = param.cmd.frms(frm_idx);
        
        fprintf('Creating surface data object for %s_%03d \n', param.day_seg,frm); 
        
        % recs: Determine the records for this frame
        if frm < length(frames.frame_idxs)
          recs = [frames.frame_idxs(frm), frames.frame_idxs(frm+1)-1];
        else
          recs = [frames.frame_idxs(frm), length(records.gps_time)];
        end
        
        
        % The following is for Debug only
        % res(2) = recs(1) + 20000;
               
        % Find gps at frame boundaries and use this to mask out FCS of a
        % frame
        frm_gps_start = records.gps_time(recs(1));
        frm_gps_stop  = records.gps_time(recs(end));
        
        % Mask used to isolate FCS for a single frame
        gps_mask = sar_coord.gps_time >= frm_gps_start & ... 
          sar_coord.gps_time <= frm_gps_stop;
        
        % GPS time, origin, Xpos, Ypos, Zpos for the frame
        frm_gps_time  = sar_coord.gps_time(gps_mask);
        frm_origin    = sar_coord.origin(:,gps_mask);
        frm_x         = sar_coord.x(:,gps_mask);
        frm_z         = sar_coord.z(:,gps_mask);
        frm_y         = cross(frm_z,frm_x);
        
        % Along track decimation factor
        delta_at      = 10;
        
        % Create surfdata
        sd = tomo.surfdata();
        sd.radar_name = param.radar_name;
        sd.season_name = param.season_name;
        sd.day_seg = param.day_seg;
        sd.frm = frm;
        sd.gps_time = frm_gps_time(1:delta_at:end);       
        sd.FCS.origin = frm_origin(:,1:delta_at:end);
        sd.FCS.x = frm_x(:,1:delta_at:end);
        sd.FCS.z = frm_z(:,1:delta_at:end);
        sd.FCS.y = frm_y(:,1:delta_at:end);
        sd.theta = param.add_surf_from_dem.theta_vec(:);
        sd.time = zeros(0,1);
        
        clear frm_gps_time frm_origin frm_x frm_y frm_z
        
        %% Load Geotiff and Ice Mask
        % =================================================================

        global gdem;
        if isempty(gdem) || ~isa(gdem,'dem_class') || ~isvalid(gdem)
          gdem = dem_class(param,param.add_surf_from_dem.dem_res);
        end
        gdem.set_res(param.add_surf_from_dem.dem_res);
        
        % Load ice mask        
        if isfield(param.add_surf_from_dem,'ice_mask_fn') && ~isempty(param.add_surf_from_dem.ice_mask_fn)
          ice_mask_fn = ct_filename_gis(param,param.add_surf_from_dem.ice_mask_fn);
          [~,ice_mask_fn_name,ice_mask_fn_ext] = fileparts(ice_mask_fn);
          if strcmpi(ice_mask_fn_ext,'.tif')
            ice_mask_all.proj = geotiffinfo(ice_mask_fn);
            [ice_mask_all.mask, R, ~] = geotiffread(ice_mask_fn);
            ice_mask_all.X = R(3,1) + R(2,1)*(1:size(ice_mask_all.mask,2));
            ice_mask_all.Y = R(3,2) + R(1,2)*(1:size(ice_mask_all.mask,1));
          else
            ice_mask_all = load(ice_mask_fn);
          end
        else
          ice_mask_all = [];
        end
        
        % sv_cal_fn: steering vector calibration filename
        sv_cal_fn = param.add_surf_from_dem.sv_cal_fn;
        if ~exist(sv_cal_fn,'file')
          sv_cal_fn = [];
        end
        
        
        %% Create DEM
        % Convert decimated origin coordinates from ECEF to geodetic and
        % use these to define boundaries of the dem
        Nx = size(sd.FCS.x,2);
        dec_idxs = round(linspace(1,Nx,min(Nx,200)));
        [frm_lat_dec, frm_lon_dec, frm_elev_dec] = ecef2geodetic(sd.FCS.origin(1,dec_idxs),sd.FCS.origin(2,dec_idxs), sd.FCS.origin(3,dec_idxs),WGS84.ellipsoid);
        frm_lat_dec = frm_lat_dec*180/pi;
        frm_lon_dec = frm_lon_dec*180/pi;
        
        [latb,lonb] = bufferm(frm_lat_dec,frm_lon_dec,param.add_surf_from_dem.dem_guard/WGS84.semimajor*180/pi);
        gdem_str = sprintf('%s:%s:%s_%03d',param.radar_name,param.season_name,param.day_seg,frm);
        
        if ~strcmpi(gdem_str,gdem.name)
          gdem.set_vector(latb,lonb,gdem_str);
        end
        
        % Generated a non-decimated version of the lat,lon for the frame -
        % used below
        [frm_lat, frm_lon, frm_elev] = ecef2geodetic(sd.FCS.origin(1,:),sd.FCS.origin(2,:), sd.FCS.origin(3,:),WGS84.ellipsoid);
        frm_lat = frm_lat * 180/pi;
        frm_lon = frm_lon * 180/pi;

        gdem.set_vector(latb,lonb,gdem_str);
        [DEM,msl,ocean_mask,proj,DEM_x,DEM_y] = gdem.get_vector_mosaic(100);
        DEM(ocean_mask) = msl(ocean_mask);
        [frm_x,frm_y] = projfwd(proj,frm_lat,frm_lon);
        
        %% Interpolate at all of the bad value locations using the good data
        bad_idxs = find(isnan(DEM));
        good_idxs = find(~isnan(DEM));
        x_idxs = repmat(1:size(DEM,2),[size(DEM,1) 1]);
        y_idxs = repmat((1:size(DEM,1))',[1 size(DEM,2)]);
        x_vals = x_idxs(good_idxs);
        y_vals = y_idxs(good_idxs);
        z_vals = DEM(good_idxs);
        x_out = x_idxs(bad_idxs);
        y_out = y_idxs(bad_idxs);
        z_out = single(griddata(x_vals,y_vals,double(z_vals),x_out,y_out));
        if ~isempty(z_out)
          DEM(bad_idxs) = z_out;
        end
        
        if 0
          figure(1);clf;
          imagesc(DEM_x,DEM_y,DEM);
          hold on;
          plot(frm_x,frm_y,'r');
          keyboard;
        end
        
        DEM_x_mesh = repmat(DEM_x,[size(DEM,1) 1]);
        DEM_y_mesh= repmat(DEM_y,[1 size(DEM,2)]);
        
        
        %% For every position along the aperture, add theta dependent TWTT to DEM from origin of the FCS   
        Nx    = size(sd.FCS.x,2);
        Nsv   = length(param.add_surf_from_dem.theta_vec);
        twtt  = zeros(Nsv,Nx);
        ice_mask = NaN(Nsv,Nx);
        
        % Handle case where the DEM is all NaNs
        if all(all(isnan(DEM)))
          warning('Input DEM contains all NaN data for Frame %d.',param.proc.frm);
          twtt(:,:) = NaN;
          Nx = 0;
        end
        
        DEM_coverage_warning = false;
        
        for rline = 1:Nx
          if ~mod(rline-1,10^floor(log10(Nx)-1))
            fprintf('  %s %d of %d (%s)\n', mfilename, rline, Nx, datestr(now));
          end
          
          dem_guard = param.add_surf_from_dem.dem_guard;
          
          DEM_mask = DEM_x_mesh > frm_x(rline)-dem_guard & DEM_x_mesh < frm_x(rline)+dem_guard ...
            & DEM_y_mesh > frm_y(rline)-dem_guard & DEM_y_mesh < frm_y(rline)+dem_guard ...
            & ~isnan(DEM);
          DEM_idxs = find(DEM_mask);
          
          if numel(DEM_idxs)==0
            warning('Range Line %d of Frame %d is not spanned by DEM.',rline,frm);
          end
          
          if 0
            set(h_img,'AlphaData',DEM_mask);
          end
          
          % Convert from projection to geodetic (lat,lon,elev)
          [DEM_lat,DEM_lon] = projinv(proj,DEM_x_mesh(DEM_idxs),DEM_y_mesh(DEM_idxs));
          DEM_elev = DEM(DEM_idxs);
          
          % Convert from geodetic (lat,lon,elev) to ECEF (x,y,z)
          physical_constants;
          [DEM_ecef_x,DEM_ecef_y,DEM_ecef_z] = geodetic2ecef(single(DEM_lat)/180*pi,single(DEM_lon)/180*pi,single(DEM_elev),WGS84.ellipsoid);
          
          origin = sd.FCS.origin(:,rline);
          
          % Define matrices to change coordinate systems from FCS to ECEF
          % and vice versa
          Tfcs_ecef = [sd.FCS.x(:,rline), sd.FCS.y(:,rline), sd.FCS.z(:,rline)];
          Tecef_fcs = inv(Tfcs_ecef);
          
          % Convert DEM coordinates from ECEF to FCS
          tmp = Tecef_fcs * [DEM_ecef_x.'-origin(1); DEM_ecef_y.'-origin(2); DEM_ecef_z.'-origin(3)];
          DEM_fcs_x = tmp(1,:);
          DEM_fcs_y = tmp(2,:);
          DEM_fcs_z = tmp(3,:);
          
          if 0
            imagesc(reshape(DEM_fcs_x,[200 200]))
            colorbar;
            
            imagesc(reshape(DEM_fcs_y,[200 200]))
            colorbar;
            
            imagesc(reshape(DEM_fcs_z,[200 200]))
            colorbar;
          end
          
          slice_mask = DEM_fcs_x > -param.add_surf_from_dem.dem_per_slice_guard & DEM_fcs_x < param.add_surf_from_dem.dem_per_slice_guard;
          
          x = DEM_fcs_x(slice_mask);
          y = DEM_fcs_y(slice_mask);
          z = DEM_fcs_z(slice_mask);
          
          if (numel(x)>=3)
            faces = delaunay(double(x),double(y));
            vertices = [double(x).' double(y).' double(z).'];  % vertices stored as Nx3 matrix
            vert1 = vertices(faces(:,1),:);
            vert2 = vertices(faces(:,2),:);
            vert3 = vertices(faces(:,3),:);
            
            orig = [0 0 0];
            
            theta_rline = param.add_surf_from_dem.theta_vec*(pi/180);
         
            intersection = zeros(3,Nsv);
            
            for theta_idx = 1:length(theta_rline)
              dir = [0 sin(theta_rline(theta_idx)) -cos(theta_rline(theta_idx))];
                            
              [Intersect, t] = TriangleRayIntersection(orig, dir, vert1, vert2, vert3);
              
              intersect_idx = find(Intersect);
              if isempty(intersect_idx)
                twtt(theta_idx,rline) = NaN;
                intersection(:,theta_idx) = NaN;
              else
                twtt(theta_idx,rline) = t(intersect_idx(1))/(3e8/2);
                intersection(:,theta_idx) = mean([vert1(intersect_idx(1),:);vert2(intersect_idx(1),:);vert3(intersect_idx(1),:)],1);
              end
              
%               if 0
%                 figure(999);clf
%                 D = c*twtt(theta_idx,rline)/2;
%                 Dvec = D*dir;
%                 Dx = Dvec(1);
%                 Dy = Dvec(2);
%                 Dz = Dvec(3);
%                 h = trisurf(faces,x,y,z, Intersect*1.0, 'FaceAlpha',0.9); 
%                 set(h,'edgecolor','none')
%                 hold on
%                 line('XData',orig(1) + [0 Dx], 'YData',orig(2) + [0 Dy], 'ZData', orig(3) + [0 Dz], 'Color', 'r', 'LineWidth',3)
%                 xlabel('X_F_C_S (m)')
%                 ylabel('Y_F_C_S (m)')
%                 zlabel('Z_F_C_S (m)')
%                 
%               end
            end
          else
            if ~DEM_coverage_warning
              DEM_coverage_warning = true;
              warning('DEM dem_per_slice_guard too small.');
            end
            clear intersection;
            twtt(:,rline) = NaN;
          end
          
          
          if exist('ice_mask_all','var') && ~isempty(ice_mask_all)
            if exist('intersection','var')
              % Convert from FCS/SAR to ECEF
              intersection_ecef = Tfcs_ecef * intersection;
              intersection_ecef_x = intersection_ecef(1,:).' + origin(1);
              intersection_ecef_y = intersection_ecef(2,:).' + origin(2);
              intersection_ecef_z = intersection_ecef(3,:).' + origin(3);
              % Convert from ECEF to geodetic
              [intersection_lat,intersection_lon,tri_h] = ecef2geodetic(intersection_ecef_x,intersection_ecef_y,intersection_ecef_z,WGS84.ellipsoid);
              intersection_lat = intersection_lat*180/pi;
              intersection_lon = intersection_lon*180/pi;
              % Convert from geodetic to projection
              [intersection_x,intersection_y] = projfwd(ice_mask_all.proj,intersection_lat,intersection_lon);
              % Get mask coordinates nearest triangle center coordinates
              intersection_x_idx = interp1(ice_mask_all.X,1:length(ice_mask_all.X),intersection_x,'nearest');
              intersection_y_idx = interp1(ice_mask_all.Y,1:length(ice_mask_all.Y),intersection_y,'nearest');
              % Find nan values and set to integer value
              nidx = find(isnan(intersection_x_idx));
              intersection_x_idx(nidx) = 1;
              intersection_y_idx(nidx) = 1;
              % Convert triangle mask coordinates to matrix indices
              mask_idx = (intersection_x_idx-1)*length(ice_mask_all.Y) + intersection_y_idx;
              % Find ice mask for triangle coordinates
              ice_mask(:,rline) = ice_mask_all.mask(mask_idx);
              % Set previously nan valued coordinates to 0 mask
%               ice_mask(nidx,rline) = 0;
            else
              %         for ice_mask_idx = 1:length(mask_idx)
              %           ice_mask(theta_rline_row_idx(ice_mask_idx),theta_rline_col_idx(ice_mask_idx),rline) = ice_mask_all.mask(mask_idx(ice_mask_idx));
              %         end
              %         % Set previously nan valued coordinates to 0 mask
              %         ice_mask(theta_rline_row_idx(nidx),theta_rline_col_idx(nidx),rline) = 0;
              %       end
              
            end
          end
          
        end
        
        %% Insert Surfaces
        
        % Insert TWTT of top surface
        surf = tomo.surfdata.empty_surf();
        surf.x = repmat((1:Nsv).',[1 size(twtt,2)]); % Nsv x Nx

        surf.y = twtt;
        surf.plot_name_values = {'color','black','marker','x'};
        surf.name = 'top twtt';
        sd.insert_surf(surf);
        
        % Insert icemask
          surf = tomo.surfdata.empty_surf();
        surf.x = repmat((1:Nsv).',[1 size(twtt,2)]); % Nsv x Nx
        surf.y = ice_mask;
        surf.plot_name_values = {'color','white','marker','x'};
        surf.name = 'ice_mask';
        sd.insert_surf(surf);
        
        %% Create output filename
        out_dir = ct_filename_out(param,param.add_surf_from_dem.surf_out_path,'surfData_sar');
        
        if ~isdir(out_dir)
          mkdir(out_dir);
        end
        
        out_fn_name = sprintf('Data_%s_%03.0f.mat',param.day_seg,frm);
        out_fn = fullfile(out_dir,out_fn_name);
        sd.save_surfdata(out_fn,doa_method_flag);
        fprintf('Done (%s)\n', datestr(now));
        
        
      end
      
      
      
      
   
    end
    
  end
  
end

