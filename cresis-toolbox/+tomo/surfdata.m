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
        try
          obj.surf = [];
          obj.load_surfdata(fn);
          
        catch ME
          error('Failed to load filename (%s):\n %s', fn, ME.getReport());
        end
      end
    end
    
    %% load_surf
    function load_surfdata(obj, fn)
      % obj.load_surfdata(fn)
      %
      % Loads a new surf struct into the object. If fn is a filename,
      % it should be a surfData file. Note that this function will overwrite
      % the existing surf struct array.
      %
      % Input:
      %   fn: The file path OR a new struct array with all the fields in a
      %   surfData file.
      
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
          obj.valid_surf(obj.surf(surf_idx));
        end
        
        obj.set_metadata(tmp);
        
      catch ME
        obj.surf = old_surf;
        if ischar(fn)
          error('Failed to load %s:\n %s', fn, ME.getReport());
        else
          error('Failed to load struct:\n %s', ME.getReport());
        end
      end
    end
    
    function set_metadata(obj, md)
      % obj.set_metadata(md)
      % 
      % Function for validating all the metadata fields
      %
      % md: struct containing gps_time, theta, FCS, radar_name,
      %   season_name, day_seg, and frm.
      
      obj.valid_metadata(md);
      
      obj.radar_name = md.radar_name;
      obj.season_name = md.season_name;
      obj.day_seg = md.day_seg;
      obj.frm = md.frm;
      obj.gps_time = md.gps_time;
      obj.theta = md.theta;
      obj.time = md.time;
      obj.FCS = md.FCS;
    end
    
    function valid_metadata(obj, md)
      % obj.valid_metadata(md)
      % 
      % Function for validating all the metadata fields
      %
      % md: struct containing gps_time, theta, FCS, radar_name,
      %   season_name, day_seg, and frm.
      
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
    function insert_surf(obj, surf_struct)
      %  obj.insert_surf(surf_struct)
      %
      %  Input: 
      %    surf_struct: A surf structure. Normally, the indexing fields
      %    such as top, active, mask, gt, quality should be cleared.
      %
      %  Result: 
      %    Adds a surf structure, surf_A, into 
      %    the surf array into the object
      %
      % See also: surfdata.clear_references
      
      % check if it is a valid surf structure
      obj.valid_surf(surf_struct);
      
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
    function [] = save_surfdata(obj, fn)
      % obj.save_surfdata(fn)
      %
      % Input: 
      %   fn: A string that will be the name of the saved surfData file.
      %
      % Result:
      %   Saves the surf struct array in the object 
      %   as a .mat file. Creates directories as needed.
      
      for surf_idx = 1:length(obj.surf)
        obj.valid_surf(obj.surf(surf_idx));
      end
      surf = obj.surf;
      version = obj.current_version;
      try
        surf = rmfield(surf,'h_plot');
      end
      valid_metadata(obj,obj);
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
      
      save(fn, 'surf', 'version', 'gps_time', 'theta', 'time', 'FCS', ...
        'radar_name', 'season_name', 'day_seg', 'frm', '-v7.3');
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
    function [] = valid_surf(obj, surf_struct)
      % obj.valid_surf(surf_struct)
      %
      % Input: A surf structure.
      %
      % Result: 
      % 1. Nothing, if the surf structure is valid.
      % 2. Throws an error if any of the condition 
      % is not satisfied.
      
      % check if it's a structure
      if ~isstruct(surf_struct)
        error('The input type is not a structure');
      end      
      
      % check if all the fields exist
      if  ~isfield(surf_struct, {'x', 'y', 'plot_name_values', ...
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
      
      surf.x = [];
      surf.y = [];
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
    
  end
  
end

