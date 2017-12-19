classdef surf_data < handle
  % 
  % A class that manages surf data in a frame. 
  % 
  % Constructor: surf_data(filepath)
  %
  % Example of a filepath:
  % X:\ct_data\rds\2014_Greenland_P3\CSARP_shane_music_surfData\20140401_03\Data_20140401_03_037.mat
  %
  % Example of how to generate a filepath:  
  %   param.radar_name   = 'rds';
  %   param.season_name  = '2014_Greenland_P3';
  %   param.day_seg = '20140401_03'
  %   surfdata_ref = 'shane_music_surfData';
  %   frame_idx = 37;
  %   fullfile(ct_filename_out(param,surfdata_ref,''),...
  %     sprintf('Data_%s_%03d.mat',param.day_seg,frame_idx));      
  %
  % Author: Shane Chu, John Paden
  
  
  properties
    surf
  end
  
  methods    
    %% constructor   
    function obj = surf_data(filename)
      % no arguments => create an empty surf struct array
      if nargin == 0
        obj.surf = [];
      % one argument => load file and create struct array accordingly
      elseif nargin == 1 && isa(filename, 'char') 
        loaded_surf_struct = load(filename);        
        obj.surf = loaded_surf_struct.surf;
        
        for surf_idx = 1:length(obj.surf)
          obj.check_if_surf_valid(obj.surf(surf_idx));
        end                 
      else
        error('Please enter a valid filename.');
      end
    end
    
    %% load_surf
    function obj = load_surf(obj, filename)
    % obj.load_surf(filename)
    %
    % Input: 
    %   filename: The file path.
    %   e.g. X:\ct_data\rds\2014_Greenland_P3\CSARP_shane_music_surfData\20140401_03\Data_20140401_03_037.mat
    %
    % Note: You can check ct_filename_out in the cresis 
    % toolbox for how to generate the file path.
    %
    % Result: 
    %   Loads the entire surf struct array from a file 
    %   that contains surface struct array.
    %
    % Note: It will overwirte the existing surface struct array.
        
      if isa(filename, 'char') 
        loaded_surf_struct = load(filename);
        obj.surf = loaded_surf_struct.surf;
        
        for surf_idx = 1:length(loaded_surf_struct.surf)
          obj.check_if_surf_valid(obj.surf(surf_idx));          
        end           
        
      else
        error('Please enter a valid filename.');
      end
    end  
    
    %% insert_surf
    function obj = insert_surf(obj, surf_struct)
      %  obj.insert_surf(surf_struct)
      %
      %  Input: 
      %    surf_struct: A surf structure. 
      %
      %  Result: 
      %    Adds a surf structure, surf_A, into 
      %    the surf array into the object
      
      % check if it is a valid surf structure
       obj.check_if_surf_valid(surf_struct);
      
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
    function surface = get_surf(obj, surface_name)
      % surf = obj.get_surf(surface_name_string)
      %
      % Input: 
      %   surface_name: A string that matches the name of the 
      %   surface in the surf struct array in the object.
      %
      % Return: 
      %   A surf structure, if found in the surf 
      %   struct array of the object.
      
      if ~isa(surface_name, 'char')
        error('Input must have type char.');
      else
        index = find(strcmpi(surface_name, {obj.surf.name}));
        if index
          surface = obj.surf(index);
        else
          surface = [];
        end
      end    
    end
    
    %% set_name
    function obj = set_name(obj, surface_old_name, surface_new_name)
      % obj.set_name(surface_old_name, surface_new_name)
      %
      % Input: 
      %   surface_old_name: A string that matches the name of 
      %   the surface in the surf struct array in the object.
      %   surface_new_name: A string that you want to replace the
      %   surface_old_name.
      %
      % Return: 
      %   1. surface_old_name changed to surface_new_name, if 
      %   surface_old_name is found.
      %   2. Throws an error if surface_old_name is not found.
      
      if ~(isa(surface_old_name, 'char') && isa(surface_new_name, 'char'))
        error('Both input must have type char.');
      else
        match_idx = obj.get_index(surface_old_name);            
        obj.surf(match_idx).name = surface_new_name;
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
      obj.check_if_surf_valid(surf_struct);
            
      % Sets an existing surface to the new field values
      match_idx = obj.get_index(surf_struct.name);
      obj.surf(match_idx) = surf_struct;     
    end
    
    %% remove_surf
    function obj = remove_surf(obj, surface_name)
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
      
      match_idx = obj.get_index(surface_name);

      for surf_idx = 1:length(obj.surf)

        if obj.surf(surf_idx).surf_layer == match_idx
          obj.surf(surf_idx).surf_layer = [];
        elseif obj.surf(surf_idx).surf_layer > match_idx
            obj.surf(surf_idx).surf_layer = obj.surf(surf_idx).surf_layer -1;            
        end

        if obj.surf(surf_idx).active_layer == match_idx
          obj.surf(surf_idx).active_layer = [];
        elseif obj.surf(surf_idx).active_layer > match_idx
          obj.surf(surf_idx).active_layer = obj.surf(surf_idx).active_layer -1;  
        end

        if obj.surf(surf_idx).mask_layer == match_idx
          obj.surf(surf_idx).mask_layer = [];
        elseif obj.surf(surf_idx).mask_layer > match_idx
          obj.surf(surf_idx).mask_layer = obj.surf(surf_idx).mask_layer -1;  
        end

        if obj.surf(surf_idx).control_layer == match_idx
          obj.surf(surf_idx).control_layer = [];
        elseif obj.surf(surf_idx).control_layer > match_idx
          obj.surf(surf_idx).control_layer = obj.surf(surf_idx).control_layer -1;  
        end

        if obj.surf(surf_idx).quality_layer == match_idx
          obj.surf(surf_idx).quality_layer = [];
        elseif obj.surf(surf_idx).quality_layer > match_idx
          obj.surf(surf_idx).quality_layer = obj.surf(surf_idx).quality_layer -1;  
        end
      end

      obj.surf = [obj.surf(1:match_idx-1) obj.surf(match_idx+1:end)];
    end
    
    %% save_surf
    function [] = save_surf(obj, filename)
      % obj.save_surf(filename)
      %
      % Input: 
      %   filename: A string that will be the name 
      %   of the saved file.
      %
      % Result:
      %   Saves the surf struct array in the object 
      %   as a .mat file.
      
      surf = obj.surf;
      save(filename, 'surf', '-v7.3');
    end

    %% get_index
    function surface_index = get_index(obj, surface_name)
      % obj.get_index(surface_name_string)
      % 
      % Input: 
      %   surface_name_string: A string that matches 
      %   the name of the surface in the surf struct 
      %   array in the object.
      %
      % Return:
      %   1. The index of that surf struct in the surf 
      %   struct array in the object, if it is found.
      %   2. Throws an error if it is not found.
      
      for index = 1:length(obj.surf) 
          if strcmp(obj.surf(index).name, surface_name)
            surface_index = index;
            return
          end
      end
      error('%s is not found in the surf struct array.', surface_name);
    end
   
    
    %% check_if_surf_valid
    function [] = check_if_surf_valid(obj, surf_struct)
      % obj.check_if_surf_valid(surf_struct)
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
                          'name', 'surf_layer', 'active_layer', ...
                          'mask_layer', 'control_layer', ...
                          'quality_layer','visible'})
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
   
      if ~isa(surf_struct.surf_layer, 'double')
        error('Invalid field type for the field surf_layer (Should be type double).');
      end

      if ~isa(surf_struct.active_layer, 'double')
        error('Invalid field type for the field active_layer (Should be type double).');
      end
      
      if ~isa(surf_struct.mask_layer, 'double')
        error('Invalid field type for the field mask_layer (Should be type double).');
      end
      
      if ~isa(surf_struct.control_layer, 'double')
        error('Invalid field type for the field control_layer (Should be type double).');
      end
           
      if ~isa(surf_struct.quality_layer, 'double')
        error('Invalid field type for the field quality_layer (Should be type double).');
      end
      
      if ~isa(surf_struct.visible, 'logical')
        error('Invalid field type for the field visible (Should be type logical).');
      end
            
      % Check that size of x and y match x and y of other layers
      if ~isequal(size(surf_struct.x), size(surf_struct.y))
        error('Size of the field x and field y does not match.');
      end
      
      % check if all the surface_struct in the object has the same dimension in x and y
      if ~isempty(obj.surf) && ...
          ~all(arrayfun(@(arr) isequal(size(surf_struct.x), size(arr.x)), obj.surf))

          error('Field x has different size than all the other surf structure in the struct array.');
          % we only check x because x and y must have the same dimension in
          % the previous check.
      end
      
      % Check that layers indices point to valid layer indices:      
      if ~isempty(surf_struct.surf_layer) ...
          && ~(surf_struct.surf_layer > 0 || surf_struct.surf_layer < length(obj.surf))
        error('surf_layer index out of range.');
      end
      
      if ~isempty(surf_struct.active_layer) ...
          && ~(surf_struct.active_layer > 0 || surf_struct.active_layer < length(obj.surf))
        error('active_layer index out of range.');
      end
      
      if ~isempty(surf_struct.mask_layer) ...
          && ~(surf_struct.mask_layer > 0 || surf_struct.mask_layer < length(obj.surf))
        error('mask_layer index out of range.');
      end
      
      if ~isempty(surf_struct.control_layer) ...
          && ~(surf_struct.control_layer > 0 || surf_struct.control_layer < length(obj.surf))
        error('control_layer index out of range.');
      end
      
      if ~isempty(surf_struct.quality_layer) ...
          && ~(surf_struct.quality_layer > 0 || surf_struct.quality_layer < length(obj.surf))
        error('quality_layer index out of range.');
      end
    end
    
  end
end

