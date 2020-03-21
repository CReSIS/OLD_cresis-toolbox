classdef layerdata < handle
  % layerdata(param,layerdata_source) < handle
  %
  % A class that manages layer data for a segment.
  %
  % param: parameter structure from parameter spreadsheet. Specific fields
  % used:
  %  param.day_seg
  %  param.radar_name
  %  param.radar.lever_arm_fh
  %  param.records.gps.time_offset
  %  param.season_name
  %  param.sw_version
  % layerdata_source: ct_filename_out input argument to locate the filepath
  %   to the layer files and the layer organization file.
  %
  % Example of how to generate a filepath:
  %   param.radar_name   = 'rds';
  %   param.season_name  = '2019_Antarctica_Ground';
  %   param.day_seg = '20200107_01'
  %   fn = fullfile(ct_filename_out(param,surfdata_ref,''),...
  %     sprintf('Data_%s_%03d.mat',param.day_seg,frm));
  %
  % Author: John Paden
  %
  % See also: run_layerdata.m, layerdata.m
  
  properties (Constant)
  end
  
  properties
    param
    
    layer_organizer % layer_organizer contents
    
    layer % Nfrm element cell array
    % layer(frm).gps_time: 1 by Nx double vector
    % layer(frm).twtt Nlayers by Nx double vector
    
    layer_organizer_modified % flag to indicate layer_organizer modified since last save
    layer_modified % 1 by Nfrm logical vector flag to indicate corresponding layer modified since last save
    
    frames
    
    records
    along_track
  end
  
  methods
    %% constructor
    function obj = layerdata(param,layerdata_source)
      if ~exist('layerdata_source','var') || isempty(layerdata_source)
        layerdata_source = 'layer';
      end
      obj.param = param;
      obj.layerdata_source = layerdata_source;
      obj.frames = [];
      obj.records = [];
      obj.along_track = [];
      obj.layer_organizer = [];
      obj.layer = [];
      obj.layer_organizer_modified = false;
      obj.layer_modified = false();
    end
    
    function gps_time = gps_time(obj,frms)
      gps_time = [];
      for frm = 1:length(frms)
        obj.check(frm);
        gps_time(end+(1:length(obj.layer{frm}.gps_time))) = obj.layer{frm}.gps_time;
      end
    end
    
    function lon(obj,frms)
      lon = [];
      for frm = 1:length(frms)
        obj.check(frm);
        lon(end+(1:length(obj.layer{frm}.lon))) = obj.layer{frm}.lon;
      end
    end
    
    function lat(obj,frms)
      lat = [];
      for frm = 1:length(frms)
        obj.check(frm);
        lat(end+(1:length(obj.layer{frm}.lat))) = obj.layer{frm}.lat;
      end
    end
    
    function lon(obj,frms)
      lon = [];
      for frm = 1:length(frms)
        obj.check(frm);
        lon(end+(1:length(obj.layer{frm}.lon))) = obj.layer{frm}.lon;
      end
    end
    
    function [twtt,quality,type] = get_layer(obj,frms,id)
        obj.check_layer_organizer();
      if ischar(id)
        % name passed in rather than id
        id = find(strcmpi(id,obj.layer_organizer.lyr_name));
        if isempty(id)
          error('Layer does not exist in layer organizer. Run insert_layer() first.');
        end
      else
        % id passed in
        if all(obj.layer_organizer.lyr_id ~= id)
          error('Layer does not exist in layer organizer. Run insert_layer() first.');
        end
      end
      twtt = [];
      quality = [];
      type = [];
      for frm = 1:length(frms)
        obj.check(frm);
        lay_idx = find(obj.layer{frm}.id == id)
        Nx = length(obj.layer{frm}.gps_time);
        if isempty(lay_idx)
          twtt(end+(1:Nx) = obj.layer{frm}.twtt(lay_idx,:);
          quality(end+(1:Nx))) = obj.layer{frm}.quality(lay_idx,:);
          type(end+(1:Nx)) = obj.layer{frm}.type(lay_idx,:);
        else
          twtt(end+(1:Nx) = NaN;
          quality(end+(1:Nx))) = 1;
          type(end+(1:Nx)) = 2;
        end
      end
    end
    
    function insert_layers(obj,layer_organizer)
      obj.check_layer_organizer();
      
      try
        Nlayers = length(layer_organizer.lyr_name);
        obj.layer_organizer.lyr_age(end+(1:Nlayers)) = layer_organizer.lyr_age;
        obj.layer_organizer.lyr_age_source(end+(1:Nlayers)) = layer_organizer.lyr_age_source;
        obj.layer_organizer.lyr_desc(end+(1:Nlayers)) = layer_organizer.lyr_desc;
        obj.layer_organizer.lyr_group_name(end+(1:Nlayers)) = layer_organizer.lyr_group_name;
        obj.layer_organizer.lyr_id(end+(1:Nlayers)) = max(obj.layer_organizer.lyr_id) + (1:Nlayers);
        obj.layer_organizer.lyr_name(end+(1:Nlayers)) = layer_organizer.lyr_name;
        obj.layer_organizer.lyr_order(end+(1:Nlayers)) = layer_organizer.lyr_order;
        obj.fix_layer_organizer();
        obj.layer_organizer_modified = true;
      catch ME
        obj.fix_layer_organizer();
      end
    end
    
    function update_layer(obj,frms,id,gps_time,twtt,quality,type)
      obj.check_layer_organizer();
      if ischar(id)
        % name passed in rather than id
        id = find(strcmpi(id,obj.layer_organizer.lyr_name));
        if isempty(id)
          error('Layer does not exist in layer organizer. Run insert_layer() first.');
        end
      else
        if all(obj.layer_organizer.lyr_id ~= id)
          error('Layer does not exist in layer organizer. Run insert_layer() first.');
        end
      end
      
      for frm = 1:length(frms)
        obj.check(frm);
        lay_idx = find(obj.layer{frm}.id == id)
        if isempty(lay_idx)
          lay_idx = length(obj.layer{frm}.id) + 1;
          obj.layer{frm}.id = id;
        end
        obj.layer{frm}.twtt(lay_idx,:) = interp1(gps_time,twtt,obj.layer{frm}.gps_time);
        quality(end+(1:Nx))) = interp1(gps_time,twtt,obj.layer{frm}.quality,'nearest');
        type(end+(1:Nx)) = interp1(gps_time,twtt,obj.layer{frm}.type,'nearest');
        obj.layer_modified(frm) = true;
      end
    end
    
    function save(obj)
      if obj.layer_organizer_modified == true
        layer_organizer_fn = fullfile(ct_filename_out(obj.param,'',obj.layerdata_source),sprintf('layer_%s.mat',obj.param.day_seg));
        layer_organizer = obj.layer_organizer;
        ct_save(layer_organizer_fn,'-struct','layer_organizer');
        obj.layer_organizer_modified = false;
      end

      for frm = 1:length(obj.layer_modified)
        if obj.layer_modified(frm) == true
          layer_fn = fullfile(ct_filename_out(obj.param,'',obj.layerdata_source),sprintf('Data_%s_%03d.mat',obj.param.day_seg,frm))
          layer = obj.layer{frm};
          ct_save(layer_fn,'-struct','layer');
          obj.layer_modified(frm) = false;
        end
      end
    end
    
    function check(obj,frm)
      if length(obj.layers) < frm || isempty(obj.layers{frm})
        % Load layer file
        layer_fn = fullfile(ct_filename_out(obj.param,'',obj.layerdata_source),sprintf('Data_%s_%03d.mat',obj.param.day_seg,frm))
        obj.load(obj,frm,layer_fn);
      end
    end
    
    function load(obj,frm,layer_fn)
      if ~exist('layer_fn','var') || isempty(layer_fn)
        layer_fn = fullfile(ct_filename_out(obj.param,'',obj.layerdata_source),sprintf('Data_%s_%03d.mat',obj.param.day_seg,frm))
      end
      if ~exist(layer_fn,'file')
        obj.create(obj,frm);
      else
        layer = load(layer_fn);
        obj.layer_modified(frm) = false;
        obj.layer{frm} = obj.fix(frm,layer);
      end
    end
    
    function layer = fix(obj,frm,layer)
      
      obj.check_layer_organizer();
      
      % Verifies all fields and fixes any problems for layer file struct
      try
        if ~isfield(layer,'gps_time')
          error('gps_time does not exist');
        end
        if ~isfield(layer,'elev')
          error('elev does not exist');
        end
        if ~isfield(layer,'lat')
          error('lat does not exist');
        end
        if ~isfield(layer,'lon')
          error('lon does not exist');
        end
        if ~isa(layer.gps_time,'double')
          error('gps_time is not double type');
        end
        if ~isa(layer.elev,'double')
          error('elev is not double type');
        end
        if ~isa(layer.lat,'double')
          error('lat is not double type');
        end
        if ~isa(layer.lon,'double')
          error('lon is not double type');
        end
        Nx = length(layer.gps_time);
        if length(layer.elev) ~= Nx
          error('elev field length does not match gps_time length');
        end
        if length(layer.lat) ~= Nx
          error('lat field length does not match gps_time length');
        end
        if length(layer.lon) ~= Nx
          error('lon field length does not match gps_time length');
        end
      catch ME
        update_gps(obj,frm);
      end

      % Check existence for each field
      if ~isfield(layer,'twtt')
        layer.twtt = zeros(1,Nx,'double');
        obj.layer_modified(frm) = true;
      end
      if ~isfield(layer,'quality')
        layer.quality = zeros(1,Nx,'double');
        obj.layer_modified(frm) = true;
      end
      if ~isfield(layer,'type')
        layer.type = zeros(1,Nx,'double');
        obj.layer_modified(frm) = true;
      end
      
      % Check that type is correct for each field
      if isa(layer.id,'double')
        layer.id = double(layer.id);
         obj.layer_modified(frm) = true;
     end
      if isa(layer.twtt,'double')
        layer.quality = double(layer.quality);
        obj.layer_modified(frm) = true;
      end
      if isa(layer.quality,'uint8')
        layer.quality = uint8(layer.quality);
        obj.layer_modified(frm) = true;
      end
      if isa(layer.type,'uint8')
        layer.type = uint8(layer.type);
        obj.layer_modified(frm) = true;
      end
      
      % Check that number of layers is correct for each field
      if size(layer.twtt,1) < length(layer.id)
        layer.id = layer.id(1:size(layer.twtt,1));
        obj.layer_modified(frm) = true;
      end
      if size(layer.twtt,1) > length(layer.id)
        Nlyrs_new = size(layer.twtt,1) - length(obj.layer_organizer.id);
        layer.id(end+(1:Nlyrs_new)) = max(layer.id) + (1:Nlyrs_new);
        obj.layer_modified(frm) = true;
      end
      Nlyrs = length(layer.id);
      if size(layer.quality,1) < Nlyrs
        layer.quality(end+1:Nlyrs) = ones(Nlyrs,size(layer.quality,2),'uint8');
        obj.layer_modified(frm) = true;
      end
      if size(layer.quality,1) > Nlyrs
        layer.quality = layer.quality(1:length(layer.id));
        obj.layer_modified(frm) = true;
      end
      if size(layer.twtt,1) < Nlyrs
        layer.twtt(end+1:Nlyrs) = 2*ones(Nlyrs,size(layer.twtt,2),'uint8');
        obj.layer_modified(frm) = true;
      end
      if size(layer.twtt,1) > Nlyrs
        layer.twtt = layer.twtt(1:length(layer.id));
        obj.layer_modified(frm) = true;
      end
      
      % Ensure that number of range lines is correct for each field
      if size(layer.twtt,2) < Nx
        Nx_new = Nx - size(layer.twtt,2);
        layer.twtt(:,end+(1:Nx_new)) = nan(1,Nx_new);
        obj.layer_modified(frm) = true;
      end
      if size(layer.type,2) < Nx
        Nx_new = Nx - size(layer.type,2);
        layer.type(:,end+(1:Nx_new)) = 2*ones(1,Nx_new,'uint8');
        obj.layer_modified(frm) = true;
      end
      if size(layer.quality,2) < Nx
        Nx_new = Nx - size(layer.quality,2);
        layer.quality(:,end+(1:Nx_new)) = ones(1,Nx_new,'uint8');
      end
      if size(layer.twtt,2) > Nx
        layer.twtt(:,Nx+1:end) = [];
        obj.layer_modified(frm) = true;
      end
      if size(layer.type,2) > Nx
        layer.type(:,Nx+1:end) = [];
        obj.layer_modified(frm) = true;
      end
      if size(layer.quality,2) > Nx
        layer.quality(:,Nx+1:end) = [];
        obj.layer_modified(frm) = true;
      end

      % Ensure all values are valid
      mask = ~isfinite(layer.quality ~= 1 & layer.quality ~= 2 & layer.quality ~= 3);
      if any(mask)
        layer.quality(mask) = 1;
        obj.layer_modified(frm) = true;
      end
      mask = ~isfinite(layer.type ~= 1 & layer.type ~= 2 & layer.type ~= 3 & layer.type ~= 4);
      if any(mask)
        layer.quality(mask) = 2;
        obj.layer_modified(frm) = true;
      end
    end
    
    function check_layer_organizer(obj)
      if isempty(layer_organizer)
        load_layer_organizer();
      end
    end
    
    function load_layer_organizer(obj)
      layer_organizer_fn = fullfile(ct_filename_out(obj.param,'',obj.layerdata_source),sprintf('layer_%s.mat',obj.param.day_seg));
      if ~exist(layer_organizer_fn,'file')
        obj.create_layer_organizer(obj,frm);
      else
        obj.layer_organizer = load(layer_organizer_fn);
        obj.layer_organizer_modified = false;
        obj.fix_layer_organizer();
      end
    end
    
    function create_layer_organizer(obj,frm)
      obj.layer_organizer.file_version = '1';
      obj.layer_organizer.file_type = 'layer_organizer';
      obj.layer_organizer.param.radar_name = obj.param.radar_name;
      obj.layer_organizer.param.season_name = obj.param.season_name;
      obj.layer_organizer.param.day_seg = obj.param.day_seg;
      obj.layer_organizer.param.sw_version = obj.param.sw_version;
      obj.layer_organizer.lyr_age = [];
      obj.layer_organizer.lyr_age_source = [];
      obj.layer_organizer.lyr_desc = [];
      obj.layer_organizer.lyr_group_name = [];
      obj.layer_organizer.lyr_id = [];
      obj.layer_organizer.lyr_name = [];
      obj.layer_organizer.lyr_order = [];
      obj.layer_organizer_modified(frm) = true;
    end

    function fix_layer_organizer(obj,layer_organizer)
      % Check the type
      if ~isdouble(obj.layer_organizer.lyr_age)
        obj.layer_organizer.lyr_age = double(obj.layer_organizer.lyr_age);
        obj.layer_organizer_modified(frm) = true;
      end
      if ~iscell(obj.layer_organizer.lyr_age_source)
        obj.layer_organizer.lyr_age_source = {};
        obj.layer_organizer_modified(frm) = true;
      else
        for cell_idx = 1:length(obj.layer_organizer.lyr_age_source)
          if ~isstruct(obj.layer_organizer.lyr_age_source{cell_idx})
            obj.layer_organizer.lyr_age_source{cell_idx} = struct('age',{},'source',{},'type',{});
            obj.layer_organizer_modified(frm) = true;
          else
            for struct_idx = 1:length(obj.layer_organizer.lyr_age_source{cell_idx})
              if ~isdouble(obj.layer_organizer.lyr_age_source{cell_idx}(struct_idx).age)
                obj.layer_organizer.lyr_age_source{cell_idx}(struct_idx).age = double(obj.layer_organizer.lyr_age_source{cell_idx}(struct_idx).age);
                obj.layer_organizer_modified(frm) = true;
              end
              if ~ischar(obj.layer_organizer.lyr_age_source{cell_idx}(struct_idx).source)
                obj.layer_organizer.lyr_age_source{cell_idx}(struct_idx).source = char(obj.layer_organizer.lyr_age_source{cell_idx}(struct_idx).source);
                obj.layer_organizer_modified(frm) = true;
              end
              if ~isdouble(obj.layer_organizer.lyr_age_source{cell_idx}(struct_idx).type)
                obj.layer_organizer.lyr_age_source{cell_idx}(struct_idx).type = double(obj.layer_organizer.lyr_age_source{cell_idx}(struct_idx).type);
                obj.layer_organizer_modified(frm) = true;
              end
              if ~any(obj.layer_organizer.lyr_age_source{cell_idx}(struct_idx).type == [0 1 2 3 4])
                obj.layer_organizer.lyr_age_source{cell_idx}(struct_idx).type = 0;
                obj.layer_organizer_modified(frm) = true;
              end
            end
          end
        end
      end
      if ~iscell(obj.layer_organizer.lyr_desc)
        obj.layer_organizer.lyr_desc = {};
        obj.layer_organizer_modified(frm) = true;
      else
        for cell_idx = 1:length(obj.layer_organizer.lyr_desc)
          if ~ischar(obj.layer_organizer.lyr_desc{cell_idx})
            obj.layer_organizer.lyr_desc{cell_idx} = char(obj.layer_organizer.lyr_desc{cell_idx});
            obj.layer_organizer_modified(frm) = true;
          end
        end
      end
      if ~iscell(obj.layer_organizer.lyr_group_name)
        obj.layer_organizer.lyr_group_name = {};
        obj.layer_organizer_modified(frm) = true;
      else
        for cell_idx = 1:length(obj.layer_organizer.lyr_group_name)
          if ~ischar(obj.layer_organizer.lyr_group_name{cell_idx})
            obj.layer_organizer.lyr_group_name{cell_idx} = char(obj.layer_organizer.lyr_group_name{cell_idx});
            obj.layer_organizer_modified(frm) = true;
          end
        end
      end
      if ~isdouble(obj.layer_organizer.lyr_id)
        obj.layer_organizer.lyr_id = double(obj.layer_organizer.lyr_id);
        obj.layer_organizer_modified(frm) = true;
      end
      if ~iscell(obj.layer_organizer.lyr_name)
        obj.layer_organizer.lyr_name = {};
        obj.layer_organizer_modified(frm) = true;
      else
        for cell_idx = 1:length(obj.layer_organizer.lyr_name)
          if ~ischar(obj.layer_organizer.lyr_name{cell_idx})
            obj.layer_organizer.lyr_name{cell_idx} = char(obj.layer_organizer.lyr_name{cell_idx});
            obj.layer_organizer_modified(frm) = true;
          end
        end
      end
      if ~isdouble(obj.layer_organizer.lyr_order)
        obj.layer_organizer.lyr_order = double(obj.layer_organizer.lyr_order);
        obj.layer_organizer_modified(frm) = true;
      end
      
      % Check the length (too long)
      Nlayers = length(obj.layer_organizer.lyr_id);
      if length(obj.layer_organizer.lyr_age) > Nlayers
        obj.layer_organizer.lyr_age(Nlayers+1:end) = [];
        obj.layer_organizer_modified(frm) = true;
      end
      if length(obj.layer_organizer.lyr_age_source) > Nlayers
        obj.layer_organizer.lyr_age_source(Nlayers+1:end) = [];
        obj.layer_organizer_modified(frm) = true;
      end
      if length(obj.layer_organizer.lyr_desc) > Nlayers
        obj.layer_organizer.lyr_desc(Nlayers+1:end) = [];
        obj.layer_organizer_modified(frm) = true;
      end
      if length(obj.layer_organizer.lyr_group_name) > Nlayers
        obj.layer_organizer.lyr_group_name(Nlayers+1:end) = [];
        obj.layer_organizer_modified(frm) = true;
      end
      if length(obj.layer_organizer.lyr_name) > Nlayers
        obj.layer_organizer.lyr_name(Nlayers+1:end) = [];
        obj.layer_organizer_modified(frm) = true;
      end
      if length(obj.layer_organizer.lyr_order) > Nlayers
        obj.layer_organizer.lyr_order(Nlayers+1:end) = [];
        obj.layer_organizer_modified(frm) = true;
      end
      
      % Check the length (too short)
      if length(obj.layer_organizer.lyr_age) < Nlayers
        obj.layer_organizer.lyr_age(end+1:Nlayers) = NaN;
        obj.layer_organizer_modified(frm) = true;
      end
      if length(obj.layer_organizer.lyr_age_source) < Nlayers
        obj.layer_organizer.lyr_age_source(end+1:Nlayers) = cellfun(@struct,cell(1,Nlayers-length(obj.layer_organizer.lyr_age_source)),'UniformOutput',false);
        obj.layer_organizer_modified(frm) = true;
      end
      if length(obj.layer_organizer.lyr_desc) < Nlayers
        obj.layer_organizer.lyr_desc(end+1:Nlayers) = cellfun(@char,cell(1,Nlayers-length(obj.layer_organizer.lyr_desc)),'UniformOutput',false);
        obj.layer_organizer_modified(frm) = true;
      end
      if length(obj.layer_organizer.lyr_group_name) < Nlayers
        obj.layer_organizer.lyr_group_name(end+1:Nlayers) = cellfun(@char,cell(1,Nlayers-length(obj.layer_organizer.lyr_group_name)),'UniformOutput',false);
        obj.layer_organizer_modified(frm) = true;
      end
      if length(obj.layer_organizer.lyr_name) < Nlayers
        obj.layer_organizer.lyr_name(end+1:Nlayers) = cellfun(@char,cell(1,Nlayers-length(obj.layer_organizer.lyr_name)),'UniformOutput',false);
        obj.layer_organizer_modified(frm) = true;
      end
      if length(obj.layer_organizer.lyr_order) < Nlayers
        obj.layer_organizer.lyr_order(end+1:Nlayers) = max(obj.layer_organizer.lyr_order) + (1:Nlayers);
        obj.layer_organizer_modified(frm) = true;
      end
      
      % Check id uniqueness
      [sorted_ids unique_idxs] = unique(obj.layer_organizer.lyr_id);
      if length(sorted_ids) < Nlayers
        obj.layer_organizer.lyr_age = obj.layer_organizer.lyr_age(unique_idxs);
        obj.layer_organizer.lyr_age_source = obj.layer_organizer.lyr_age_source(unique_idxs);
        obj.layer_organizer.lyr_desc = obj.layer_organizer.lyr_desc(unique_idxs);
        obj.layer_organizer.lyr_group_name = obj.layer_organizer.lyr_group_name(unique_idxs);
        obj.layer_organizer.lyr_id = sorted_ids;
        obj.layer_organizer.lyr_name = obj.layer_organizer.lyr_name(unique_idxs);
        obj.layer_organizer.lyr_order = obj.layer_organizer.lyr_order(unique_idxs);
        Nlayers = length(obj.layer_organizer.lyr_id);
        obj.layer_organizer_modified(frm) = true;
      end
      
      % Check name uniqueness
      for idx = 2:Nlayers
        base_name = obj.layer_organizer.lyr_name{idx};
        test_name = obj.layer_organizer.lyr_name{idx};
        test_idx = 0;
        while any(strcmp(test_name,obj.layer_organizer.lyr_name(1:idx-1)))
          test_idx = test_idx + 1;
          test_name = sprintf('%s_%03d', base_name, test_idx);
        end
        if ~strcmp(test_name, obj.layer_organizer.lyr_name{idx})
          obj.layer_organizer.lyr_name{idx} = test_name;
          obj.layer_organizer_modified(frm) = true;
        end
      end
      
    end
    
    function update_gps(obj,frm)
      obj.check_records();

      if frm < length(obj.frame_idxs)
        frm_mask = find(obj.along_track >= obj.records.along_track(obj.frames.frame_idxs(frm)) ...
          & obj.along_track < obj.records.along_track(obj.frames.frame_idxs(frm+1)))
      else
        frm_mask = find(obj.along_track >= obj.records.along_track(obj.frames.frame_idxs(frm)));
      end
      obj.layer{frm}.gps_time = interp1(obj.records.along_track,obj.records.gps_time,obj.along_track(frm_mask));
      obj.layer{frm}.elev = interp1(obj.records.along_track,obj.records.elev,obj.along_track(frm_mask));
      obj.layer{frm}.lat = interp1(obj.records.along_track,obj.records.lat,obj.along_track(frm_mask));
      obj.layer{frm}.lon = interp1(obj.records.along_track,obj.records.lon,obj.along_track(frm_mask));
      obj.layer_modified(frm) = true;
    end
    
    function create(obj,frm)
      obj.check_records();

      if frm < length(obj.frame_idxs)
        frm_mask = find(obj.along_track >= obj.records.along_track(obj.frames.frame_idxs(frm)) ...
          & obj.along_track < obj.records.along_track(obj.frames.frame_idxs(frm+1)))
      else
        frm_mask = find(obj.along_track >= obj.records.along_track(obj.frames.frame_idxs(frm)));
      end
      obj.layer{frm}.file_version = '1';
      obj.layer{frm}.file_type = 'layer';
      obj.layer{frm}.param.radar_name = obj.param.radar_name;
      obj.layer{frm}.param.season_name = obj.param.season_name;
      obj.layer{frm}.param.day_seg = obj.param.day_seg;
      obj.layer{frm}.param.records.gps.time_offset = obj.param.records.gps.time_offset;
      obj.layer{frm}.param.radar.lever_arm_fh = obj.param.radar.lever_arm_fh;
      obj.layer{frm}.param.sw_version = obj.param.sw_version;
      obj.layer{frm}.gps_time = interp1(obj.records.along_track,obj.records.gps_time,obj.along_track(frm_mask));
      obj.layer{frm}.elev = interp1(obj.records.along_track,obj.records.elev,obj.along_track(frm_mask));
      obj.layer{frm}.lat = interp1(obj.records.along_track,obj.records.lat,obj.along_track(frm_mask));
      obj.layer{frm}.lon = interp1(obj.records.along_track,obj.records.lon,obj.along_track(frm_mask));
      obj.layer{frm}.id = [];
      Nx = length(obj.layer{frm}.gps_time);
      obj.layer{frm}.twtt = zeros(0,size(1,Nx));
      obj.layer{frm}.quality = zeros(0,size(1,Nx),'uint8');
      obj.layer{frm}.type = zeros(0,size(1,Nx),'uint8');
      obj.layer_modified(frm) = true;
    end
    
    function check_records(obj,frm,layer_fn)
      if isempty(obj.records)
        obj.load_records();
      end
    end
    
    function load_records(obj,records)
      
      if ~exist('records','var') || isempty(records)
        % Load records file
        records_fn = ct_filename_support(param,'','records');
        obj.records = load(records_fn);
      end
      
      % Load frames file
      obj.frames = load(ct_filename_support(param,'','frames'));
      if ~isfield(obj.frames,'frame_idxs')
        warning('Old frames file format. frames_update.m should be run on this segment %s.', obj.param.day_seg);
        % Convert loaded frames file to new format
        obj.frames = obj.frames.frames;
        obj.frames.gps_time = [obj.records.gps_time(obj.frames.frame_idxs), obj.records.gps_time(end)];
      end
      
      % Create reference trajectory (rx_path == 0, tx_weights = []). Update
      % the records field with this information.
      trajectory_param = struct('gps_source',obj.records.gps_source, ...
        'season_name',obj.param.season_name,'radar_name',obj.param.radar_name,'rx_path', 0, ...
        'tx_weights', [], 'lever_arm_fh', obj.param.radar.lever_arm_fh);
      obj.records = trajectory_with_leverarm(obj.records,trajectory_param);
      
      obj.records.along_track = geodetic_to_along_track(obj.records.lat,obj.records.lon);
      
      sys = ct_output_dir(obj.param.radar_name);
      switch (sys)
        case {'rds','accum'}
          sample_spacing = 15;
        otherwise
          sample_spacing = 5;
      end
      
      obj.along_track = [0:sample_spacing:obj.records.along_track(end)]
    end
    
  end
  
end




