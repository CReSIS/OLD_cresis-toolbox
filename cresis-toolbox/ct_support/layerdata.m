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
  % Author: John Paden
  %
  % See also: run_layerdata.m, layerdata.m
  
  properties (Constant)
  end
  
  properties
    param
    layerdata_source
    
    layer_organizer % layer_organizer contents
    
    layer % Nfrm element cell array
    % layer(frm).gps_time: 1 by Nx double vector
    % layer(frm).twtt Nlayers by Nx double vector
    
    layer_organizer_modified % flag to indicate layer_organizer modified since last save
    layer_modified % 1 by Nfrm logical vector flag to indicate corresponding layer modified since last save
    
    frames
    
    records
    
    along_track
    sample_spacing_method % 'records' or 'along_track' (default)
    along_track_spacing % either 5 m (non-rds) or 15 m (rds) by default
    record_spacing % 200 by default

  end
  
  %% Private Methods
  methods (Access = private)
    
    %% create: create layer file
    function create(obj,frm)
      obj.check_records();

      obj.layer{frm}.file_version = '1';
      obj.layer{frm}.file_type = 'layer';
      obj.layer{frm}.param.radar_name = obj.param.radar_name;
      obj.layer{frm}.param.season_name = obj.param.season_name;
      obj.layer{frm}.param.day_seg = obj.param.day_seg;
      obj.layer{frm}.param.records.gps.time_offset = obj.param.records.gps.time_offset;
      obj.layer{frm}.param.radar.lever_arm_fh = obj.param.radar.lever_arm_fh;
      obj.layer{frm}.param.sw_version = obj.param.sw_version;
      
      if strcmp(obj.sample_spacing_method,'records')
        if frm < length(obj.frames.frame_idxs)
          recs = obj.frames.frame_idxs(frm):obj.record_spacing:obj.frames.frame_idxs(frm+1);
        else
          recs = obj.frames.frame_idxs(frm):obj.record_spacing:obj.frames.Nx;
        end
        obj.layer{frm}.gps_time = obj.records.gps_time(recs);
        obj.layer{frm}.elev = obj.records.elev(recs);
        obj.layer{frm}.lat = obj.records.lat(recs);
        obj.layer{frm}.lon = obj.records.lon(recs);
      else
        if frm < length(obj.frames.frame_idxs)
          frm_mask = find(obj.along_track >= obj.records.along_track(obj.frames.frame_idxs(frm)) ...
            & obj.along_track < obj.records.along_track(obj.frames.frame_idxs(frm+1)));
        else
          frm_mask = find(obj.along_track >= obj.records.along_track(obj.frames.frame_idxs(frm)));
        end
        good_idxs = [true diff(obj.records.along_track)>0];
        obj.layer{frm}.gps_time = interp1(obj.records.along_track(good_idxs),obj.records.gps_time(good_idxs),obj.along_track(frm_mask));
        obj.layer{frm}.elev = interp1(obj.records.along_track(good_idxs),obj.records.elev(good_idxs),obj.along_track(frm_mask));
        obj.layer{frm}.lat = interp1(obj.records.along_track(good_idxs),obj.records.lat(good_idxs),obj.along_track(frm_mask));
        obj.layer{frm}.lon = interp1(obj.records.along_track(good_idxs),obj.records.lon(good_idxs),obj.along_track(frm_mask));
      end
      
      obj.layer{frm}.id = [];
      Nx = length(obj.layer{frm}.gps_time);
      obj.layer{frm}.twtt = zeros(0,size(1,Nx));
      obj.layer{frm}.quality = zeros(0,size(1,Nx),'uint8');
      obj.layer{frm}.type = zeros(0,size(1,Nx),'uint8');
      obj.layer_modified(frm) = true;
    end
    
    %% create_layer_organizer: create layer organizer if it does not exist
    function create_layer_organizer(obj,frm)
      obj.layer_organizer.file_version = '1';
      obj.layer_organizer.file_type = 'layer_organizer';
      obj.layer_organizer.param.radar_name = obj.param.radar_name;
      obj.layer_organizer.param.season_name = obj.param.season_name;
      obj.layer_organizer.param.day_seg = obj.param.day_seg;
      obj.layer_organizer.param.sw_version = obj.param.sw_version;
      obj.layer_organizer.lyr_age = [];
      obj.layer_organizer.lyr_age_source = {};
      obj.layer_organizer.lyr_desc = {};
      obj.layer_organizer.lyr_group_name = {};
      obj.layer_organizer.lyr_id = [];
      obj.layer_organizer.lyr_name = {};
      obj.layer_organizer.lyr_order = [];
      obj.layer_organizer_modified = true;
    end

  end
  
  %% Public Methods
  methods
    %% constructor
    function obj = layerdata(param,layerdata_source)
      if ~exist('layerdata_source','var') || isempty(layerdata_source)
        layerdata_source = 'layer';
      end
      obj.param = param;
      obj.param.radar_name = ct_output_dir(obj.param.radar_name);
      obj.layerdata_source = layerdata_source;
      obj.frames = [];
      obj.records = [];
      obj.layer_organizer = [];
      obj.layer = [];
      obj.layer_organizer_modified = false;
      obj.layer_modified = false();
      
      obj.along_track = [];
      obj.sample_spacing_method = '';
      obj.along_track_spacing = [];
      obj.record_spacing = 200;
    end
    
    %% check: check that layer is loaded
    function check(obj,frm)
      if length(obj.layer) < frm || isempty(obj.layer{frm})
        % Load layer file
        obj.load(frm);
      end
    end
    
    %% check_all: check to make sure layer organizer and all layers loaded
    function check_all(obj)
      obj.check_records();
      obj.check_layer_organizer();
      for frm = 1:length(obj.frames.frame_idxs)
        obj.check(frm);
      end
    end
    
    %% check_layer_organizer: check that layer organizer is loaded
    function check_layer_organizer(obj)
      if isempty(obj.layer_organizer)
        obj.load_layer_organizer();
      end
    end
    
    %% check_records: check that records and frames are loaded
    % records: optionally records may be passed in from external source
    function check_records(obj,records)
      if exist('records','var')
        obj.records = records;
      elseif isempty(obj.records)
        obj.load_records();
      end
    end
    
    %% delete: delete layer file
    function delete_layer_file(obj,frms)
      for frm = frms(:).'
        layer_fn = fullfile(ct_filename_out(obj.param,'',obj.layerdata_source),sprintf('Data_%s_%03d.mat',obj.param.day_seg,frm));
        try
          delete(layer_fn);
        end
      end
    end
    
    %% delete_all: delete all layers and layer organizer file
    function delete_all(obj)
      obj.check_records();
      obj.delete_layer_file(1:length(obj.frames.frame_idxs));
      obj.delete_layer_organizer();
    end
    
    %% delete_layer_organizer: delete layer organizer file
    function delete_layer_organizer(obj)
      layer_organizer_fn = fullfile(ct_filename_out(obj.param,'',obj.layerdata_source),sprintf('layer_%s.mat',obj.param.day_seg));
      try
        delete(layer_organizer_fn);
      end
    end
    
    %% elev: get elev
    function elev = elev(obj,frms)
      elev = [];
      for frm = 1:frms(:).'
        obj.check(frm);
        elev(end+(1:length(obj.layer{frm}.elev))) = obj.layer{frm}.elev;
      end
    end
    
    %% fix: fix any problems with a layer
    function fix(obj,frm)
      obj.check_layer_organizer();
    
      if ~isfield(obj.layer{frm},'file_version')
        % Old file format
        warning('Old layer file format frame %d\n', frm);
        lay = obj.layer{frm};
        obj.check_records();
        obj.layer_modified(frm) = true;
        
        % Remove data that is not contained within frame boundaries
        frms_mask = false(size(lay.GPS_time));
        if frm < length(obj.frames.frame_idxs)
          frms_mask(lay.GPS_time >= obj.frames.gps_time(frm)...
            & lay.GPS_time < obj.frames.gps_time(frm+1)) = true;
        else
          frms_mask(lay.GPS_time >= obj.frames.gps_time(frm)...
            & lay.GPS_time <= obj.frames.gps_time(frm+1)) = true;
        end
        Nx_old = length(lay.GPS_time);
        
        % Recreate lat, lon, elev data using records file
        lay.GPS_time = lay.GPS_time(frms_mask);
        lay.Elevation = interp1(obj.records.gps_time,obj.records.elev,lay.GPS_time);
        lay.Latitude = interp1(obj.records.gps_time,obj.records.lat,lay.GPS_time);
        lay.Longitude = interp1(obj.records.gps_time,obj.records.lon,lay.GPS_time);
        
        Nx = length(lay.GPS_time);
        
        new_layer_ids = [];
        for lay_idx = 1:length(lay.layerData)
          if ~isfield(lay.layerData{lay_idx},'value')
            lay.layerData{lay_idx}.value = {};
          end
          if isempty(lay.layerData{lay_idx}.value)
            lay.layerData{lay_idx}.value{1}.data = nan(1,Nx);
          end
          if length(lay.layerData{lay_idx}.value) < 2
            lay.layerData{lay_idx}.value{2}.data = nan(1,Nx);
          end
          if ~isfield(lay.layerData{lay_idx},'quality')
            lay.layerData{lay_idx}.quality = ones(1,Nx);
          end
          if isfield(lay.layerData{lay_idx},'name') && ~isfield(lay.layerData{lay_idx},'id')
            % Old file format, switch name to id
            match_idx = find(strcmp(lay.layerData{lay_idx}.name,obj.layer_organizer.lyr_name));
            
            if isempty(match_idx)
              id = obj.insert_layers(struct('lyr_name',{{lay.layerData{lay_idx}.name}}));
              lay.layerData{lay_idx}.id = id;
            else
              lay.layerData{lay_idx}.id = obj.layer_organizer.lyr_id(match_idx);
            end
            
          else
            if ~isfield(lay.layerData{lay_idx},'id')
              lay.layerData{lay_idx}.id = lay_idx;
            end
            match_idx = find(lay.layerData{lay_idx}.id == obj.layer_organizer.lyr_id);
            if isempty(match_idx)
              % Add the layer to the layer_organizer
              layer_organizer = [];
              if lay_idx == 1
                % Enforce standard name/group for old file format
                layer_organizer.lyr_name{1} = 'surface';
                layer_organizer.lyr_group_name{1} = 'standard';
              elseif lay_idx == 2
                % Enforce standard name/group for old file format
                layer_organizer.lyr_name{1} = 'bottom';
                layer_organizer.lyr_group_name{1} = 'standard';
              else
                layer_organizer.lyr_name{1} = '';
                layer_organizer.lyr_group_name{1} = sprintf('auto_%03d',lay_idx);
              end
              id = obj.insert_layers(layer_organizer);
              lay.layerData{lay_idx}.id = id;
            else
              lay.layerData{lay_idx}.id = obj.layer_organizer.lyr_id(match_idx);
            end
          end
          
          % Ensure all fields are consistent in length
          % -------------------------------------------------------------------
          % Too short:
          if length(lay.layerData{lay_idx}.quality) < Nx_old
            lay.layerData{lay_idx}.quality(end+1:Nx_old) = NaN;
          end
          if length(lay.layerData{lay_idx}.value{1}.data) < Nx_old
            lay.layerData{lay_idx}.value{1}.data(end+1:Nx_old) = NaN;
          end
          if length(lay.layerData{lay_idx}.value{2}.data) < Nx_old
            lay.layerData{lay_idx}.value{2}.data(end+1:Nx_old) = NaN;
          end
          % Too long:
          if length(lay.layerData{lay_idx}.quality) > Nx_old
            lay.layerData{lay_idx}.quality = lay.layerData{lay_idx}.quality(1:Nx_old);
          end
          if length(lay.layerData{lay_idx}.value{1}.data) > Nx_old
            lay.layerData{lay_idx}.value{1}.data = lay.layerData{lay_idx}.value{1}.data(1:Nx_old);
          end
          if length(lay.layerData{lay_idx}.value{2}.data) > Nx_old
            lay.layerData{lay_idx}.value{2}.data = lay.layerData{lay_idx}.value{2}.data(1:Nx_old);
          end
          lay.layerData{lay_idx}.quality = lay.layerData{lay_idx}.quality(frms_mask);
          lay.layerData{lay_idx}.value{1}.data = lay.layerData{lay_idx}.value{1}.data(frms_mask);
          lay.layerData{lay_idx}.value{2}.data = lay.layerData{lay_idx}.value{2}.data(frms_mask);
        end
        
        new_lay.elev = lay.Elevation;
        new_lay.file_type = 'layer';
        new_lay.file_version = '1';
        new_lay.gps_time = lay.GPS_time;
        new_lay.lat = lay.Latitude;
        new_lay.lon = lay.Longitude;
        new_lay.quality = zeros(length(lay.layerData),Nx,'uint8');
        new_lay.twtt = zeros(length(lay.layerData),Nx);
        new_lay.type = zeros(length(lay.layerData),Nx,'uint8');
        for lay_idx = 1:length(lay.layerData)
          new_lay.id(lay_idx,1) = lay.layerData{lay_idx}.id;
          new_lay.quality(lay_idx,1:Nx) = lay.layerData{lay_idx}.quality;
          new_lay.twtt(lay_idx,1:Nx) = lay.layerData{lay_idx}.value{2}.data;
          new_lay.type(lay_idx,1:Nx) = 1 + ~isfinite(lay.layerData{lay_idx}.value{1}.data);
        end
        obj.layer{frm} = new_lay;
      end
      
      % Verifies all fields and fixes any problems for layer file struct
      try
        if ~isfield(obj.layer{frm},'gps_time')
          error('gps_time does not exist');
        end
        if ~isfield(obj.layer{frm},'elev')
          error('elev does not exist');
        end
        if ~isfield(obj.layer{frm},'lat')
          error('lat does not exist');
        end
        if ~isfield(obj.layer{frm},'lon')
          error('lon does not exist');
        end
        if ~isa(obj.layer{frm}.gps_time,'double')
          error('gps_time is not double type');
        end
        if ~isa(obj.layer{frm}.elev,'double')
          error('elev is not double type');
        end
        if ~isa(obj.layer{frm}.lat,'double')
          error('lat is not double type');
        end
        if ~isa(obj.layer{frm}.lon,'double')
          error('lon is not double type');
        end
        Nx = length(obj.layer{frm}.gps_time);
        if length(obj.layer{frm}.elev) ~= Nx
          error('elev field length does not match gps_time length');
        end
        if length(obj.layer{frm}.lat) ~= Nx
          error('lat field length does not match gps_time length');
        end
        if length(obj.layer{frm}.lon) ~= Nx
          error('lon field length does not match gps_time length');
        end
      catch ME
        update_gps(obj,frm);
      end

      % Check existence for each field
      if ~isfield(obj.layer{frm},'twtt')
        obj.layer{frm}.twtt = zeros(1,Nx,'double');
        obj.layer_modified(frm) = true;
      end
      if ~isfield(obj.layer{frm},'quality')
        obj.layer{frm}.quality = zeros(1,Nx,'double');
        obj.layer_modified(frm) = true;
      end
      if ~isfield(obj.layer{frm},'type')
        obj.layer{frm}.type = zeros(1,Nx,'double');
        obj.layer_modified(frm) = true;
      end
      
      % Check that type is correct for each field
      if ~isa(obj.layer{frm}.id,'double')
        obj.layer{frm}.id = double(obj.layer{frm}.id);
        obj.layer_modified(frm) = true;
      end
      if ~isa(obj.layer{frm}.twtt,'double')
        obj.layer{frm}.quality = double(obj.layer{frm}.quality);
        obj.layer_modified(frm) = true;
      end
      if ~isa(obj.layer{frm}.quality,'uint8')
        obj.layer{frm}.quality = uint8(obj.layer{frm}.quality);
        obj.layer_modified(frm) = true;
      end
      if ~isa(obj.layer{frm}.type,'uint8')
        obj.layer{frm}.type = uint8(obj.layer{frm}.type);
        obj.layer_modified(frm) = true;
      end
      
      % Check that number of layers is correct for each field
      if size(obj.layer{frm}.twtt,1) < length(obj.layer{frm}.id)
        obj.layer{frm}.id = obj.layer{frm}.id(1:size(obj.layer{frm}.twtt,1));
        obj.layer_modified(frm) = true;
      end
      if size(obj.layer{frm}.twtt,1) > length(obj.layer{frm}.id)
        Nlayers = size(obj.layer{frm}.twtt,1) - length(obj.layer{frm}.id);
        if isempty(obj.layer_organizer.lyr_id)
          new_ids = 1:Nlayers;
        else
          new_ids = max(obj.layer_organizer.lyr_id) + (1:Nlayers);
        end
        obj.layer{frm}.id(end+(1:Nlayers)) = new_ids;
        obj.layer_modified(frm) = true;
      end
      Nlayers = length(obj.layer{frm}.id);
      if size(obj.layer{frm}.quality,1) < Nlayers
        obj.layer{frm}.quality(end+1:Nlayers) = ones(Nlayers,size(obj.layer{frm}.quality,2),'uint8');
        obj.layer_modified(frm) = true;
      end
      if size(obj.layer{frm}.quality,1) > Nlayers
        obj.layer{frm}.quality = obj.layer{frm}.quality(1:length(obj.layer{frm}.id));
        obj.layer_modified(frm) = true;
      end
      if size(obj.layer{frm}.twtt,1) < Nlayers
        obj.layer{frm}.twtt(end+1:Nlayers) = 2*ones(Nlayers,size(obj.layer{frm}.twtt,2),'uint8');
        obj.layer_modified(frm) = true;
      end
      if size(obj.layer{frm}.twtt,1) > Nlayers
        obj.layer{frm}.twtt = obj.layer{frm}.twtt(1:length(obj.layer{frm}.id));
        obj.layer_modified(frm) = true;
      end
      
      % Ensure that number of range lines is correct for each field
      if size(obj.layer{frm}.twtt,2) < Nx
        Nx_new = Nx - size(obj.layer{frm}.twtt,2);
        obj.layer{frm}.twtt(:,end+(1:Nx_new)) = nan(size(obj.layer{frm}.twtt,1),Nx_new);
        obj.layer_modified(frm) = true;
      end
      if size(obj.layer{frm}.type,2) < Nx
        Nx_new = Nx - size(obj.layer{frm}.type,2);
        obj.layer{frm}.type(:,end+(1:Nx_new)) = 2*ones(size(obj.layer{frm}.twtt,1),Nx_new,'uint8');
        obj.layer_modified(frm) = true;
      end
      if size(obj.layer{frm}.quality,2) < Nx
        Nx_new = Nx - size(obj.layer{frm}.quality,2);
        obj.layer{frm}.quality(:,end+(1:Nx_new)) = ones(size(obj.layer{frm}.twtt,1),Nx_new,'uint8');
      end
      if size(obj.layer{frm}.twtt,2) > Nx
        obj.layer{frm}.twtt(:,Nx+1:end) = [];
        obj.layer_modified(frm) = true;
      end
      if size(obj.layer{frm}.type,2) > Nx
        obj.layer{frm}.type(:,Nx+1:end) = [];
        obj.layer_modified(frm) = true;
      end
      if size(obj.layer{frm}.quality,2) > Nx
        obj.layer{frm}.quality(:,Nx+1:end) = [];
        obj.layer_modified(frm) = true;
      end

      % Ensure all values are valid
      mask = ~isfinite(obj.layer{frm}.quality < 1 | obj.layer{frm}.quality > 3);
      if any(mask)
        obj.layer{frm}.quality(mask) = 1;
        obj.layer_modified(frm) = true;
      end
      mask = ~isfinite(obj.layer{frm}.type < 1 | obj.layer{frm}.type > 4);
      if any(mask)
        obj.layer{frm}.type(mask) = 2;
        obj.layer_modified(frm) = true;
      end

      % Ensure all layers are present in the layer_organizer
      for lay_idx = 1:length(obj.layer{frm}.id)
        % Check to see if this is a duplicate obj.layer{frm} ID in this obj.layer{frm} file
        match_idx = find(obj.layer{frm}.id(lay_idx) == obj.layer{frm}.id(1:lay_idx-1));
        if ~isempty(match_idx)
          % This is a duplicate obj.layer{frm} ID: Create new obj.layer{frm} in obj.layer{frm}
          % organizer for this duplicate obj.layer{frm}
          obj.layer{frm}.id(lay_idx) = obj.insert_layers(struct('lyr_name',{{''}}));
          obj.layer_modified(frm) = true;
        else
          match_idx = find(obj.layer_organizer.lyr_id == obj.layer{frm}.id(lay_idx));
          if isempty(match_idx)
            lyr_organizer = [];
            if lay_idx == 1
              lyr_organizer.lyr_name = {'surface'};
              lyr_organizer.lyr_group_name = {'standard'};
            elseif lay_idx == 2
              lyr_organizer.lyr_name = {'bottom'};
              lyr_organizer.lyr_group_name = {'standard'};
            else
              lyr_organizer.lyr_name = {'auto'};
              lyr_organizer.lyr_group_name = {''};
            end
            obj.layer{frm}.id(lay_idx) = obj.insert_layers(lyr_organizer);
            obj.layer_modified(frm) = true;
          end
        end
      end
      
      % Check param fields
      if ~isfield(obj.layer{frm},'param') || ~isstruct(obj.layer{frm}.param)
        obj.layer{frm}.param = [];
        obj.layer_modified(frm) = true;
      end
      if ~isfield(obj.layer{frm}.param,'radar_name') || ~strcmp(obj.layer{frm}.param.radar_name,obj.param.radar_name)
        obj.layer{frm}.param.radar_name = obj.param.radar_name;
        obj.layer_modified(frm) = true;
      end
      if ~isfield(obj.layer{frm}.param,'season_name') || ~strcmp(obj.layer{frm}.param.season_name,obj.param.season_name)
        obj.layer{frm}.param.season_name = obj.param.season_name;
        obj.layer_modified(frm) = true;
      end
      if ~isfield(obj.layer{frm}.param,'day_seg') || ~strcmp(obj.layer{frm}.param.day_seg,obj.param.day_seg)
        obj.layer{frm}.param.day_seg = obj.param.day_seg;
        obj.layer_modified(frm) = true;
      end
      if ~isfield(obj.layer{frm}.param,'sw_version')
        obj.layer{frm}.param.sw_version = obj.param.sw_version;
        obj.layer_modified(frm) = true;
      end
      if ~isfield(obj.layer{frm}.param,'records') || ~isstruct(obj.layer{frm}.param.records)
        obj.layer{frm}.param.records = [];
        obj.layer_modified(frm) = true;
      end
      if ~isfield(obj.layer{frm}.param.records,'gps') || ~isstruct(obj.layer{frm}.param.records.gps)
        obj.layer{frm}.param.records.gps = [];
        obj.layer_modified(frm) = true;
      end
      if ~isfield(obj.layer{frm}.param.records.gps,'time_offset') || isempty(obj.layer{frm}.param.records.gps.time_offset)
        obj.layer{frm}.param.records.gps.time_offset = obj.param.records.gps.time_offset;
        obj.layer_modified(frm) = true;
      end
      if ~isfield(obj.layer{frm}.param,'radar') || ~isstruct(obj.layer{frm}.param.radar)
        obj.layer{frm}.param.radar = [];
        obj.layer_modified(frm) = true;
      end
      if ~isfield(obj.layer{frm}.param.radar,'lever_arm_fh') || isempty(obj.layer{frm}.param.radar.lever_arm_fh)
        obj.layer{frm}.param.radar.lever_arm_fh = obj.param.radar.lever_arm_fh;
        obj.layer_modified(frm) = true;
      end
      
    end

    %% fix_layer_organizer: fix any problems with layer organizer
    function fix_layer_organizer(obj)
      % Check for field existence
      if ~isfield(obj.layer_organizer,'lyr_age')
        obj.layer_organizer.lyr_age = [];
      end
      if ~isfield(obj.layer_organizer,'lyr_age_source')
        obj.layer_organizer.lyr_age_source = {};
      end
      if ~isfield(obj.layer_organizer,'lyr_desc')
        obj.layer_organizer.lyr_desc = {};
      end
      if ~isfield(obj.layer_organizer,'lyr_group_name')
        obj.layer_organizer.lyr_group_name = {};
      end
      if ~isfield(obj.layer_organizer,'lyr_id')
        obj.layer_organizer.lyr_id = [];
      end
      if ~isfield(obj.layer_organizer,'lyr_name')
        obj.layer_organizer.lyr_name = {};
      end
      if ~isfield(obj.layer_organizer,'lyr_order')
        obj.layer_organizer.lyr_order = [];
      end
      % Check the type
      if ~isa(obj.layer_organizer.lyr_age,'double')
        obj.layer_organizer.lyr_age = double(obj.layer_organizer.lyr_age);
        obj.layer_organizer_modified = true;
      end
      if ~iscell(obj.layer_organizer.lyr_age_source)
        obj.layer_organizer.lyr_age_source = {};
        obj.layer_organizer_modified = true;
      else
        for cell_idx = 1:length(obj.layer_organizer.lyr_age_source)
          if ~isstruct(obj.layer_organizer.lyr_age_source{cell_idx})
            obj.layer_organizer.lyr_age_source{cell_idx} = struct('age',{},'source',{},'type',{});
            obj.layer_organizer_modified = true;
          else
            for struct_idx = 1:length(obj.layer_organizer.lyr_age_source{cell_idx})
              if ~isa(obj.layer_organizer.lyr_age_source{cell_idx}(struct_idx).age,'double')
                obj.layer_organizer.lyr_age_source{cell_idx}(struct_idx).age = double(obj.layer_organizer.lyr_age_source{cell_idx}(struct_idx).age);
                obj.layer_organizer_modified = true;
              end
              if ~ischar(obj.layer_organizer.lyr_age_source{cell_idx}(struct_idx).source)
                obj.layer_organizer.lyr_age_source{cell_idx}(struct_idx).source = char(obj.layer_organizer.lyr_age_source{cell_idx}(struct_idx).source);
                obj.layer_organizer_modified = true;
              end
              if ~isa(obj.layer_organizer.lyr_age_source{cell_idx}(struct_idx).type,'double')
                obj.layer_organizer.lyr_age_source{cell_idx}(struct_idx).type = double(obj.layer_organizer.lyr_age_source{cell_idx}(struct_idx).type);
                obj.layer_organizer_modified = true;
              end
              if ~any(obj.layer_organizer.lyr_age_source{cell_idx}(struct_idx).type == [0 1 2 3 4])
                obj.layer_organizer.lyr_age_source{cell_idx}(struct_idx).type = 0;
                obj.layer_organizer_modified = true;
              end
            end
          end
        end
      end
      if ~iscell(obj.layer_organizer.lyr_desc)
        obj.layer_organizer.lyr_desc = {};
        obj.layer_organizer_modified = true;
      else
        for cell_idx = 1:length(obj.layer_organizer.lyr_desc)
          if ~ischar(obj.layer_organizer.lyr_desc{cell_idx})
            obj.layer_organizer.lyr_desc{cell_idx} = char(obj.layer_organizer.lyr_desc{cell_idx});
            obj.layer_organizer_modified = true;
          end
        end
      end
      if ~iscell(obj.layer_organizer.lyr_group_name)
        obj.layer_organizer.lyr_group_name = {};
        obj.layer_organizer_modified = true;
      else
        for cell_idx = 1:length(obj.layer_organizer.lyr_group_name)
          if ~ischar(obj.layer_organizer.lyr_group_name{cell_idx})
            obj.layer_organizer.lyr_group_name{cell_idx} = char(obj.layer_organizer.lyr_group_name{cell_idx});
            obj.layer_organizer_modified = true;
          end
        end
      end
      if ~isa(obj.layer_organizer.lyr_id,'double')
        obj.layer_organizer.lyr_id = double(obj.layer_organizer.lyr_id);
        obj.layer_organizer_modified = true;
      end
      if ~iscell(obj.layer_organizer.lyr_name)
        obj.layer_organizer.lyr_name = {};
        obj.layer_organizer_modified = true;
      else
        for cell_idx = 1:length(obj.layer_organizer.lyr_name)
          if ~ischar(obj.layer_organizer.lyr_name{cell_idx})
            obj.layer_organizer.lyr_name{cell_idx} = char(obj.layer_organizer.lyr_name{cell_idx});
            obj.layer_organizer_modified = true;
          end
        end
      end
      if ~isa(obj.layer_organizer.lyr_order,'double')
        obj.layer_organizer.lyr_order = double(obj.layer_organizer.lyr_order);
        obj.layer_organizer_modified = true;
      end
      
      % Check the length (too long)
      Nlayers = length(obj.layer_organizer.lyr_id);
      if length(obj.layer_organizer.lyr_age) > Nlayers
        obj.layer_organizer.lyr_age(Nlayers+1:end) = [];
        obj.layer_organizer_modified = true;
      end
      if length(obj.layer_organizer.lyr_age_source) > Nlayers
        obj.layer_organizer.lyr_age_source(Nlayers+1:end) = [];
        obj.layer_organizer_modified = true;
      end
      if length(obj.layer_organizer.lyr_desc) > Nlayers
        obj.layer_organizer.lyr_desc(Nlayers+1:end) = [];
        obj.layer_organizer_modified = true;
      end
      if length(obj.layer_organizer.lyr_group_name) > Nlayers
        obj.layer_organizer.lyr_group_name(Nlayers+1:end) = [];
        obj.layer_organizer_modified = true;
      end
      if length(obj.layer_organizer.lyr_name) > Nlayers
        obj.layer_organizer.lyr_name(Nlayers+1:end) = [];
        obj.layer_organizer_modified = true;
      end
      if length(obj.layer_organizer.lyr_order) > Nlayers
        obj.layer_organizer.lyr_order(Nlayers+1:end) = [];
        obj.layer_organizer_modified = true;
      end
      
      % Check the length (too short)
      if length(obj.layer_organizer.lyr_age) < Nlayers
        obj.layer_organizer.lyr_age(end+1:Nlayers) = NaN;
        obj.layer_organizer_modified = true;
      end
      if length(obj.layer_organizer.lyr_age_source) < Nlayers
        obj.layer_organizer.lyr_age_source(end+1:Nlayers) = cellfun(@struct,cell(1,Nlayers-length(obj.layer_organizer.lyr_age_source)),'UniformOutput',false);
        obj.layer_organizer_modified = true;
      end
      if length(obj.layer_organizer.lyr_desc) < Nlayers
        obj.layer_organizer.lyr_desc(end+1:Nlayers) = cellfun(@char,cell(1,Nlayers-length(obj.layer_organizer.lyr_desc)),'UniformOutput',false);
        obj.layer_organizer_modified = true;
      end
      if length(obj.layer_organizer.lyr_group_name) < Nlayers
        obj.layer_organizer.lyr_group_name(end+1:Nlayers) = cellfun(@char,cell(1,Nlayers-length(obj.layer_organizer.lyr_group_name)),'UniformOutput',false);
        obj.layer_organizer_modified = true;
      end
      if length(obj.layer_organizer.lyr_name) < Nlayers
        obj.layer_organizer.lyr_name(end+1:Nlayers) = cellfun(@char,cell(1,Nlayers-length(obj.layer_organizer.lyr_name)),'UniformOutput',false);
        obj.layer_organizer_modified = true;
      end
      if length(obj.layer_organizer.lyr_order) < Nlayers
        if isempty(obj.layer_organizer.lyr_order)
          new_orders = 1:Nlayers;
        else
          new_orders = max(obj.layer_organizer.lyr_order) + (1:Nlayers);
        end
        obj.layer_organizer.lyr_order(end+1:Nlayers) = new_orders;
        obj.layer_organizer_modified = true;
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
        obj.layer_organizer_modified = true;
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
          obj.layer_organizer_modified = true;
        end
      end
      
      % Check param fields
      if ~isfield(obj.layer_organizer,'param') || ~isstruct(obj.layer_organizer.param)
        obj.layer_organizer.param = [];
        obj.layer_organizer_modified = true;
      end
      if ~isfield(obj.layer_organizer.param,'radar_name') || ~strcmp(obj.layer_organizer.param.radar_name,obj.param.radar_name)
        obj.layer_organizer.param.radar_name = obj.param.radar_name;
        obj.layer_organizer_modified = true;
      end
      if ~isfield(obj.layer_organizer.param,'season_name') || ~strcmp(obj.layer_organizer.param.season_name,obj.param.season_name)
        obj.layer_organizer.param.season_name = obj.param.season_name;
        obj.layer_organizer_modified = true;
      end
      if ~isfield(obj.layer_organizer.param,'day_seg') || ~strcmp(obj.layer_organizer.param.day_seg,obj.param.day_seg)
        obj.layer_organizer.param.day_seg = obj.param.day_seg;
        obj.layer_organizer_modified = true;
      end
      if ~isfield(obj.layer_organizer.param,'sw_version')
        obj.layer_organizer.param.sw_version = obj.param.sw_version;
        obj.layer_organizer_modified = true;
      end
    end
    
    %% get_id: get layer id from layer name
    function id = get_id(obj,name)
      obj.check_layer_organizer();
      match_idx = find(strcmp(name,obj.layer_organizer.lyr_name));
      id = obj.layer_organizer.lyr_id(match_idx);
    end
    
    %% get_layer: get layer twtt, quality, type
    function [twtt,quality,type] = get_layer(obj,frms,id)
      obj.check_layer_organizer();
      if ischar(id)
        % name passed in rather than id
        match_idx = find(strcmpi(id,obj.layer_organizer.lyr_name));
        if isempty(id)
          error('Layer does not exist in layer organizer. Run insert_layers() first.');
        end
        id = obj.layer_organizer.lyr_id(match_idx);
      else
        % id passed in
        if all(obj.layer_organizer.lyr_id ~= id)
          error('Layer does not exist in layer organizer. Run insert_layers() first.');
        end
      end
      twtt = [];
      quality = [];
      type = [];
      for frm = 1:frms(:).'
        obj.check(frm);
        lay_idx = find(obj.layer{frm}.id == id);
        Nx = length(obj.layer{frm}.gps_time);
        if isempty(lay_idx)
          % Layer does not exist in file, set to defaults
          twtt(end+(1:Nx)) = NaN;
          quality(end+(1:Nx)) = 1;
          type(end+(1:Nx)) = 2;
        else
          % Layer exists, get values
          twtt(end+(1:Nx)) = obj.layer{frm}.twtt(lay_idx,:);
          quality(end+(1:Nx)) = obj.layer{frm}.quality(lay_idx,:);
          type(end+(1:Nx)) = obj.layer{frm}.type(lay_idx,:);
        end
      end
    end
    
    %% gps_time: get gps_time
    function gps_time = gps_time(obj,frms)
      gps_time = [];
      for frm = 1:frms(:).'
        obj.check(frm);
        gps_time(end+(1:length(obj.layer{frm}.gps_time))) = obj.layer{frm}.gps_time;
      end
    end
    
    %% insert_layers: into layer organizer
    function ids = insert_layers(obj,layer_organizer)
      obj.check_layer_organizer();
      
      Nlayers = length(layer_organizer.lyr_name);
      if isempty(obj.layer_organizer.lyr_id)
        ids = 1:Nlayers;
      else
        ids = max(obj.layer_organizer.lyr_id) + (1:Nlayers);
      end
      obj.layer_organizer.lyr_id(end+(1:Nlayers)) = ids;
      obj.layer_organizer.lyr_name(end+(1:Nlayers)) = layer_organizer.lyr_name;
      
      if isfield(layer_organizer,'lyr_age') && length(layer_organizer.lyr_age) == Nlayers
        obj.layer_organizer.lyr_age(end+(1:Nlayers)) = layer_organizer.lyr_age;
      else
        obj.layer_organizer.lyr_age(end+(1:Nlayers)) = NaN;
      end
      if isfield(layer_organizer,'lyr_age_source') && length(layer_organizer.lyr_age_source) == Nlayers
        obj.layer_organizer.lyr_age_source(end+(1:Nlayers)) = layer_organizer.lyr_age_source;
      else
        obj.layer_organizer.lyr_age_source(end+(1:Nlayers)) = cellfun(@struct,cell(1,Nlayers),'UniformOutput',false);
      end
      if isfield(layer_organizer,'lyr_desc') && length(layer_organizer.lyr_desc) == Nlayers
        obj.layer_organizer.lyr_desc(end+(1:Nlayers)) = layer_organizer.lyr_desc;
      else
        obj.layer_organizer.lyr_desc(end+(1:Nlayers)) = cellfun(@char,cell(1,Nlayers),'UniformOutput',false);
      end
      if isfield(layer_organizer,'lyr_group_name') && length(layer_organizer.lyr_group_name) == Nlayers
        obj.layer_organizer.lyr_group_name(end+(1:Nlayers)) = layer_organizer.lyr_group_name;
      else
        obj.layer_organizer.lyr_group_name(end+(1:Nlayers)) = cellfun(@char,cell(1,Nlayers),'UniformOutput',false);
      end
      if isfield(layer_organizer,'lyr_order') && length(layer_organizer.lyr_order) == Nlayers
        obj.layer_organizer.lyr_order(end+(1:Nlayers)) = layer_organizer.lyr_order;
      else
        if isempty(obj.layer_organizer.lyr_order)
          new_orders = 1:Nlayers;
        else
          new_orders = max(obj.layer_organizer.lyr_order) + (1:Nlayers);
        end
        obj.layer_organizer.lyr_order(end+(1:Nlayers)) = new_orders;
      end
      obj.fix_layer_organizer();
      obj.layer_organizer_modified = true;
    end
    
    %% lat: get lat
    function lat = lat(obj,frms)
      lat = [];
      for frm = 1:frms(:).'
        obj.check(frm);
        lat(end+(1:length(obj.layer{frm}.lat))) = obj.layer{frm}.lat;
      end
    end
    
    %% layer_fn: get layer filename
    function fn = layer_fn(obj,frm)
      fn = fullfile(ct_filename_out(obj.param,'',obj.layerdata_source),sprintf('Data_%s_%03d.mat',obj.param.day_seg,frm));
    end
    
    %% layer_organizer_fn: get layer organizer filename
    function fn = layer_organizer_fn(obj,frm)
      fn = fullfile(ct_filename_out(obj.param,'',obj.layerdata_source),sprintf('layer_%s.mat',obj.param.day_seg));
    end
    
    %% load: load layer
    function load(obj,frm)
      layer_fn = fullfile(ct_filename_out(obj.param,'',obj.layerdata_source),sprintf('Data_%s_%03d.mat',obj.param.day_seg,frm));
      if ~exist(layer_fn,'file')
        obj.create(frm);
      else
        obj.layer{frm} = load(layer_fn);
        obj.layer_modified(frm) = false;
        obj.fix(frm);
      end
    end
    
    %% load_layer_organizer: load layer organizer
    function load_layer_organizer(obj)
      layer_organizer_fn = fullfile(ct_filename_out(obj.param,'',obj.layerdata_source),sprintf('layer_%s.mat',obj.param.day_seg));
      if ~exist(layer_organizer_fn,'file')
        obj.create_layer_organizer();
      else
        obj.layer_organizer = load(layer_organizer_fn);
        obj.layer_organizer_modified = false;
        obj.fix_layer_organizer();
      end
    end
    
    %% load_records: load records and frames
    function load_records(obj,records)
      
      if ~exist('records','var') || isempty(records)
        % Load records file
        records_fn = ct_filename_support(obj.param,'','records');
        obj.records = load(records_fn);
      end
      
      % Load frames file
      obj.frames = frames_load(obj.param);
      
      % Create reference trajectory (rx_path == 0, tx_weights = []). Update
      % the records field with this information.
      trajectory_param = struct('gps_source',obj.records.gps_source, ...
        'season_name',obj.param.season_name,'radar_name',obj.param.radar_name,'rx_path', 0, ...
        'tx_weights', [], 'lever_arm_fh', obj.param.radar.lever_arm_fh);
      obj.records = trajectory_with_leverarm(obj.records,trajectory_param);
      
      obj.records.along_track = geodetic_to_along_track(obj.records.lat,obj.records.lon);
      
      if isempty(obj.along_track_spacing)
        sys = ct_output_dir(obj.param.radar_name);
        switch (sys)
          case {'rds','accum'}
            obj.along_track_spacing = 15;
          otherwise
            obj.along_track_spacing = 5;
        end
      end
      obj.along_track = [0:obj.along_track_spacing:obj.records.along_track(end)];

      if isempty(obj.sample_spacing_method)
        if length(obj.along_track) < 2 && length(obj.records.gps_time) > obj.record_spacing
          obj.sample_spacing_method = 'records';
        else
          obj.sample_spacing_method = 'along_track';
        end
      end

    end
    
    %% lon: get lon
    function lon = lon(obj,frms)
      lon = [];
      for frm = 1:frms(:).'
        obj.check(frm);
        lon(end+(1:length(obj.layer{frm}.lon))) = obj.layer{frm}.lon;
      end
    end
    
    %% save: save changes to layers
    function save(obj, layerdata_source)
      if exist('layerdata_source','var') && ~isempty(layerdata_source)
        obj.layer_organizer_modified = true;
        obj.layerdata_source = layerdata_source;
        obj.check_records();
        obj.layer_modified(1:length(obj.frames.frame_idxs)) = true;
      end
      
      if obj.layer_organizer_modified == true
        layer_organizer_fn = fullfile(ct_filename_out(obj.param,'',obj.layerdata_source),sprintf('layer_%s.mat',obj.param.day_seg));
        layer_organizer = obj.layer_organizer;
        fprintf('Saving %s\n', layer_organizer_fn);
        ct_save(layer_organizer_fn,'-struct','layer_organizer');
        obj.layer_organizer_modified = false;
      end

      for frm = 1:length(obj.layer_modified)
        if obj.layer_modified(frm) == true
          layer_fn = fullfile(ct_filename_out(obj.param,'',obj.layerdata_source),sprintf('Data_%s_%03d.mat',obj.param.day_seg,frm));
          layer = obj.layer{frm};
          fprintf('Saving %s\n', layer_fn);
          ct_save(layer_fn,'-struct','layer');
          obj.layer_modified(frm) = false;
        end
      end
    end
    
    %% update_gps: update gps for a frame
    function update_gps(obj,frms)
      obj.check_records();

      for frm = frms(:).'
        if strcmp(obj.sample_spacing_method,'records')
          if frm < length(obj.frames.frame_idxs)
            recs = obj.frames.frame_idxs(frm):obj.record_spacing:obj.frames.frame_idxs(frm+1);
          else
            recs = obj.frames.frame_idxs(frm):obj.record_spacing:obj.frames.Nx;
          end
          obj.layer{frm}.gps_time = obj.records.gps_time(recs);
          obj.layer{frm}.elev = obj.records.elev(recs);
          obj.layer{frm}.lat = obj.records.lat(recs);
          obj.layer{frm}.lon = obj.records.lon(recs);
        else
          if frm < length(obj.frame_idxs)
            frm_mask = find(obj.along_track >= obj.records.along_track(obj.frames.frame_idxs(frm)) ...
              & obj.along_track < obj.records.along_track(obj.frames.frame_idxs(frm+1)));
          else
            frm_mask = find(obj.along_track >= obj.records.along_track(obj.frames.frame_idxs(frm)));
          end
          good_idxs = [true diff(obj.records.along_track)>0];
          obj.layer{frm}.gps_time = interp1(obj.records.along_track(good_idxs),obj.records.gps_time(good_idxs),obj.along_track(frm_mask));
          obj.layer{frm}.elev = interp1(obj.records.along_track(good_idxs),obj.records.elev(good_idxs),obj.along_track(frm_mask));
          obj.layer{frm}.lat = interp1(obj.records.along_track(good_idxs),obj.records.lat(good_idxs),obj.along_track(frm_mask));
          obj.layer{frm}.lon = interp1(obj.records.along_track(good_idxs),obj.records.lon(good_idxs),obj.along_track(frm_mask));
        end
        obj.layer{frm}.param.radar.lever_arm_fh = obj.param.radar.lever_arm_fh;
        obj.layer{frm}.param.records.gps.time_offset = obj.param.records.gps.time_offset;
        obj.layer{frm}.param.sw_version = obj.param.sw_version;
        obj.layer_modified(frm) = true;
      end
      
    end
    
    %% update_gps_all: update gps for a frame
    function update_gps_all(obj)
      obj.check_records();
      obj.update_gps(1:length(obj.frames.frame_idxs));
    end
    
    %% update_layers: update layer's twtt, quality, type
    function update_layer(obj,frms,id,gps_time,twtt,quality,type)
      obj.check_layer_organizer();
      if ischar(id)
        % name passed in rather than id
        id = find(strcmpi(id,obj.layer_organizer.lyr_name));
        if isempty(id)
          error('Layer does not exist in layer organizer. Run insert_layers() first.');
        end
      else
        if all(obj.layer_organizer.lyr_id ~= id)
          error('Layer does not exist in layer organizer. Run insert_layers() first.');
        end
      end
      
      for frm = 1:frms(:).'
        obj.check(frm);
        lay_idx = find(obj.layer{frm}.id == id);
        if isempty(lay_idx)
          lay_idx = length(obj.layer{frm}.id) + 1;
          obj.layer{frm}.id(lay_idx) = id;
        end
        Nx = length(obj.layer{frm}.gps_time);
        obj.layer{frm}.twtt(lay_idx,1:Nx) = interp1(gps_time,twtt,obj.layer{frm}.gps_time);
        if ~exist('quality','var') || isempty(quality)
          quality = ones(size(gps_time),'uint8');
        end
        quality = uint8(quality);
        quality(quality<1 | quality>3) = 1;
        obj.layer{frm}.quality(lay_idx,1:Nx) = uint8(interp1(gps_time,double(quality),double(obj.layer{frm}.gps_time),'nearest'));
        if ~exist('type','var') || isempty(type)
          type = 2*ones(size(gps_time),'uint8');
        end
        type = uint8(type);
        type(type<1 | type>4) = 2;
        obj.layer{frm}.type(lay_idx,1:Nx) = uint8(interp1(gps_time,double(type),double(obj.layer{frm}.gps_time),'nearest'));
        obj.layer_modified(frm) = true;
      end
    end
    
  end
  
end
