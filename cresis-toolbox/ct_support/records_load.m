function records = records_load(param,varargin)
% records = records_load(param,varargin)
%
% records = records(param)
%   Loads records file.
%
% param: parameter spreadsheet structure
%
% records = records(param,recs)
%   Loads only the specified records from the records file.
%
% param: parameter spreadsheet structure
% recs: 2x1 vector specifying the start/stop records [start_record, stop_record]
%   Records are one-indexed (first records is record "1")
%   Set stop record to inf to read to end of file
%   "inf" can be used for start or stop record to mean last record
%
% records = records(param,FIELD_STRING_ARGS ...)
%   Loads only the specified records from the records file. This
%   functionality is slower, but reduces memory usage.
%
% param: parameter spreadsheet structure
% FIELD_STRING_ARGS: Argument list of strings containing specific fields to load from the records file.
%
% Example
%   param = read_param_xls(ct_filename_param('rds_param_2019_Antarctica_Ground.xls'),'20191231_04');
%   records = records_load(param);
%
%   records = records_load(param,[1 100]);
%
%   records = records_load(param,'lat','lon');
%
% Author: John Paden

if ischar(param)
  records = load(param,'param_records');
  param = records.param_records;
end
records_fn = ct_filename_support(param,'','records');
if ~exist(records_fn,'file')
  error('Records file does not exist: %s (%s).', records_fn, datestr(now));
end

%% Correct old files
% =========================================================================
update_records_flag = false;
records = load(records_fn,'settings');
if isfield(records,'settings') && isfield(records.settings,'wfs') && isfield(records.settings.wfs,'wfs')
  update_records_flag = true;
end
if isfield(records,'settings') && isfield(records.settings,'nyquist_zone_hw')
  update_records_flag = true;
end
if isfield(records,'settings') && isfield(records.settings,'nyquist_zone')
  update_records_flag = true;
end

mat_vars = whos('-file',records_fn);

if any(strcmp('vectors',{mat_vars.name}))
  update_records_flag = true;
end

if any(strcmp('surface',{mat_vars.name}))
%   update_records_flag = true;
end

if ~update_records_flag
  try
    records = load(records_fn,'param_records');
    delta_offset = max(param.records.gps.time_offset) - max(records.param_records.records.gps.time_offset);
    if delta_offset ~= 0
      update_records_flag = true;
    end
  end
end

if update_records_flag
  records_update(param,[]);
end

if nargin == 1
  %% Load all fields
  % =======================================================================
  records = load(records_fn);
  if ~isfield(records,'bit_mask')
    records.bit_mask = zeros(size(records.offset),'uint8');
  end
  
elseif nargin == 2 && isnumeric(varargin{1})
  %% Load specific records
  % =======================================================================

  recs = varargin{1};
  match_idx = find(strcmp('gps_time',{mat_vars.name}));
  num_recs = mat_vars(match_idx).size(2);
  if recs(2) == inf
    recs(2) = num_recs;
  end
  if recs(2) > num_recs
    error('Requested records beyond the end of the records file: recs(2)=%d > num_recs=%d. There are %d records in the records file: %s.', recs(2), num_recs, num_recs, records_fn);
  end
  if recs(2) < recs(1)
    error('Requested records beyond the end of the records file or requested zero records: recs(1)=%d > recs(2)=%d. There are %d records in the records file: %s.', recs(1), recs(2), num_recs, records_fn);
  end

  records = load(records_fn,'gps_source','param_records','relative_filename','relative_rec_num','settings');
  
  records_mat = matfile(records_fn);
  records.elev = records_mat.elev(:,recs(1):recs(2));
  records.gps_time = records_mat.gps_time(:,recs(1):recs(2));
  records.heading = records_mat.heading(:,recs(1):recs(2));
  records.lat = records_mat.lat(:,recs(1):recs(2));
  records.lon = records_mat.lon(:,recs(1):recs(2));
  records.offset = records_mat.offset(:,recs(1):recs(2));
  records.pitch = records_mat.pitch(:,recs(1):recs(2));
  records.roll = records_mat.roll(:,recs(1):recs(2));
  
  if any(strcmp('bit_mask',{mat_vars.name}))
    records.bit_mask = records_mat.bit_mask(:,recs(1):recs(2));
  else
    records.bit_mask = zeros(size(records.offset),'uint8');
  end
  
  if any(strcmp('phase_correction',{mat_vars.name}))
    records.phase_correction = records_mat.phase_correction(:,recs(1):recs(2));
  else
    records.phase_correction = zeros(size(records.gps_time));
  end
  
  if any(strcmp('Tadc_correction',{mat_vars.name}))
    records.Tadc_correction = records_mat.Tadc_correction(:,recs(1):recs(2));
  else
    records.Tadc_correction = zeros(size(records.gps_time));
  end
  
  if any(strcmp('nyquist_zone_sig',{mat_vars.name}))
    records.nyquist_zone_sig = records_mat.nyquist_zone_sig(:,recs(1):recs(2));
  else
    records.nyquist_zone_sig = zeros(size(records.gps_time));
  end
  
  if any(strcmp('nyquist_zone_hw',{mat_vars.name}))
    records.nyquist_zone_hw = records_mat.nyquist_zone_hw(:,recs(1):recs(2));
  else
    records.nyquist_zone_hw = zeros(size(records.gps_time));
  end

elseif nargin > 1 && all(cellfun(@ischar,varargin))
  %% Load specific fields
  records = load(records_fn,varargin{:});
  
else
  error('Invalid arguments to records_load.');
end
