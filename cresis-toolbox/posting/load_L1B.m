function mdata = load_L1B(fn)
% function load_L1B(fn)
%
% fn = filename string and must contain correct extension ('.nc' or '.mat')
%
% Loads L1B cresis echogram files.
% L1B files can be in netcdf or mat format.
% They can be "compressed" (compression is used for Snow and Kuband radar data
% to make them smaller... main difference is that the surface is tracked
% so that the range gate stored in the file changes on each range line).
%
% Example:
%  fn = 'IRMCR1B_V01_20130408_01_020.nc';
%  mdata = load_L1B(fn);
%
%  fn = 'Data_20111107_02_191.mat';
%  mdata = load_L1B(fn);
%
% Author: John Paden
%
% See also: plot_L1B.m, uncompress_echogram.m

[fn_dir,fn_name,fn_ext] = fileparts(fn);

if strcmpi(fn_ext,'.nc')
  mdata = netcdf_to_mat(fn);
  
  mdata.Latitude = mdata.lat;
  mdata = rmfield(mdata,'lat');
  mdata.Longitude = mdata.lon;
  mdata = rmfield(mdata,'lon');
  mdata.Elevation = mdata.altitude;
  mdata = rmfield(mdata,'altitude');
  
  mdata.Data = 10.^(mdata.amplitude/10);
  mdata = rmfield(mdata,'amplitude');
  
  mdata.GPS_time = mdata.time;
  epoch_date_str = ncreadatt(fn, 'time', 'units');
  mdata.GPS_time = mdata.GPS_time + datenum_to_epoch(datenum(epoch_date_str(15:end)));
  mdata.GPS_time = mdata.GPS_time + utc_leap_seconds(mdata.GPS_time(1))
  mdata = rmfield(mdata,'time');
  
  mdata.Time = mdata.fasttime/1e6;
  mdata = rmfield(mdata,'fasttime');
  
  mdata.Roll = mdata.roll/180*pi;
  mdata = rmfield(mdata,'roll');
  mdata.Pitch = mdata.pitch/180*pi;
  mdata = rmfield(mdata,'pitch');
  mdata.Heading = mdata.heading/180*pi;
  mdata = rmfield(mdata,'heading');
  
elseif strcmpi(fn_ext,'.mat')
  mdata = load(fn);
  
else
  error('Unsupported extention %s', fn_ext);
end

mdata = uncompress_echogram(mdata);

if isfield(mdata,'param_get_heights')
  if isequal({'get_heights'},fieldnames(mdata.param_get_heights))
    % param_get_heights is incomplete structure, supplement with param_records
    mdata.param_qlook = mdata.param_records;
    mdata.param_qlook = rmfield(mdata.param_qlook,'get_heights');
    mdata.param_qlook.qlook = mdata.param_get_heights;
    mdata.file_version = '0';
    mdata.file_type = 'echo';
  else
    % param_get_heights is complete structure, just rename to qlook
    mdata.param_qlook = mdata.param_get_heights;
    mdata.param_qlook.qlook = mdata.param_qlook.get_heights;
    mdata.param_qlook = rmfield(mdata.param_qlook,'get_heights');
    % Ensure param_records has standard 3 fields (some old files do not
    % have)
    mdata.param_records.day_seg = mdata.param_qlook.day_seg;
    mdata.param_records.season_name = mdata.param_qlook.season_name;
    mdata.param_records.radar_name = mdata.param_qlook.radar_name;
    if isfield(mdata.param_records,'gps')
      mdata.param_records.records.gps = mdata.param_records.gps;
      mdata.param_qlook.records.gps = mdata.param_records.gps;
      mdata.param_records = rmfield(mdata.param_records,'gps');
    end
    if isfield(mdata.param_records,'file')
      mdata.param_records.records.file = mdata.param_records.file;
      mdata.param_records = rmfield(mdata.param_records,'file');
    end
    if isfield(mdata.param_records,'gps_source')
      mdata.param_records.records.gps_source = mdata.param_records.gps_source;
      mdata.param_records = rmfield(mdata.param_records,'gps_source');
    end
    if isfield(mdata.param_records,'records_fn')
      mdata.param_records = rmfield(mdata.param_records,'records_fn');
    end
    % Ensure both parameter structs have sw_version
    if ~isfield(mdata.param_records,'sw_version')
      mdata.param_records.sw_version.ver = '';
      mdata.param_records.sw_version.date_time = '';
      mdata.param_records.sw_version.cur_date_time = '';
      mdata.param_records.sw_version.rev = '';
      mdata.param_records.sw_version.URL = '';
    end
    if ~isfield(mdata.param_qlook,'sw_version')
      mdata.param_qlook.sw_version = mdata.param_records.sw_version;
    end
    % Ensure both parameter structs have radar.lever_arm_fh
    if ~isfield(mdata.param_records,'radar')
      mdata.param_records.radar = [];
    end
    if ~isfield(mdata.param_qlook,'radar')
      mdata.param_qlook.radar = [];
    end
    if ~isfield(mdata.param_records.radar,'lever_arm_fh')
      mdata.param_records.radar.lever_arm_fh = [];
    end
    if ~isfield(mdata.param_qlook.radar,'lever_arm_fh')
      mdata.param_qlook.radar.lever_arm_fh = [];
    end
    % Add file_version
    mdata.file_version = '0';
    mdata.file_type = 'echo';
  end
  mdata = rmfield(mdata,'param_get_heights');
  if isfield(mdata.param_qlook,'proc')
    mdata.param_qlook.load.frm = mdata.param_qlook.proc.frm;
  end
end

if isfield(mdata,'param_csarp')
  mdata.param_sar = mdata.param_csarp;
  mdata = rmfield(mdata,'param_csarp');
  mdata.param_sar.sar = mdata.param_sar.csarp;
  mdata.param_sar = rmfield(mdata.param_sar,'csarp');
  mdata.file_version = '0';
  mdata.file_type = 'echo';
end

if isfield(mdata,'param_combine')
  mdata.param_array = mdata.param_combine;
  mdata = rmfield(mdata,'param_combine');
  mdata.file_version = '0';
  mdata.file_type = 'echo';
end

if isfield(mdata,'param_combine_wf_chan')
  mdata.param_array.array = mdata.param_combine_wf_chan;
  mdata = rmfield(mdata,'param_combine_wf_chan');
  
  % Ensure param_records has standard 3 fields (some old files do not
  % have)
  mdata.param_records.day_seg = mdata.param_sar.day_seg;
  mdata.param_records.season_name = mdata.param_sar.season_name;
  mdata.param_records.radar_name = mdata.param_sar.radar_name;
  mdata.param_array.day_seg = mdata.param_sar.day_seg;
  mdata.param_array.season_name = mdata.param_sar.season_name;
  mdata.param_array.radar_name = mdata.param_sar.radar_name;
  if isfield(mdata.param_records,'gps')
    mdata.param_records.records.gps = mdata.param_records.gps;
    mdata.param_array.records.gps = mdata.param_records.gps;
    mdata.param_records = rmfield(mdata.param_records,'gps');
  end
  if isfield(mdata.param_records,'file')
    mdata.param_records.records.file = mdata.param_records.file;
    mdata.param_records = rmfield(mdata.param_records,'file');
  end
  if isfield(mdata.param_records,'gps_source')
    mdata.param_records.records.gps_source = mdata.param_records.gps_source;
    mdata.param_records = rmfield(mdata.param_records,'gps_source');
  end
  if isfield(mdata.param_records,'records_fn')
    mdata.param_records = rmfield(mdata.param_records,'records_fn');
  end
  % Ensure both parameter structs have sw_version
  if ~isfield(mdata.param_records,'sw_version')
    mdata.param_records.sw_version.ver = '';
    mdata.param_records.sw_version.date_time = '';
    mdata.param_records.sw_version.cur_date_time = '';
    mdata.param_records.sw_version.rev = '';
    mdata.param_records.sw_version.URL = '';
  end
  if ~isfield(mdata.param_array,'sw_version')
    mdata.param_array.sw_version = mdata.param_records.sw_version;
  end
  if ~isfield(mdata.param_array,'sw_version')
    mdata.param_sar.sw_version = mdata.param_records.sw_version;
  end
  % Ensure both parameter structs have radar.lever_arm_fh
  if ~isfield(mdata.param_records,'radar')
    mdata.param_records.radar = [];
  end
  if ~isfield(mdata.param_array,'radar')
    mdata.param_array.radar = [];
  end
  if ~isfield(mdata.param_sar,'radar')
    mdata.param_sar.radar = [];
  end
  if ~isfield(mdata.param_records.radar,'lever_arm_fh')
    mdata.param_records.radar.lever_arm_fh = [];
  end
  if ~isfield(mdata.param_array.radar,'lever_arm_fh')
    mdata.param_array.radar.lever_arm_fh = [];
  end
  if ~isfield(mdata.param_sar.radar,'lever_arm_fh')
    mdata.param_sar.radar.lever_arm_fh = [];
  end
  
  mdata.file_version = '0';
  mdata.file_type = 'echo';
end

if ~isfield(mdata,'Roll') && isfield(mdata,'GPS_time')
  mdata.Roll = zeros(size(mdata.GPS_time));
end

if ~isfield(mdata,'Pitch') && isfield(mdata,'GPS_time')
  mdata.Pitch = zeros(size(mdata.GPS_time));
end

if ~isfield(mdata,'Heading') && isfield(mdata,'GPS_time')
  mdata.Heading = zeros(size(mdata.GPS_time));
end

