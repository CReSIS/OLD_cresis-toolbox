function create_records(param, param_override)
% hfrds.create_records(param, param_override)
% create_records_mcords2(param, param_override)
%
% Function for creating records file for HFRDS data. The function is
% usually called from master.m but can also be called from
% hfrds.run_create_records.m.
%
% This function should be run after the GPS file has been created.
% For example, cresis-toobox/gps/missions/make_gps_2016_greenland_G1XB.m
%
% This function's output file is used by all other parts of the processing.
%
% param = struct with processing parameters
%           -- OR --
%         function handle to script which creates a structure with the
%           processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Author: Carl Leuschen, John Paden
%
% See also: run_master.m, master.m, hfrds.run_create_records.m, hfrds.create_records.m

%% General Setup
% =====================================================================

if ~isstruct(param)
  % Functional form
  param();
end
param = merge_structs(param, param_override);

dbstack_info = dbstack;
fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', dbstack_info(1).name, param.day_seg, datestr(now,'HH:MM:SS'));
fprintf('=====================================================================\n');

%% Prep work
% =====================================================================

% Get WGS84 ellipsoid parameters
physical_constants

%% Get the file
% =====================================================================
[base_dir,adc_folder_name,fns,file_idxs] = get_segment_file_list(param);

fn = fns{1};
[~,fn_name,fn_ext] = fileparts(fn);
fprintf('  Parsing file %s (%s)\n', fn_name, datestr(now));

%% Load and parse the file
% =====================================================================
fdat = fopen(fn,'r');
tmp_hdr_fn = ct_filename_ct_tmp(rmfield(param,'day_seg'),'','headers', ...
        fullfile(adc_folder_name, [fn_name '.bin']));
tmp_hdr_fn_dir = fileparts(tmp_hdr_fn);
if ~exist(tmp_hdr_fn_dir,'dir')
  mkdir(tmp_hdr_fn_dir);
end
fhdr = fopen(tmp_hdr_fn,'w');
nrec = hfrds.strip_hdr(fdat,fhdr);
fclose(fhdr);
fhdr = fopen(tmp_hdr_fn,'r');
[hdr,wave] = hfrds.fread_data(fhdr,fdat,0,nrec);

utc_time_sod = 3600*(10*(hdr.gps(:,1)-'0')+(hdr.gps(:,2)-'0')) ...
  + 60*(10*(hdr.gps(:,3)-'0')+(hdr.gps(:,4)-'0')) ...
  + (10*(hdr.gps(:,5)-'0')+(hdr.gps(:,6)-'0'));

% Segment specific corrections to GPS/UTC time stamp in the radar header
if strcmpi(param.day_seg,'20160413_01')
  utc_time_sod(1:6) = utc_time_sod(1:6) + 1;
  gps_notes = 'Added 1 sec to utc_time_sod(1:6)\n';
elseif strcmpi(param.day_seg,'20160413_02')
  utc_time_sod(23813:23832) = utc_time_sod(23813:23832) + 60;
  gps_notes = 'Added 60 secs to utc_time_sod(23813:23832)\n';
else
  gps_notes = '';
end

records = [];
records.raw.epri = hdr.epri(:).';
records.raw.seconds = utc_time_sod(:).';
records.raw.fraction = hdr.frac(:).';

records.settings = [];

records.settings.wfs_records = 1;
records.settings.wfs.presums = hdr.pre(1) + 1;
records.settings.wfs.bit_shifts = hdr.bshft(1);
records.settings.wfs.num_sam = hdr.stop(1)-hdr.start(1)-1;

records.settings.wfs.t0 = (hdr.start(1) - hdr.delay(1)) / param.radar.fs;
if records.settings.wfs.t0 ~= -2.108e-5
  % Not the expected setting... need to actually look into getting t0 set
  % properly
  keyboard;
end
records.settings.wfs.t0 = 0;

% param.radar.prf = 2*param.radar.fs / (hdr.pri(1)+1);  %% NOT SURE
PRI = (hdr.pre(1)+1) / param.radar.prf; % hdr.pre = presums
radar_time = hdr.epri * PRI;

utc_time_sod_epri = radar_time - radar_time(1) ...
  + utc_time_sod(1) + hdr.frac(1)/param.radar.fs;

utc_time_sod_corrected = utc_time_sod_epri(:).';

%% Debug check
if 1
  figure(1); clf;
  plot(utc_time_sod+hdr.frac/param.radar.fs)
  hold on;
  plot(utc_time_sod_epri);
  xlabel('Record');
  ylabel('Time (sec)');
  title('UTC Time Seconds of Day');
  legend('Header Time','Corrected Time');
  figure(2); clf;
  plot(utc_time_sod+hdr.frac/param.radar.fs - utc_time_sod_epri)
  xlabel('Record');
  ylabel('Time Error (sec)');
  title('Header Time Minus Corrected Time');
  tmp_fn = ct_filename_ct_tmp(param,'','records','UTC_time.fig');
  tmp_fn_dir = fileparts(tmp_fn);
  if ~exist(tmp_fn_dir,'dir')
    mkdir(tmp_fn_dir);
  end
  saveas(1,tmp_fn);
  tmp_fn = ct_filename_ct_tmp(param,'','records','UTC_time_correction.fig');
  saveas(2,tmp_fn);
end

%% Correlate GPS with radar data
% ===================================================================
fprintf('Loading GPS data (%s)\n', datestr(now));

if param.records.gps.en
  records = sync_radar_to_gps(param,records,utc_time_sod_corrected);
  
else
  records.lat = NaN*zeros(size(utc_time_sod_corrected));
  records.lon = NaN*zeros(size(utc_time_sod_corrected));
  records.elev = NaN*zeros(size(utc_time_sod_corrected));
  records.gps_time = NaN*zeros(size(utc_time_sod_corrected));
  records.roll = NaN*zeros(size(utc_time_sod_corrected));
  records.pitch = NaN*zeros(size(utc_time_sod_corrected));
  records.heading = NaN*zeros(size(utc_time_sod_corrected));
  records.gps_source = 'NA';
end

%% Save records files
% =====================================================================

% Save concatenated files in records directories after time fix
records_fn = ct_filename_support(param,'','records');
[records_fn_dir records_fn_name] = fileparts(records_fn);
if ~exist(records_fn_dir,'dir')
  fprintf('Output directory %s does not exist, creating...\n', records_fn_dir);
  mkdir(records_fn_dir);
end

% Standard Fields
% records
%  .lat
%  .lon
%  .elev
%  .roll
%  .pitch
%  .heading
%  .gps_time
%  .gps_source
records.surface = NaN*zeros(size(records.lat));
records.relative_filename = {{[fn_name fn_ext]}};
records.relative_rec_num = {1};
records.offset = hdr.loc(:).';

records.radar_name = param.radar_name;
records.ver = 3;

records.notes = cat(2,sprintf('\nGPS NOTES\n%s',gps_notes));

records.param_records = param;

fprintf('Saving records file %s (%s)\n',records_fn,datestr(now));
save(records_fn,'-v7.3','-struct','records'); % Handle large file sizes, so use v7.3

% =====================================================================
% Create record aux files for faster loading times
% =====================================================================

fprintf('Creating auxiliary records files %s (%s)\n',records_fn,datestr(now));
create_records_aux_files(records_fn);

fprintf('Done (%s)\n\n', datestr(now));

return;
