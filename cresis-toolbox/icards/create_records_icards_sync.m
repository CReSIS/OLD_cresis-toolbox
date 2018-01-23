% create_records_icards_sync
% used templated fmcw_accum_sync
%
% Script called from create_records_fmcw_accum
%
% To run this in debug mode, you need to set the debug setup section in
% this file and then just run this as a script (must be run from ">>"
% prompt and not a debug "K>>" prompt).
%
% Author: John Paden

dbstack_info = dbstack;
global gRadar;%QQQQ
if ~exist('param','var') || isempty(param) || length(dbstack_info) == 1
  % =====================================================================
  % Debug Setup
  % =====================================================================
  
  new_param = read_param_xls(ct_filename_param('rds_param_2002_Greenland_P3.xls'),'20020520_01');
  new_param.vectors.file.base_dir='Z:\ICARDS\2002\';
  new_param.date=new_param.day_seg(1:8);
  fn = ct_filename_ct_tmp(new_param,'','records','workspace');
  fn = [fn '.mat'];
  fprintf('Loading workspace %s (%s)\n', fn, datestr(now));
  
  if exist(fn,'file')
    load(fn);
  else
    error('Temporary records file does not exist');
  end
  param = new_param;
  clear new_param;
  clear('param_override');
  param_override.sched.type = 'no scheduler';
  
  % Input checking
  if ~exist('param','var')
    error('A struct array of parameters must be passed in\n');
  end
  if exist('param_override','var')
    param_override = merge_structs(gRadar,param_override);
  else
    param_override = gRadar;
  end
  param = merge_structs(param, param_override);
end

fprintf('Running %s correction and gps sync (%s)\n', param.day_seg, datestr(now));

% ======================================================================
%See if need gps_time_offset  
if isempty(param.vectors.gps.time_offset)
  param.vectors.gps.time_offset=0;
end 

if ~isfield(param.records,'manual_time_correct') || isempty(param.records.manual_time_correct)
  param.records.manual_time_correct = 0;
end
% % % % % % % radar_time_notes = '';
%% Create hdr.utc_time_sod vector
% NMEA TIME=UTC TIME?
% hdr.utc_time_sod = double(hdr.seconds) + double(hdr.fraction)/param.radar.fs;
hdr.radar_time = hdr.nmea_time + param.vectors.gps.time_offset;
records=[];

%% Find EPRI (effective pulse repetition interval) by radar name
%NO EPRI IN ICARDS

%% Check to see if 1 PPS signal was connected?
%NEED EPRI TO CHECK?

% =====================================================================
%% Find bad UTC time SOD and EPRI entries
% =====================================================================
% % % % % % % % orig_hdr = hdr;
% % % % % % % % epri_notes = [];

%% Quick Load using just first and last record


%% Correct for digital errors in epri
if param.records.gps.en        
   utc_time_sod = epoch_to_sod(hdr.radar_time,param.day_seg(1:8));
   utc_time_sod = utc_time_sod + param.vectors.gps.time_offset;
   records = sync_radar_to_gps(param,records,utc_time_sod);
else
  records.lat = NaN*zeros(size(hdr.utc_time_sod));
  records.lon = NaN*zeros(size(hdr.utc_time_sod));
  records.elev = NaN*zeros(size(hdr.utc_time_sod));
  records.gps_time = NaN*zeros(size(hdr.utc_time_sod));
  records.roll = NaN*zeros(size(hdr.utc_time_sod));
  records.pitch = NaN*zeros(size(hdr.utc_time_sod));
  records.heading = NaN*zeros(size(hdr.utc_time_sod));
  records.gps_source = 'NA';
end

% =====================================================================
% Save records files
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
records.relative_filename{1}=hdr.filenames;
records.relative_rec_num{1}=hdr.file_rec_offset;
%  .relative_filename{board_idx == 1}{1 .. Nf(1)}
%  .relative_rec_num{board_idx == 1}(1 .. Nf(1))
%OFFSET: Nb by Nx uint32 array of file byte offset for each record (rows select board/file-groups, columns select record)
%%records.offset(1,:) = hdr.offset{1}(:);
%Offset in records is the position in the file for each data record.
%i.e. fseek(fid,records.offset{1}(N),-1)  would take you to the position in the file where record N starts
% % % % % % % records.offset=hdr.offset;
records.radar_name = param.radar_name;
records.ver = 3;
%-----------add settings field---------------------------------------------
records.settings = [];
records.settings.wfs_records = 1;
records.settings.wfs_file=1;
wf=1;
records.settings.wfs(wf).num_sam = num_rec_sample; 
records.settings.wfs(wf).t0 = 0;
records.settings.wfs(wf).presums = 1;
records.settings.wfs(wf).bit_shifts = 0;
records.settings.wfs(wf).fs = 1.88e7;                        %get from parameter spread sheet
records.settings.wfs(wf).f0 = 1.8e8;                         %get from parameter spread sheet
records.settings.wfs(wf).f1 = 2.1e8;                         %get from parameter spread sheet
records.settings.wfs(wf).Tpd=3e-6;                           %get from parameter spread sheet
records.settings.wfs(wf).adc_gains=10.^((52-0*ones(1,7))/20);%get from parameter spread sheet
%--------------------------------------------------------------------------
records.notes = '';
records.param_records = param;
   
% % % % if length(hdr.offset)~=records.gps_time
% % % %     records.offset=cat(2,hdr.offset,hdr.offset(end)+[1:1:(length(records.gps_time)-length(hdr.offset))]*hdr.offset(1));
% % % % else
% % % % records.offset=hdr.offset;
% % % % end
records.offset=hdr.offset;

fprintf('Saving records file %s (%s)\n',records_fn,datestr(now));
save(records_fn,'-v6','-struct','records');

% =====================================================================
% Create record aux files for faster loading times
% =====================================================================

fprintf('Creating auxiliary records files %s (%s)\n',records_fn,datestr(now));
create_records_aux_files(records_fn);


fprintf('Done (%s)\n\n', datestr(now));

