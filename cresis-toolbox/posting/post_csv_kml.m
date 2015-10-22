% script post_csv_kml
%
% Creates a CSV and KML file for each segment file.
%
% Outputs:
%   gRadar.out_path/param.radar_name/param.season_name/CSARP_post/csv/
%     Data_YYYYMMDD_SS.csv
%     Browse_Data_YYYYMMDD_SS.kml
%
% Author: John Paden

% ======================================================================
% User Settings
% ======================================================================

% Parameters spreadsheet to use for updating
%   1. Segment and frame list are taken from the parameter sheet
%   2. For GPS update, GPS time offsets are pulled from the parameter sheet
param_fn = '/users/paden/scripts/branch/params-cr1/kuband_param_2012_Antarctica_DC8.xls';

% param.post_dir: occasionally you may want to operate in the post directory
%   set this to CSARP_post, otherwise leave empty
post_dir = '';

% param.skip_phrase: All segments will be skipped with this phrase in their
%   verification field. "do not process" is the standard. Leave this field
%   blank to do all segments.
skip_phrase = 'do not process';

% point_spacing: specifies the spacing between geographic points in the CSV
% and KML files
point_spacing = 1000;

% ======================================================================
% Automated Section
% ======================================================================

fprintf('=============================================================\n');
fprintf('post_csv_kml (%s)\n\n', datestr(now));

% params: structure of parameters for each segment
params = read_param_xls(param_fn);

for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isempty(skip_phrase) ...
      && ~isempty(strfind(lower(param.cmd.notes),skip_phrase)) ...
      || ~param.cmd.generic
    continue;
  end
  
  csv_dir = fullfile(ct_filename_out(param, ...
    post_dir, 'CSARP_post', true),'csv');
  if ~exist(csv_dir,'dir')
    mkdir(csv_dir)
  end
  kml_dir = fullfile(ct_filename_out(param, ...
    post_dir, 'CSARP_post', true),'kml');
  if ~exist(kml_dir,'dir')
    mkdir(kml_dir)
  end
  
  records_fn = ct_filename_support(param, '', 'records');
  fprintf('Loading %s\n', records_fn);
  records = load(records_fn);
  UTC_sod = epoch_to_sod(records.gps_time,param.day_seg(1:8));
  
  frames_fn = ct_filename_support(param, '', 'frames');
  fprintf('Loading %s\n', frames_fn);
  load(frames_fn);
  
  % Create new filename
  csv_fn = sprintf('Browse_Data_%s.csv',param.day_seg);
  csv_fn = fullfile(csv_dir,csv_fn);
  
  fid_csv = fopen(csv_fn,'wt');
  fprintf(fid_csv,'%s,%s,%s,%s,%s\n',...
    'LAT','LON','ELEVATION','UTCTIMESOD','FRAME');
  % Make a point every point_spacing and include at least the first and last point in each frame
  along_track = geodetic_to_along_track(records.lat,records.lon,records.elev);
  for frm = 1:length(frames.frame_idxs)
    first_rec = frames.frame_idxs(frm);
    if frm < length(frames.frame_idxs)
      last_rec = frames.frame_idxs(frm+1)-1;
    else
      last_rec = length(records.lat);
    end
    
    recs = first_rec;
    rec = first_rec + 1;
    while rec < last_rec
      if along_track(rec) > along_track(recs(end)) + point_spacing
        recs(end+1) = rec;
      end
      rec = rec + 1;
    end
    recs = [recs last_rec];
    
    frm_id_num = [param.day_seg(1:8) param.day_seg(10:11) sprintf('%03d',frm)];
    for rec = recs
      fprintf(fid_csv,'%2.6f,%2.6f,%5.4f,%4.4f,%s\n',...
        records.lat(rec),records.lon(rec),...
        records.elev(rec),UTC_sod(rec),frm_id_num);
    end
  end
  fclose(fid_csv);
  fprintf('  Writing %s\n',csv_fn);
  
  % Create KML browse files for each segment
  % Extract day_seg from filename
  in_fn = csv_fn;
  [in_fn_dir in_fn_name] = fileparts(in_fn);
  kml_out_fn = fullfile(kml_dir, [in_fn_name '.kml']);
  day_seg = in_fn_name(13:end);
  kml_write_cresis(in_fn, kml_out_fn, day_seg,'vectors',1);
  fprintf('  Writing %s\n',kml_out_fn);
end

% Create a concatenated CSV browse file of the whole season
in_search_str = fullfile(csv_dir,'*.csv');
delim = ',';
season_fn_name = sprintf('Browse_%s',param.season_name);
out_fn = fullfile(csv_dir, [season_fn_name '.csv']);
fprintf('Creating concatenated CSV browse file of the whole season %s\n',out_fn);
sys_cmd = sprintf('cat %s  | grep -v LAT | sed 1i"LAT%sLON%sELEVATION%sUTCTIMESOD%sFRAME" >%s', ...
  in_search_str,delim,delim,delim,delim,out_fn);
system(sys_cmd);

return;

