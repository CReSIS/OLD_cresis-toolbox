% script nsidc_data_manifest
%
% Creates the legacy data manifest for OIB data. This function is defunct and
% should no longer be used.
%
% See also: type "nsidc_help.m"
%
% Author: John Paden

error('This function is defunct and should no longer be used.');

% =======================================================================
% User Settings
% =======================================================================

in_dir = 'L:\CSARP_post';

out_fn_dir = 'L:\';

% .release = Release 1.0 (field or level 1B) or 2.0 (post processed or level 1B/2)
param.release = 1.0;

param.data_type = 'Level 1B';
% param.data_type = 'Level 2';
% param.data_type = 'Level 1B/Level 2';
% param.data_type = 'Level 3';
% param.data_type = 'GPS/INS';

% param.radar = 'accum';
% param.radar = 'gps_ins';
param.radar = 'kuband';
% param.radar = 'icards';
% param.radar = 'acords';
% param.radar = 'mcrds';
% param.radar = 'mcords';
% param.radar = 'snow';
param.season_name = '2012_Greenland_P3';

% .transfer_type = Method Of Data Transfer (Ftp, HardDrive, etc)
param.transfer_type = 'HardDrive';

% =======================================================================
% Automated Section
% =======================================================================
nsidc_data_manifest_tstart = tic;
radar = param.radar;

release_str = sprintf('%03.1f',param.release);
release_str(2) = '_';
data_type_str = param.data_type(param.data_type ~= ' ' & param.data_type ~= '/');
out_fn = fullfile(out_fn_dir, sprintf('%s_%s_%s_release_v%s_manifest.txt', param.radar, data_type_str, param.season_name, release_str));
[fid msg] = fopen(out_fn, 'w');
if fid < 1
  fprintf('Could not open file %s\n', out_fn);
  error(msg);
end

fprintf('Opened %s manifest file\n', out_fn);

fprintf(fid, 'Organization Name: CReSIS\n');
fprintf(fid, 'Contact Email Address: cresis_data@cresis.ku.edu\n');
fprintf(fid, '\n');
if param.release == 1
  release_str = sprintf('%s release %.1f', param.data_type, param.release);
elseif param.release > 1
  release_str = sprintf('%s release %.1f', param.data_type, param.release);
end
[year remain] = strtok(param.season_name,'_');
year = str2double(year);
[location remain] = strtok(remain,'_');
platform = strtok(remain,'_');

fprintf('Date Sent: %s\n', datestr(now));
%fprintf('Data Sent: %s %04d %s %s %s\n', radar, year, location, platform, release_str);
fprintf('Data Set Title: %s %04d %s %s %s\n', radar, year, location, platform, release_str);
fprintf('Data Set ShortName (if known):\n');
fprintf('Campaign(Region/year): %s/%04d\n', location, year);

fprintf(fid, 'Date Sent: %s\n', datestr(now));
%fprintf(fid, 'Data Sent: %s %04d %s %s %s\n', radar, year, location, platform, release_str);
fprintf(fid, 'Data Set Title: %s %04d %s %s %s\n', radar, year, location, platform, release_str);
fprintf(fid,'Data Set ShortName (if known):\n');
fprintf(fid, 'Campaign(Region/year): %s/%04d\n', location, year);

if param.release == 1
  fprintf(fid, 'Is this the first delivery of this Data Set?: Yes\n');
elseif param.release > 1
  fprintf(fid, 'Is this the first delivery of this Data Set?: No\n');
end
fprintf(fid, 'Has the data format changed since the most recent delivery?: No\n');
fprintf(fid, 'Has the file naming convention changed since the most recent delivery?: No\n');
fprintf(fid, 'If users will require an updated reader to view these files, is it included in the delivery?: No\n');

% fprintf(fid, 'Date Sent: MCoRDS 2009 Antarctica DC8 Level 1B/2 release 2.0\n');
% fprintf(fid, 'Data Set Title: MCoRDS 2009 Antarctica DC8 Level 1B/2 release 2.0\n');
% fprintf(fid, 'Campaign(Region/year): Antarctica/2009\n');
fprintf(fid, '\n');
fprintf(fid, 'Method Of Data Transfer (Ftp, HardDrive, etc): %s\n', param.transfer_type);

fprintf('Getting file list (%.1f sec)\n', toc(nsidc_data_manifest_tstart));

% Get the file list and the file sizes simultaneously
if ispc
  dirs = get_filenames(in_dir,'','','',struct('type','d','recursive',1));
  fns = {};
  file_size = [];
  for dir_idx = 1:length(dirs)
    dir_fns = dir(dirs{dir_idx});
    for file_idx = 1:length(dir_fns)
      if ~dir_fns(file_idx).isdir
        fns{end+1} = fullfile(dirs{dir_idx},dir_fns(file_idx).name);
        file_size(end+1) = dir_fns(file_idx).bytes;
      end
    end
  end
else
  sys_cmd = sprintf('du -ab %s', in_dir)
  [tmp du_output] = system(sys_cmd);
  % Parse the output
  C = textscan(du_output,'%f%s');
  % Remove the directories from the output, keeping just the files
  is_file_mask = ones(size(C{2}));
  for fn_idx = 1:length(C{2})
    if exist(C{2}{fn_idx},'dir')
      is_file_mask(fn_idx) = 0;
    end
  end
  fns = C{2}(logical(is_file_mask));
  file_size = C{1}(logical(is_file_mask));
end

% Extract extensions
fns_ext = fns;
for fn_idx = 1:length(fns)
  [fn_dir fn_name fn_ext] = fileparts(fns{fn_idx});
  fns_ext{fn_idx} = fn_ext;
end
fns_ext = unique(fns_ext);
fprintf(fid, 'Data Type(s) (.mat, .pdf, .jpeg, etc): ');
for ext_idx = 1:length(fns_ext)
  if ext_idx == 1
    fprintf(fid, '%s', fns_ext{ext_idx});
  else
    fprintf(fid, ', %s', fns_ext{ext_idx});
  end
end
fprintf(fid, '\n');

fprintf(fid, 'Total Number Of Files: %d\n', length(fns));
fprintf(fid, 'Total Size Of Data (mb): %.3f Megabytes\n', sum(file_size)/2^20);
fprintf(fid, 'Type Of Checksum (MD5 or CKSUM): MD5\n');
fprintf(fid, '\n');

fprintf('Writing MD5Sums to manifest file (%.1f sec)\n', toc(nsidc_data_manifest_tstart));
for fn_idx = 1:length(fns)
  fprintf('  File %d of %d (%.1f sec)\n', fn_idx, length(fns), toc(nsidc_data_manifest_tstart));

  fn = fns{fn_idx};
  if ispc
    [tmp md5sum] = system(sprintf('md5sums -u %s',fns{fn_idx}));
    md5sum = md5sum(1:32);
  else
    [tmp md5sum] = system(sprintf('md5sum -b %s | cut -b 1-32', ...
      fns{fn_idx}));
  end
  if tmp ~= 0
    error('MD5 sum failed');
  end
  
  fprintf(fid, '%s %d %s\n', fn(length(in_dir)+1:end), file_size(fn_idx), ...
    md5sum);
end

fclose(fid);
fprintf('Done\n');

return;
