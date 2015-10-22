% Function for copying backups in the field from the scratch drive to the
% NAS (or some other archive).

param = [];
param.season_name = '2013_Greenland_P3';
radar_names = {'accum'};

backup_dir = '/mnt/products/products/';

copy_support_files = false;

copy_type = 1; % 1=rsync, 0=md5_copy

qc_level = -1; % md5_copy: -1=no md5, 0=md5 on copied files only, 1=md5 on checks/copies
rsync_options = '-rtv'; % rsync: command line options for rsync command
% rsync_options = '-av'; % rsync: command line options for rsync command

% ========================================================================
% Automated Section
% ========================================================================

param = merge_structs(param, gRadar);

% ========================================================================
%% Copy output data products
% ========================================================================
for radar_idx = 1:length(radar_names)
  param.radar_name = radar_names{radar_idx};
  
  out_path = ct_filename_out(param,'','CSARP_post',1);
  
  backup_param = param;
  backup_param.out_path = backup_dir;
  archive_path = ct_filename_out(backup_param,'','CSARP_post',1);
  
  if copy_type
    mkdir(archive_path);
    cmd = sprintf('rsync %s %s %s',rsync_options, [out_path filesep],archive_path)
    system(cmd);
  else
    md5_copy(out_path,archive_path,qc_level)
  end
  
  %   if any(strcmpi(ct_output_dir(param.radar_name),{'kuband','snow'}))
  %     out_path = ct_filename_out(param,'','CSARP_qlook',1);
  %
  %     archive_path = ct_filename_out(backup_param,'','CSARP_qlook',1);
  %
  %     if copy_type
  %       mkdir(archive_path);
  %       cmd = sprintf('rsync %s %s %s',rsync_options, [out_path filesep],archive_path)
  %       system(cmd);
  %     else
  %       md5_copy(out_path,archive_path,qc_level)
  %     end
  %   end
  
end

% ========================================================================
%% Copy support files
% ========================================================================

if copy_support_files
  % ======================================================================
  % Backup csarp-support
  % ======================================================================
  [out_dir out_fn] = fileparts(gRadar.support_path);
  if isempty(out_fn)
    [out_dir out_fn] = fileparts(out_dir);
  end
  archive_path = fullfile(backup_dir,out_fn);
  
  if copy_type
    mkdir(archive_path);
    cmd = sprintf('rsync %s %s %s',rsync_options, [gRadar.support_path filesep],archive_path)
    system(cmd);
  else
    md5_copy(gRadar.support_path,archive_path,qc_level)
  end
  
  % ======================================================================
  % Backup GIS_data
  % ======================================================================
  [out_dir out_fn] = fileparts(gRadar.gis_path);
  if isempty(out_fn)
    [out_dir out_fn] = fileparts(out_dir);
  end
  archive_path = fullfile(backup_dir,out_fn);
  
  if copy_type
    mkdir(archive_path);
    cmd = sprintf('rsync %s %s %s',rsync_options, [gRadar.gis_path filesep],archive_path)
    system(cmd);
  else
    md5_copy(gRadar.gis_path,archive_path,qc_level)
  end
  
  % ======================================================================
  % Backup metadata
  % ======================================================================
  [out_dir out_fn] = fileparts(gRadar.data_support_path);
  if isempty(out_fn)
    [out_dir out_fn] = fileparts(out_dir);
  end
  archive_path = fullfile(backup_dir,out_fn);
  
  if copy_type
    mkdir(archive_path);
    cmd = sprintf('rsync %s %s %s',rsync_options, [gRadar.data_support_path filesep],archive_path)
    system(cmd);
  else
    md5_copy(gRadar.data_support_path,archive_path,qc_level)
  end
end

return;


