function success = tomo_collate_task(param)
% success = tomo_collate_task(param)
%
% Description. Usually this function is called from tomo.collate.
%   Calls data_loader_prep, DEM_alignment, and surface_extractor.
%
% Inputs:
%   param: struct with processing parameters
%    .tomo_collate
%     .fuse_images_flag
%     .create_surfData_flag
%     .create_surfData
%   param_override: parameters in this struct will override parameters
%     in param.
%
% See also: tomo.run_collate, tomo.collate, tomo_collate_task,
%   tomo.fuse_images, tomo.add_icemask_surfacedem, tomo.create_surfdata,
%
% Author: John Paden, Jordan Sprick, and Mingze Xu

%% Combine beams into a single data file
if param.tomo_collate.fuse_images_flag
  fprintf('Fusing Images (%s)\n', datestr(now));
  mdata = tomo.fuse_images(param);
  fprintf('  Done (%s)\n', datestr(now));
end

%% Add ice mask and ice surface DEM into data file
if param.tomo_collate.add_icemask_surfacedem_flag
  fprintf('Adding ice mask and ice surface DEM to data file (%s)\n', datestr(now));
  
  if ~exist('mdata','var')
    % in_dir: Directory where 3D image files are at
    in_dir = ct_filename_out(param,param.tomo_collate.in_dir);
    
    % combined_fn: Filename with 3D data
    combined_fn = fullfile(in_dir,sprintf('Data_%s_%03.0f.mat',param.day_seg,param.proc.frm));
    mdata = load(combined_fn);
  end
  
  mdata = tomo.add_icemask_surfacedem(param, mdata);
  
  fprintf('  Done (%s)\n', datestr(now));
end

%% Create surfData file
if param.tomo_collate.create_surfData_flag
  fprintf('Creating surfData file (%s)\n', datestr(now));
  
  if ~exist('mdata','var')
    % in_dir: Directory where 3D image files are at
    in_dir = ct_filename_out(param,param.tomo_collate.in_dir);
    
    % combined_fn: Filename with 3D data
    combined_fn = fullfile(in_dir,sprintf('Data_%s_%03.0f.mat',param.day_seg,param.proc.frm));
    fprintf('Loading Data\n');
    mdata = load(combined_fn);
  end
  
  tomo.create_surfdata(param,mdata);
  fprintf('  Done (%s)\n', datestr(now));
end

success = true;

end
