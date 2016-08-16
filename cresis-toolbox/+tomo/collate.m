function collate(param, param_override)
% tomo.collate.m
%
% Description. Usually this function is called from tomo.run_collate.
%   Calls data_loader_prep, DEM_alignment, and surface_extractor.
%
% Inputs:
%   param = struct with processing parameters
%   param_override = parameters in this struct will override parameters
%     in param.
%
% See also: tomo.run_collate, tomo.data_loader_prep, tomo.DEM_alignment,
%   tomo.surface_extractor
%
% Author: John Paden, Jordan Sprick, and Mingze Xu

param = merge_structs(param,param_override);

% fn_dir: Directory where 3D image files are at
fn_dir = ct_filename_out(param,param.surf_extract.out_dir);

if ~isfield(param.records,'records_fn')
  param.records.records_fn = '';
end
if ~isfield(param.records,'frames_fn')
  param.records.frames_fn = '';
end

% Load frames file
load(ct_filename_support(param,param.records.frames_fn,'frames'));
% Load records file
records_fn = ct_filename_support(param,param.records.records_fn,'records');
records = load(records_fn);

if isempty(param.cmd.frms)
  param.cmd.frms = 1:length(frames.frame_idxs);
end
% Remove frames that do not exist from param.cmd.frms list
[valid_frms,keep_idxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
if length(valid_frms) ~= length(param.cmd.frms)
  bad_mask = ones(size(param.cmd.frms));
  bad_mask(keep_idxs) = 0;
  warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
    param.cmd.frms(find(bad_mask,1)));
  param.cmd.frms = valid_frms;
end

if 0
  % Compile C++ functions
  % If you get a C++11 option error, you may be using pre-G++ 4.7. You can
  % check the g++ version with system('g++ --version');
  % To fix this, add -v option to mex function and look for a line like this:
  %   Options file: ~/.matlab/R2015b/mex_C++_glnxa64.xml
  % Replace -std=c++11 with -std=c++0x (should occur in two places)
  % Reference: http://stackoverflow.com/questions/14674597/cc1plus-error-unrecognized-command-line-option-std-c11-with-g
  mex -largeArrayDims fuse.cpp
  mex -largeArrayDims train_params.cpp
  mex -largeArrayDims detect.cpp
  mex -largeArrayDims extract_flag.cpp
end

for frm_idx = 1:length(param.cmd.frms)
  frm = param.cmd.frms(frm_idx);

  % Load Data
  fprintf('Loading frame data...\n');
  mdata = {};
  for img=1:3
    fn = fullfile(fn_dir,sprintf('Data_img_%02.0f_%s_%03.0f.mat',img, ...
      param.day_seg,frm));
    mdata{img} = load(fn);
    mdata{img}.frm = frm;
  end
  
  if param.surf_extract.add_layers_flag
    mdata = tomo.data_loader_prep(param,mdata);
  end
  
  if param.surf_extract.ice_twtt_flag
    mdata = tomo.DEM_alignment(param,mdata);
  end
    
  if param.surf_extract.extract_flag
    mdata_combined = tomo.surface_extractor(param,mdata);
  end
  
end