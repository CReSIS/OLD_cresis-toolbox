%% Setup

params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),'','post');
% params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),'20140401_03','post');
% params.cmd.generic = 1;
% params.cmd.frms = 22;
ct_set_params(params,'cmd.generic',0);
ct_set_params(params,'cmd.generic',1,'20140401_03|20140506_01|20140325_05|20140325_06|20140325_07');
param_override.update.input = 'CSA_music_surfData_no_QC';
param_override.update.output = 'surfData_v2_no_MC';
param_override.update.echogram = 'CSA_music';

% params = read_param_xls(ct_filename_param('rds_param_2009_Antarctica_TO.xls'),'','post');
% params.cmd.generic = 1;
% params.cmd.frms = [];
% param_override.update.input = 'surfData';
% param_override.update.output = 'surfData_v2';
% param_override.update.echogram = 'music3D';

%% Automated section
global gRadar;

% Input checking
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  param = merge_structs(param, param_override);

  
  % Load frames file
  load(ct_filename_support(param,'','frames'));
  
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

  for frm = param.cmd.frms
    
    fn = fullfile(ct_filename_out(param,param.update.input,''),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
    fn_cur_ver = fullfile(ct_filename_out(param,param.update.output,''),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
    echogram_fn = fullfile(ct_filename_out(param,param.update.echogram,''),sprintf('Data_%s_%03d.mat',param.day_seg,frm));

    fprintf('Convert\n  %s\n  %s\n  %s\n', fn, fn_cur_ver,echogram_fn);
    tomo.surfdata.update_file(fn,fn_cur_ver,echogram_fn);
  end
  
end

