% script update_surface_twtt_delta
%
% Updates the Surface variable with a delay offset.
% Only works on single waveform/single adc files.

%% User Settings
% ----------------------------------------------------------------------
params = read_param_xls(ct_filename_param('kuband_param_2014_Greenland_P3.xls'));
data_type = 'qlook';

%% Automated Section
% ----------------------------------------------------------------------
for param_idx = 1:length(params)
  param = params(param_idx);
  
  if ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  fprintf('Updating surface %s (%s)\n', param.day_seg, datestr(now,'HH:MM:SS'));
  
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
  
  for frm_idx = 1:length(param.cmd.frms)
    frm = param.cmd.frms(frm_idx);
    
    data_fn = fullfile(ct_filename_out(param,data_type,''), ...
      sprintf('Data_%s_%03d.mat', param.day_seg, frm));
    
    warning off;
    mdata = load(data_fn,'Time','Surface','Bottom','param_get_heights','param_csarp');
    warning on;
    
    fields_to_update = {'Time'};
    
    % Should be only one of these types of fields:
    if isfield(mdata,'param_get_heights')
      delta_offset = param.radar.wfs(1).Tsys - mdata.param_get_heights.radar.wfs(1).Tsys;
      mdata.param_get_heights.radar.wfs(1).Tsys = param.radar.wfs(1).Tsys;
      fields_to_update{end+1} = 'param_get_heights';
    elseif isfield(mdata,'param_csarp')
      delta_offset = param.radar.wfs(1).Tsys - mdata.param_get_heights.radar.wfs(1).Tsys;
      mdata.param_csarp.radar.wfs(1).Tsys = param.radar.wfs(1).Tsys;
      fields_to_update{end+1} = 'param_csarp';
    end
    
    if delta_offset ~= 0
      fprintf('  Offset %g %s (%s)\n', delta_offset, data_fn, datestr(now,'HH:MM:SS'));
      
      
      mdata.Time = mdata.Time - delta_offset;
      if isfield(mdata,'Surface')
        mdata.Surface = mdata.Surface - delta_offset;
        fields_to_update{end+1} = 'Surface';
      end
      if isfield(mdata,'Bottom')
        mdata.Bottom = mdata.Bottom - delta_offset;
        fields_to_update{end+1} = 'Bottom';
      end
      save(data_fn,'-append','-struct','mdata',fields_to_update{:});
      
    end
  end
  
end

return;



%% Test Code
% -------------------------------------------------------------------------

sd = load('/scratch/snow/2014_Greenland_P3/CSARP_qlook/20140516_02/Data_20140516_02_035.mat');
kd = load('/scratch/kuband/2014_Greenland_P3/CSARP_qlook/20140516_02/Data_20140516_02_035.mat');

% sd = load('/scratch/snow/2014_Greenland_P3/CSARP_qlook/20140516_02/Data_20140516_02_300.mat');
% kd = load('/scratch/kuband/2014_Greenland_P3/CSARP_qlook/20140516_02/Data_20140516_02_300.mat');

surface_correction = 6.1170e-09
surface_correction = 0
kd.Time = kd.Time-surface_correction;
kd.Surface = kd.Surface-surface_correction;

figure(1); clf;
imagesc([],sd.Time,lp(sd.Data));
hold on;
plot(sd.Surface);

figure(2); clf;
imagesc([],kd.Time,lp(kd.Data))
hold on;
plot(kd.Surface);

figure(3); clf;
plot(sd.Surface - kd.Surface);
surface_correction = surface_correction-mean(sd.Surface - kd.Surface)



