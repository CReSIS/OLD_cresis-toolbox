function update_surface_twtt_delta(param,param_override)
% update_surface_twtt_delta(param,param_override)
%
% Updates the fast time related variables for changes in Tadc_adjust and Tsys.
% Can be used on any echogram data product (not CSARP_out). However, for
% SAR processed outputs, the data should be reprocessed if the offset is
% very large because large time system time delay errors affect the focussing.
%
% Updates Time, param_{qlook,sar}.radar.wfs(wf).Tadc_adjust,
% param_{qlook,sar}radar.wfs(wf).Tsys, Surface, and Bottom variables.
%
% Note that only the param_{qlook,sar} field that matters gets
% updated. For example, param_records and param_combine will not be used or
% updated because Tadc_adjust and Tsys are not used during these processes
% and so their value does not matter when those processes are applied.
%
% Tsys, t_ref and Tadc_adjust changes can only be corrected when the same
% change is applied to all waveform-adc pairs used in the data product.
%
% Note that Tadc_adjust and t_ref are opposite sign to Tsys
%  - Tsys is subtracted away from time
%  - Tadc_adjust and t_ref are added on to time
%
% Examples:
%   See run_update_surface_twtt_delta.m to run
%
% Author: John Paden

%% Setup

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

% Remove frames that do not exist from param.cmd.frms list
load(ct_filename_support(param,'','frames'));
if isempty(param.cmd.frms)
  param.cmd.frms = 1:length(frames.frame_idxs);
end
[valid_frms,keep_idxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
if length(valid_frms) ~= length(param.cmd.frms)
  bad_mask = ones(size(param.cmd.frms));
  bad_mask(keep_idxs) = 0;
  warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
    param.cmd.frms(find(bad_mask,1)));
  param.cmd.frms = valid_frms;
end

% Rename a few variables for readability
data_types = param.update_surface_twtt_delta.data_types;
imgs = param.update_surface_twtt_delta.imgs;

%% Loop through each frame and update
for frm_idx = 1:length(param.cmd.frms)
  frm = param.cmd.frms(frm_idx);
  fprintf('%s_%03d\n', param.day_seg, frm);
  
  %% Loop through each data type
  for data_type_idx = 1:length(data_types)
    data_type = data_types{data_type_idx};
    
    %% Loop through each image type
    for data_img = imgs
      
      %% Load the data
      if data_img == 0
        data_fn = fullfile(ct_filename_out(param,data_type,''), ...
          sprintf('Data_%s_%03d.mat', param.day_seg, frm));
      else
        data_fn = fullfile(ct_filename_out(param,data_type,''), ...
          sprintf('Data_img_%02d_%s_%03d.mat', data_img, param.day_seg, frm));
      end
      
      if ~exist(data_fn,'file')
        warning('Missing %s\n', data_fn);
        continue;
      end
      
      warning off;
      mdata = load(data_fn,'Time','Surface','Bottom','param_qlook','param_sar');
      warning on;
      
      fields_to_update = {'Time'};
      
      if isfield(mdata,'param_qlook')
        param_field = 'param_qlook';
      elseif isfield(mdata,'param_sar')
        param_field = 'param_sar';
      else
        error('Data file contained neither param_qlook or param_sar, but should contain one or the other.');
      end
      fields_to_update{end+1} = param_field;
      
      %% Determine which fields changed
      delta_offset_t_ref = [];
      delta_offset_t_ref_mask = logical([]);
      delta_offset_Tsys = [];
      delta_offset_Tsys_mask = logical([]);
      delta_offset_Tadc_adjust = [];
      delta_offset_Tadc_adjust_mask = logical([]);
      for img = 1:length(mdata.(param_field).(param_field(7:end)).imgs)
        for wf_adc_pair = 1:size(mdata.(param_field).(param_field(7:end)).imgs{img},1)
          wf = abs(mdata.(param_field).(param_field(7:end)).imgs{img}(wf_adc_pair,1));
          adc = abs(mdata.(param_field).(param_field(7:end)).imgs{img}(wf_adc_pair,2));
          
          % t_ref correction
          if ~isfield(param.radar.wfs(wf),'t_ref') ...
              || isempty(param.radar.wfs(wf).t_ref)
            param.radar.wfs(wf).t_ref = 0;
          end
          if ~isfield(mdata.(param_field).radar.wfs(wf),'t_ref') ...
              || isempty(mdata.(param_field).radar.wfs(wf).t_ref)
            mdata.(param_field).radar.wfs(wf).t_ref = 0;
          end
          delta_offset_t_ref(wf) = param.radar.wfs(wf).t_ref ...
            - mdata.(param_field).radar.wfs(wf).t_ref;
          delta_offset_t_ref_mask(wf) = true;
          
          % Tsys correction
          if ~isfield(param.radar.wfs(wf),'Tsys') ...
              || numel(param.radar.wfs(wf).Tsys) < adc
            param.radar.wfs(wf).Tsys(adc) = 0;
          end
          if ~isfield(mdata.(param_field).radar.wfs(wf),'Tsys') ...
              || numel(mdata.(param_field).radar.wfs(wf).Tsys) < adc
            mdata.(param_field).radar.wfs(wf).Tsys(adc) = 0;
          end
          
          delta_offset_Tsys(img,wf_adc_pair) = param.radar.wfs(wf).Tsys(adc) ...
            - mdata.(param_field).radar.wfs(wf).Tsys(adc);
          delta_offset_Tsys_mask(img,wf_adc_pair) = true;
          
          % Tadc_adjust correction
          if ~isfield(param.radar.wfs(wf),'Tadc_adjust') ...
              || isempty(param.radar.wfs(wf).Tadc_adjust)
            param.radar.wfs(wf).Tadc_adjust = 0;
          end
          if ~isfield(mdata.(param_field).radar.wfs(wf),'Tadc_adjust') ...
              || isempty(mdata.(param_field).radar.wfs(wf).Tadc_adjust)
            mdata.(param_field).radar.wfs(wf).Tadc_adjust = 0;
          end
          delta_offset_Tadc_adjust(wf) = param.radar.wfs(wf).Tadc_adjust ...
            - mdata.(param_field).radar.wfs(wf).Tadc_adjust;
          delta_offset_Tadc_adjust_mask(wf) = true;
        end
      end
      
      first_idx = find(delta_offset_t_ref_mask,1);
      delta_offset_t_ref(~delta_offset_t_ref_mask) = NaN;
      if ~all(delta_offset_t_ref(delta_offset_t_ref_mask) == delta_offset_t_ref(first_idx))
        delta_offset_t_ref
        error('Different t_ref delta offsets for each waveform, cannot proceed: reprocess data.');
      end
      delta_offset_t_ref = delta_offset_t_ref(first_idx);
      
      first_idx = find(delta_offset_Tsys_mask,1);
      delta_offset_Tsys(~delta_offset_Tsys_mask) = NaN;
      if ~all(delta_offset_Tsys(delta_offset_Tsys_mask) == delta_offset_Tsys(first_idx))
        delta_offset_Tsys
        error('Different Tsys delta offsets for each waveform-adc-pair, cannot proceed: reprocess data.');
      end
      delta_offset_Tsys = delta_offset_Tsys(first_idx);
      
      first_idx = find(delta_offset_Tadc_adjust_mask,1);
      delta_offset_Tadc_adjust(~delta_offset_Tadc_adjust_mask) = NaN;
      if ~all(delta_offset_Tadc_adjust(delta_offset_Tadc_adjust_mask) == delta_offset_Tadc_adjust(first_idx))
        delta_offset_Tadc_adjust
        error('Different Tadc_adjust delta offsets for each waveform, cannot proceed: reprocess data.');
      end
      delta_offset_Tadc_adjust = delta_offset_Tadc_adjust(first_idx);
      
      for img = 1:length(mdata.(param_field).(param_field(7:end)).imgs)
        for wf_adc_pair = 1:size(mdata.(param_field).(param_field(7:end)).imgs{img},1)
          wf = abs(mdata.(param_field).(param_field(7:end)).imgs{img}(wf_adc_pair,1));
          adc = abs(mdata.(param_field).(param_field(7:end)).imgs{img}(wf_adc_pair,2));
          
          mdata.(param_field).radar.wfs(wf).t_ref = param.radar.wfs(wf).t_ref;
          mdata.(param_field).radar.wfs(wf).Tsys(adc) = param.radar.wfs(wf).Tsys(adc);
          mdata.(param_field).radar.wfs(wf).Tadc_adjust = param.radar.wfs(wf).Tadc_adjust;
        end
      end
      
      %% Update Data File
      if delta_offset_t_ref ~= 0 || delta_offset_Tsys ~= 0 || delta_offset_Tadc_adjust ~= 0
        fprintf('  t_ref Offset  %g %s (%s)\n', delta_offset_t_ref, data_fn, datestr(now,'HH:MM:SS'));
        fprintf('  Tsys Offset %g %s (%s)\n', delta_offset_Tsys, data_fn, datestr(now,'HH:MM:SS'));
        fprintf('  Tadc_adjust Offset %g %s (%s)\n', delta_offset_Tadc_adjust, data_fn, datestr(now,'HH:MM:SS'));
        
        % Note that Tadc_adjust is opposite sign to Tsys
        %  - Tsys is subtracted away
        %  - Tadc_adjust is added on
        delta_offset = delta_offset_Tsys - delta_offset_Tadc_adjust - delta_offset_t_ref;
        
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
  
end
