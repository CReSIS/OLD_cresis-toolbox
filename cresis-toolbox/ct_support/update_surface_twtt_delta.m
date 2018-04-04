% script update_surface_twtt_delta
%
% Updates the fast time related variables for changes in Tadc_adjust and Tsys.
% Can be used on any echogram data product (not CSARP_out). However, for
% SAR processed outputs, the data should be reprocessed if the offset is
% very large because large time system time delay errors affect the focussing.
%
% Updates Time, param_{get_heights,csarp}.radar.wfs(wf).Tadc_adjust,
% param_{get_heights,csarp}radar.wfs(wf).Tsys, Surface, and Bottom variables.
%
% Note that only the param_{get_heights,csarp} field that matters gets
% updated. For example, param_records and param_combine will not be used or
% updated because Tadc_adjust and Tsys are not used during these processes
% and so their value does not matter when those processes are applied.
%
% Tsys and Tadc_adjust changes can only be corrected when the same change is applied to all
% waveform-adc pairs used in the data product.
%
% Note that Tadc_adjust is opposite sign to Tsys
%  - Tsys is subtracted away
%  - Tadc_adjust is added on
%
% Examples:
%   See run_update_surface_twtt_delta.m to run
%
% Author: John Paden

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
    
    for data_type_idx = 1:length(data_types)
      data_type = data_types{data_type_idx};
      
      for data_img = imgs
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
        mdata = load(data_fn,'Time','Surface','Bottom','param_get_heights','param_csarp');
        warning on;
        
        fields_to_update = {'Time'};
        
        % Should be only one of these types of fields:
        if isfield(mdata,'param_get_heights')
          delta_offset_Tsys = [];
          delta_offset_Tsys_mask = logical([]);
          delta_offset_Tadc_adjust = [];
          delta_offset_Tadc_adjust_mask = logical([]);
          for img = 1:length(mdata.param_get_heights.get_heights.imgs)
            for wf_adc_pair = 1:size(mdata.param_get_heights.get_heights.imgs{img},1)
              wf = abs(mdata.param_get_heights.get_heights.imgs{img}(wf_adc_pair,1));
              adc = abs(mdata.param_get_heights.get_heights.imgs{img}(wf_adc_pair,2));
              
              delta_offset_Tsys(img,wf_adc_pair) = param.radar.wfs(wf).Tsys(adc) - mdata.param_get_heights.radar.wfs(wf).Tsys(adc);
              delta_offset_Tsys_mask(img,wf_adc_pair) = true;
              
              if ~isfield(param.radar.wfs(wf),'Tadc_adjust') || isempty(param.radar.wfs(wf).Tadc_adjust)
                param.radar.wfs(wf).Tadc_adjust = 0;
              end
              if ~isfield(mdata.param_get_heights.radar.wfs(wf),'Tadc_adjust') ...
                  || isempty(mdata.param_get_heights.radar.wfs(wf).Tadc_adjust)
                mdata.param_get_heights.radar.wfs(wf).Tadc_adjust = 0;
              end
              delta_offset_Tadc_adjust(wf) = param.radar.wfs(wf).Tadc_adjust - mdata.param_get_heights.radar.wfs(wf).Tadc_adjust;
              delta_offset_Tadc_adjust_mask(wf) = true;
            end
          end
          
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
          
          for img = 1:length(mdata.param_get_heights.get_heights.imgs)
            for wf_adc_pair = 1:size(mdata.param_get_heights.get_heights.imgs{img},1)
              wf = abs(mdata.param_get_heights.get_heights.imgs{img}(wf_adc_pair,1));
              adc = abs(mdata.param_get_heights.get_heights.imgs{img}(wf_adc_pair,2));
              
              mdata.param_get_heights.radar.wfs(wf).Tsys(adc) = param.radar.wfs(wf).Tsys(adc);
              mdata.param_get_heights.radar.wfs(wf).Tadc_adjust = param.radar.wfs(wf).Tadc_adjust;
            end
          end
          
          fields_to_update{end+1} = 'param_get_heights';
          
        elseif isfield(mdata,'param_csarp')
          delta_offset_Tsys = [];
          delta_offset_Tsys_mask = logical([]);
          delta_offset_Tadc_adjust = [];
          delta_offset_Tadc_adjust_mask = logical([]);
          for img = 1:length(mdata.param_csarp.csarp.imgs)
            for wf_adc_pair = 1:size(mdata.param_csarp.csarp.imgs{img},1)
              wf = abs(mdata.param_csarp.csarp.imgs{img}(wf_adc_pair,1));
              adc = abs(mdata.param_csarp.csarp.imgs{img}(wf_adc_pair,2));
              
              delta_offset_Tsys(img,wf_adc_pair) = param.radar.wfs(wf).Tsys(adc) - mdata.param_csarp.radar.wfs(wf).Tsys(adc);
              delta_offset_Tsys_mask(img,wf_adc_pair) = true;
              
              if ~isfield(param.radar.wfs(wf),'Tadc_adjust') || isempty(param.radar.wfs(wf).Tadc_adjust)
                param.radar.wfs(wf).Tadc_adjust = 0;
              end
              if ~isfield(mdata.param_csarp.radar.wfs(wf),'Tadc_adjust') ...
                  || isempty(mdata.param_csarp.radar.wfs(wf).Tadc_adjust)
                mdata.param_csarp.radar.wfs(wf).Tadc_adjust = 0;
              end
              delta_offset_Tadc_adjust(wf) = param.radar.wfs(wf).Tadc_adjust - mdata.param_csarp.radar.wfs(wf).Tadc_adjust;
              delta_offset_Tadc_adjust_mask(wf) = true;
            end
          end
          
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
          
          for img = 1:length(mdata.param_csarp.csarp.imgs)
            for wf_adc_pair = 1:size(mdata.param_csarp.csarp.imgs{img},1)
              wf = abs(mdata.param_csarp.csarp.imgs{img}(wf_adc_pair,1));
              adc = abs(mdata.param_csarp.csarp.imgs{img}(wf_adc_pair,2));
              
              mdata.param_csarp.radar.wfs(wf).Tsys(adc) = param.radar.wfs(wf).Tsys(adc);
              mdata.param_csarp.radar.wfs(wf).Tadc_adjust = param.radar.wfs(wf).Tadc_adjust;
            end
          end
          
          fields_to_update{end+1} = 'param_csarp';
        end
        
        if delta_offset_Tsys ~= 0 || delta_offset_Tadc_adjust ~= 0
          fprintf('  Tsys Offset %g %s (%s)\n', delta_offset_Tsys, data_fn, datestr(now,'HH:MM:SS'));
          fprintf('  Tadc_adjust Offset %g %s (%s)\n', delta_offset_Tadc_adjust, data_fn, datestr(now,'HH:MM:SS'));
          
          % Note that Tadc_adjust is opposite sign to Tsys
          %  - Tsys is subtracted away
          %  - Tadc_adjust is added on
          delta_offset = delta_offset_Tsys - delta_offset_Tadc_adjust;
          
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
  
end

return;
