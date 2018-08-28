% script update_mult_factor
%
% Updates the mult_factor field (adc_gains) in specified files.
%
% Examples:
%   See run_update_mult_factor.m to run
%
% Author: John Paden

%% Automated Section
% ----------------------------------------------------------------------
for param_idx = 1:length(params)
  param = params(param_idx);
  
  if ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  fprintf('Updating mult_factor %s (%s)\n', param.day_seg, datestr(now,'HH:MM:SS'));
  
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
  
  for data_type_idx = 1:length(data_types)
    data_type = data_types(data_type_idx);
    
    if strcmpi(data_type.type,'coh_noise')
      
      
      for img = data_type.imgs
        
        if isempty(data_type.wf_adcs)
          wf_adcs = 1:size(param.analysis.imgs{img},1);
        else
          wf_adcs = data_type.wf_adcs;
        end
        for wf_adc = wf_adcs
          wf = param.analysis.imgs{img}(wf_adcs,1);
          adc = param.analysis.imgs{img}(wf_adcs,2);
          
          % Update analysis_combine_task output
          fn_dir = fileparts(ct_filename_out(param,data_type.dir));
          fn = fullfile(fn_dir,sprintf('coh_noise_%s_wf_%d_adc_%d.mat', param.day_seg, wf, adc));
          if ~exist(fn,'file')
            warning('Missing %s\n', fn);
          else
            fprintf('Loading %s (%s)\n', fn, datestr(now));
            noise = load(fn);
            d_adc_gain = param.radar.wfs(wf).adc_gains(adc) / noise.param_analysis.radar.wfs(wf).adc_gains(adc);
            if abs(d_adc_gain - 1) > 1e-6
              fprintf('  Gain change of %.3f\n', d_adc_gain);
              
              % Update fields
              for idx=1:length(noise.coh_ave)
                noise.coh_ave{idx} = noise.coh_ave{idx} / d_adc_gain;
              end
              noise.doppler = noise.doppler / d_adc_gain.^2;
              noise.param_analysis.radar.wfs(wf).adc_gains(adc) = param.radar.wfs(wf).adc_gains(adc);
              
              % Save result
              fprintf('  Saving change\n');
              ct_file_lock_check(fn);
              save(fn,'-v7.3','-struct','noise');
            end
          end
          
          % Update collate_coh_noise output
          fn_dir = fileparts(ct_filename_out(param,data_type.dir));
          fn = fullfile(fn_dir,sprintf('coh_noise_simp_%s_wf_%d_adc_%d.nc', param.day_seg, wf, adc));
          if ~exist(fn,'file')
            warning('Missing %s\n', fn);
          else
            fprintf('Loading %s (%s)\n', fn, datestr(now));
            noise = netcdf_to_mat(fn);
            
            d_adc_gain = param.radar.wfs(wf).adc_gains(adc) / noise.param_analysis.radar.wfs(wf).adc_gains(adc);
            if abs(d_adc_gain - 1) > 1e-6
              fprintf('  Gain change of %.3f\n', d_adc_gain);
              
              % Update fields
              noise.dftI = noise.dftI / d_adc_gain;
              noise.dftQ = noise.dftQ / d_adc_gain;
              noise.param_analysis.radar.wfs(wf).adc_gains(adc) = param.radar.wfs(wf).adc_gains(adc);
              
              % Save result
              fprintf('  Saving change\n');
              ct_file_lock_check(fn);
              netcdf_from_mat(fn,noise);
            end
          end
        end
      end
      
      
    elseif strcmpi(data_type.type,'echogram')
      
      for frm_idx = 1:length(param.cmd.frms)
        frm = param.cmd.frms(frm_idx);
        
        for data_img = data_type.imgs
          if data_img == 0
            fn = fullfile(ct_filename_out(param,data_type.dir,''), ...
              sprintf('Data_%s_%03d.mat', param.day_seg, frm));
          else
            fn = fullfile(ct_filename_out(param,data_type.dir,''), ...
              sprintf('Data_img_%02d_%s_%03d.mat', data_img, param.day_seg, frm));
          end
          
          if ~exist(fn,'file')
            warning('Missing %s\n', fn);
            continue;
          end
          
          fprintf('Loading %s (%s)\n', fn, datestr(now));
          mdata = load(fn);
          
          % Should be only one of these types of fields:
          if isfield(mdata,'param_qlook')
            
            d_adc_gain = param.radar.wfs(wf).adc_gains(adc) / mdata.param_qlook.radar.wfs(wf).adc_gains(adc);
            if abs(d_adc_gain - 1) > 1e-6
              fprintf('  Gain change of %.3f\n', d_adc_gain);
              
              % Update fields
              if mdata.param_qlook.qlook.inc_dec > 0
                % Incoherent data
                mdata.Data = mdata.Data / d_adc_gain.^2;
              else
                % Coherent data
                mdata.Data = mdata.Data / d_adc_gain;
              end
              mdata.param_qlook.radar.wfs(wf).adc_gains(adc) = param.radar.wfs(wf).adc_gains(adc);
              
              % Save result
              fprintf('  Saving change\n');
              ct_file_lock_check(fn);
              save(fn,'-v7.3','-struct','mdata');
            end
          end
          
        end
        
      end
    end
    
  end
  
end

return;
