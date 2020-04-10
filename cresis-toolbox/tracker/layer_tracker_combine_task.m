function success = layer_tracker_combine_task(param)
%function to combine layer data information received from the tracker. To
%store in layerData files. Use opsCopyLayers
%global gRadar;
if (param.layer_tracker.save_ops_copy_layers)
  layer_params = param.layer_tracker.cmds.layer_params;
  for layer_idx = 1:length(layer_params)
    TWTT{layer_idx} = [];%lay.layerData{layer_idx}.quality = [];
    %     lay.layerData{layer_idx}.value{2}.data =[];
  end
  if param.layer_tracker.track.save_layerData
    
    for iter = 1:16
      for layer_idx = 1:length(layer_params)
        TWTT{layer_idx,iter} = [];
      end
    end
  end
  
  GPS_time = [];
  numFrms = [];
  
  
  for idx = 1:length(param.tmp_out_fn)
    tmp_out_fn = load(param.tmp_out_fn{idx});
    if param.layer_tracker.track.save_layerData
      for iter = 1:16
        for layer_idx = 1:length(layer_params)
          TWTT{layer_idx,iter} = cat(2,TWTT{layer_idx,iter},tmp_out_fn.twtt{layer_idx,iter});;
        end
      end
    else
      %for frm = param.totalfrms{idx}(1):param.totalfrms{idx}(2)
      for layer_idx = 1:length(layer_params)
        %         layer_param = layer_params(layer_idx);
        %         layer_fn = fullfile(ct_filename_out(param,layer_param{1}.source,''),sprintf('Data_%s_%03d.mat', param.day_seg, frm));
        %         lay = load(layer_fn);
        %         temp = find(tmp_out_fn.indexes{layer_idx}==frm);
        %         tmp_out_fn.TWTT = interp1(1:length(tmp_out_fn.gps_time),tmp_out_fn.gps_time,tmp_out_fn.twtt{layer_idx});
        TWTT{layer_idx} = cat(2,TWTT{layer_idx},tmp_out_fn.twtt{layer_idx});
        %         tmp_out_fn.TWTT = interp_finite(tmp_out_fn.TWTT,NaN);
        %         tmp_out_fn.GPS_time = tmp_out_fn.gps_time;
        %         lay.layerData{layer_idx}.quality = interp1(tmp_out_fn.gps_time(temp),ones(1,length(temp)),lay.GPS_time,'nearest');
        %         lay.layerData{layer_idx}.value{2}.data = interp1(tmp_out_fn.gps_time(temp),tmp_out_fn.twtt{layer_idx}(length(temp)),lay.GPS_time);
        %tmp_info = Ldata.tracker_data{layer_idx}.(sprintf('%s_%s_%03d',param.layer_tracker.track.method,param.layer_tracker.name,layer_idx));
      end
    end
    %layerData = lay.layerData;
    GPS_time = cat(2,GPS_time,tmp_out_fn.gps_time);
    %save(layer_fn,'-append','layerData');
    %end
    numFrms = cat(2,numFrms,(param.totalfrms{idx}(1):param.totalfrms{idx}(2)));
    param.cmd.frms = numFrms;
  end
  
  
  % lay.GPS_time+1:lay.GPS_time*idx(which ranges from 1, 2, 3... end)
  % for idx = 1:length(numfolders)
  %   layer_fn = fullfile(ct_filename_out(param,layer_param{1}.source,''),sprintf('Data_%s_%03d.mat', param.day_seg, frm));
  %   lay = load(layer_fn);
  %   fname = fullfile(fpath,sprintf('layer_tracker_%03d_%03d',idx,idx),sprintf('layer_%s_%s',param.layer_tracker.track.method,param.layer_tracker.name));
  %   Ldata  = load(fname);
  %
  % end
  % lay.layerData{layer_idx}.value{2}.data = interp1(tmp_info.GPS_time,tmp_info.layerData{layer_idx}.value{2}.data,lay.GPS_time);
  %
  %keyboard;
  copy_param = [];
  copy_param.layer_source.existence_check = false;
  copy_param.layer_dest.existence_check = false;
  
  % Set the source
  copy_param.layer_source.source = 'custom';
  copy_param.layer_dest.source = param.layer_tracker.track.layer_dest_source;
  
  if strcmp(copy_param.layer_dest.source, 'layerdata')
    copy_param.layer_dest.layerdata_source = param.layer_tracker.track.layer_dest_layerdata_source; %layerData_test
    copy_param.layer_dest.echogram_source = param.layer_tracker.track.layer_dest_echogram_source; %standard
  end
  
  copy_param.copy_method = 'overwrite';
  
  %if strcmpi(copy_param.layer_dest.name,'surface')
  
  %else
  
  %end
  
  %keyboard;
  %for surface
  if strcmp(param.layer_tracker.track.method,'viterbi')
    copy_param.layer_source.gps_time = GPS_time;
    copy_param.layer_source.twtt = TWTT;
    
    copy_param.gaps_fill.method = 'preserve_gaps';
    copy_param.gaps_fill.method_args = [40 20];
    
    copy_param.layer_dest.name = 'bottom';%param.layer_tracker.track.(sprintf('%s',param.layer_tracker.track.method)).lyrbot;
    fprintf('opsCopyLayers %s (%s)\n', param.day_seg, datestr(now));
    opsCopyLayers(param,copy_param);
  else
    if param.layer_tracker.track.save_layerData
      for iter = 1:16
        copy_param.layer_source.gps_time = GPS_time;
        copy_param.layer_source.twtt = TWTT{1,iter};
        copy_param.gaps_fill.method = 'interp_finite';
        copy_param.layer_dest.name = sprintf('lsm_%s_surf_%03d',param.layer_tracker.name,iter);
        fprintf('opsCopyLayers %s (%s)\n', param.day_seg, datestr(now));
        opsCopyLayers(param,copy_param);
      end
      for iter = 1:16
        copy_param.layer_source.gps_time = GPS_time;
        copy_param.layer_source.twtt = TWTT{2,iter};
        
        copy_param.gaps_fill.method = 'preserve_gaps';
        copy_param.gaps_fill.method_args = [40 20];
        
        copy_param.layer_dest.name = sprintf('lsm_%s_bot_%03d',param.layer_tracker.name,iter);%param.layer_tracker.track.(sprintf('%s',param.layer_tracker.track.method)).lyrbot;
        fprintf('opsCopyLayers %s (%s)\n', param.day_seg, datestr(now));
        opsCopyLayers(param,copy_param);
      end
    else
      
      copy_param.layer_source.gps_time = GPS_time;
      copy_param.layer_source.twtt = TWTT{1};
      
      copy_param.gaps_fill.method = 'interp_finite';
      
      copy_param.layer_dest.name = 'surface';%param.layer_tracker.track.(sprintf('%s',param.layer_tracker.track.method)).lyrtop;
      fprintf('opsCopyLayers %s (%s)\n', param.day_seg, datestr(now));
      opsCopyLayers(param,copy_param);
      
      if any(strcmp(param.layer_tracker.track.method,'lsm')) || any(strcmp(param.layer_tracker.track.method,'stereo')) || any(strcmp(param.layer_tracker.track.method,'mcmc'))
        % for bottom

        copy_param.layer_source.gps_time = GPS_time;
        copy_param.layer_source.twtt = TWTT{2};
        
        copy_param.gaps_fill.method = 'preserve_gaps';
        copy_param.gaps_fill.method_args = [40 20];
        
        copy_param.layer_dest.name = 'bottom';%param.layer_tracker.track.(sprintf('%s',param.layer_tracker.track.method)).lyrbot;
        fprintf('opsCopyLayers %s (%s)\n', param.day_seg, datestr(now));
        opsCopyLayers(param,copy_param);
      end
    end
  end
  
  fprintf('Complete (%s)\n', datestr(now));
  warning('on');
end
success=true;
end
% surface and bottom can be named as the layer name. Change param.cmd.frms
% to get the number of frames that you are working with. Look into saving
% the layerData info and viewing the results in imb.picker or by plotting a
% graph. Group,description,source,echogram_source, set the propoerties by
% User. Look into other ways to save layername of the destination.
