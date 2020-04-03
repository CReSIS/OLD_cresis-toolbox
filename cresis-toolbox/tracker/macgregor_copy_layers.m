% script macgregor_copy_layers.m
%
% Example script for running opsCopyLayers.m to copy internal layers based
% on MacGregor et al layer tracking tool.
%
% Authors: Jilu Li

% =====================================================================
%% User Settings
% =====================================================================
% clear gris_strat
global gRadar;
param = [];
gris_fns = {'/cresis/snfs1/dataproducts/metadata/macgregor_gris_layers/Greenland_radiostratigraphy.mat'};
if 1
  load(gris_fns{1});
end
for yr_idx = 19:19 %1:gris_strat.num_campaign
  year = gris_strat.campaign(yr_idx).name(1:4);
  aircraft = gris_strat.campaign(yr_idx).name(6:7);
  param_fn = sprintf('rds_param_%4s_Greenland_%2s.xls',year,aircraft);
  params = read_param_xls(ct_filename_param(param_fn));
  params = ct_set_params(params,'cmd.generic',0);
  params = ct_set_params(params,'cmd.frms',[]);
  
  copy_param = [];
  copy_param.layer_source.name = 'macgregor';
  copy_param.layer_source.source = 'custom';
  copy_param.layer_dest.existence_check = false;
  copy_param.layer_dest.source = 'ops';
  copy_param.layer_dest.group = 'macgregor';
  copy_param.copy_method = 'overwrite';
  copy_param.gaps_fill.method = 'preserve_gaps';
  copy_param.gaps_fill.method_args = [40 20];
  for seg_idx = 5:5 %1:gris_strat.campaign(yr_idx).num_segment
    seg = gris_strat.campaign(yr_idx).segment(seg_idx).name;
    params = ct_set_params(params,'cmd.generic',1,'day_seg',seg);
    for param_idx = 1:length(params)
      if params(param_idx).cmd.generic
        param = params(param_idx);
      end
      param = merge_structs(param,gRadar);
    end
    for div_idx = 1:gris_strat.campaign(yr_idx).segment(seg_idx).num_division
      for lyr_idx = 257:gris_strat.campaign(yr_idx).segment(seg_idx).division(div_idx).num_layer
        % skip empty layers without any twtt values
        if ~all(isnan(gris_strat.campaign(yr_idx).segment(seg_idx).division(div_idx).layer(lyr_idx).traveltime))
          % remove points without twtt values
          good_idxs = find(~isnan(gris_strat.campaign(yr_idx).segment(seg_idx).division(div_idx).layer(lyr_idx).traveltime));
          copy_param.layer_source.gps_time = gris_strat.campaign(yr_idx).segment(seg_idx).division(div_idx).time_gps(good_idxs);
          copy_param.layer_source.twtt = gris_strat.campaign(yr_idx).segment(seg_idx).division(div_idx).layer(lyr_idx).traveltime(good_idxs);
          copy_param.layer_source.quality = ones(size(copy_param.layer_source.twtt));
          layer_age = gris_strat.campaign(yr_idx).segment(seg_idx).division(div_idx).layer(lyr_idx).age;
          if ~isnan(layer_age)
            copy_param.layer_dest.name = sprintf('layer_1950-%06d',round(layer_age));
          else
            copy_param.layer_dest.name = sprintf('L%03d-%s-D%02d',lyr_idx,param.day_seg,div_idx);
          end
          copy_param.layer_dest.description = sprintf('The %3dth MacGregor picked layer in %s_%02d.',lyr_idx,seg,div_idx);
          fprintf('opsCopyLayers %s (%s)\n', param.day_seg, datestr(now));
          opsCopyLayers(param,copy_param);
          fprintf('  layer %3d complete (%s)\n', lyr_idx,datestr(now));
        end
      end
      if 0 % debug plot
        param.season_name = params(1).season_name;
        param.radar_name = params(1).radar_name;
        param.day_seg = seg;
        frm = 2;
        fn = fullfile(ct_filename_out(param,'CSARP_post/standard',''),...
          sprintf('Data_%s_%03d.mat', seg, frm));
        tmp = load(fn);
        good_idxs = find(copy_param.layer_source.gps_time>=tmp.GPS_time(1) & copy_param.layer_source.gps_time<=tmp.GPS_time(end));
        figure(1);clf;imagesc(tmp.GPS_time,tmp.Time*1e6,lp(tmp.Data));
        colormap(1-gray)
        for lyr_idx = 1:gris_strat.campaign(yr_idx).segment(seg_idx).division(div_idx).num_layer
          figure(1);hold on;plot(copy_param.layer_source.gps_time(good_idxs),...
            gris_strat.campaign(yr_idx).segment(seg_idx).division(div_idx).layer(lyr_idx).traveltime(good_idxs)*1e6);
        end
      end
    end
  end
end
