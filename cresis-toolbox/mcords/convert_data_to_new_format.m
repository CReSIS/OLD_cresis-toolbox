function convert_data_to_new_format
% convert_data_to_new_format
%
% Function for converting from old data files (~Fall 2010) to V1.0 of toolbox
%
% Convert from:
% /cresis/scratch2/mdce/old/mcords/2009_Antarctica_DC8/PTSK
%   Data_20091107_seg1_06_01.mat
%          Data: [625x4840 single]
%         Depth: [1x625 double]
%     Elevation: [1x4840 double]
%      GPS_time: [1x4840 double]
%      Latitude: [1x4840 double]
%     Longitude: [1x4840 double]
%       Surface: [1x4840 double]
%          Time: [1x625 double]
% /cresis/scratch2/mdce/old/mcords/2009_Antarctica_DC8/CSARP_layerData/PTSK
%    layerData []
%
% To:
% /cresis/scratch2/mdce/mcords/2009_Antarctica_DC8/CSARP_standard/20091016_01/
%   Data_20091016_01_001.mat
%   Data_img_01_20091016_01_001.mat
%   Data_img_02_20091016_01_001.mat
%            Bottom: [1x3706 double]
%              Data: [436x3706 double]
%             Depth: [1x436 double]
%         Elevation: [1x3706 double]
%          GPS_time: [1x3706 double]
%          Latitude: [1x3706 double]
%         Longitude: [1x3706 double]
%           Surface: [1x3706 double]
%              Time: [1x436 double]
%       array_param: [1x1 struct]
%       param_csarp: [1x1 struct]
%       param_radar: [1x1 struct]
%     param_records: [1x1 struct]
% /cresis/scratch2/mdce/mcords/2009_Antarctica_DC8/CSARP_layerData/20091016_01/
%    GPS_time
%    Latitude
%    Longitude
%    Elevation
%    layerData
%
% Author: John Paden

in_path = {}; out_path = {}; out_type = {};

if 0
  % Getz, Abbott, Arc86, Peninsula, PTSK
  idx = 1;
  in_path{idx} = '/cresis/scratch2/mdce/old/mcords/2009_Antarctica_DC8/Peninsula';
  out_path{idx} = '/cresis/scratch2/mdce/mcords/2009_Antarctica_DC8/CSARP_standard/';
  out_type{idx} = '';
  idx = idx + 1;
  in_path{idx} = '/cresis/scratch2/mdce/old/mcords/2009_Antarctica_DC8/Peninsula_wf_01';
  out_path{idx} = '/cresis/scratch2/mdce/mcords/2009_Antarctica_DC8/CSARP_standard/';
  out_type{idx} = 'img_01_';
  idx = idx + 1;
  in_path{idx} = '/cresis/scratch2/mdce/old/mcords/2009_Antarctica_DC8/Peninsula_MVDR';
  out_path{idx} = '/cresis/scratch2/mdce/mcords/2009_Antarctica_DC8/CSARP_mvdr/';
  out_type{idx} = '';
  idx = idx + 1;
  in_path{idx} = '/cresis/scratch2/mdce/old/mcords/2009_Antarctica_DC8/Peninsula_QLook';
  out_path{idx} = '/cresis/scratch2/mdce/mcords/2009_Antarctica_DC8/CSARP_qlook/';
  out_type{idx} = 'img_01_';

  in_layer_path = '/cresis/scratch2/mdce/old/mcords/2009_Antarctica_DC8/CSARP_layerData/Peninsula';
  out_layer_path = '/cresis/scratch2/mdce/mcords/2009_Antarctica_DC8/CSARP_layerData/';
elseif 0
  % Northeast, Northwest, Rink, West, Southeast
  idx = 1;
  in_path{idx} = '/cresis/scratch2/mdce/old/mcords/2010_Greenland_DC8/Southeast';
  out_path{idx} = '/cresis/scratch2/mdce/mcords/2010_Greenland_DC8/CSARP_standard/';
  out_type{idx} = '';
  idx = idx + 1;
  in_path{idx} = '/cresis/scratch2/mdce/old/mcords/2010_Greenland_DC8/Southeast_wf_01';
  out_path{idx} = '/cresis/scratch2/mdce/mcords/2010_Greenland_DC8/CSARP_standard/';
  out_type{idx} = 'img_01_';
  idx = idx + 1;
  in_path{idx} = '/cresis/scratch2/mdce/old/mcords/2010_Greenland_DC8/Southeast_MVDR';
  out_path{idx} = '/cresis/scratch2/mdce/mcords/2010_Greenland_DC8/CSARP_mvdr/';
  out_type{idx} = '';
  idx = idx + 1;
  in_path{idx} = '/cresis/scratch2/mdce/old/mcords/2010_Greenland_DC8/Southeast_QLook';
  out_path{idx} = '/cresis/scratch2/mdce/mcords/2010_Greenland_DC8/CSARP_qlook/';
  out_type{idx} = 'img_01_';

  in_layer_path = '/cresis/scratch2/mdce/old/mcords/2010_Greenland_DC8/CSARP_layerData/Southeast';
  out_layer_path = '/cresis/scratch2/mdce/mcords/2010_Greenland_DC8/CSARP_layerData/';
elseif 1
  % Separated by day_seg
  idx = 1;
  in_path{idx} = '/cresis/scratch2/mdce/old/mcords/2010_Greenland_P3/CSARP_standard/';
  out_path{idx} = '/cresis/scratch2/mdce/mcords/2010_Greenland_P3/CSARP_standard/';
  out_type{idx} = '';
  idx = idx + 1;
  in_path{idx} = '/cresis/scratch2/mdce/old/mcords/2010_Greenland_P3/CSARP_mvdr';
  out_path{idx} = '/cresis/scratch2/mdce/mcords/2010_Greenland_P3/CSARP_mvdr/';
  out_type{idx} = '';

  in_layer_path = '/cresis/scratch2/mdce/old/mcords/2010_Greenland_P3/CSARP_layerData/';
  out_layer_path = '/cresis/scratch2/mdce/mcords/2010_Greenland_P3/CSARP_layerData/';
end

in_fns = get_filenames(in_path{1},'Data_2','seg','.mat',struct('recursive',1));

fprintf('Converting files\n');
for file_idx = 1:length(in_fns)
  fn = in_fns{file_idx};
  
  [fn_dir fn_name] = fileparts(fn);
  fprintf('  %s\n', fn_name);
  
  date_str = fn_name(6:13);
  seg_num = sscanf(fn_name(18:19),'%f');
  day_seg = sprintf('%s_%02d', date_str, seg_num);
  day_seg_old = sprintf('%s_seg%d', date_str, seg_num);
  
  if seg_num < 10
    frm = str2double(fn_name(20:21));
  else
    frm = str2double(fn_name(21:22));
  end

  for in_idx = 1:length(in_path)
    out_fn_name = sprintf('Data_%s%s_%03d.mat', out_type{in_idx}, day_seg, frm);
    in_fn = fullfile(in_path{in_idx},day_seg_old,[fn_name '.mat']);
    out_fn = fullfile(out_path{in_idx},day_seg,out_fn_name);
    out_fn_dir = fileparts(out_fn);
    if ~exist(out_fn_dir,'dir')
      mkdir(out_fn_dir);
    end
    cmd = sprintf('cp %s %s', in_fn, out_fn);
    system(cmd);
  end

  out_fn_name = sprintf('Data_%s_%03d.mat', day_seg, frm);
  in_fn_name = sprintf('Data_%s_%02d.mat', day_seg_old, frm);
  in_fn = fullfile(in_layer_path,[fn_name '.mat']);
  out_fn = fullfile(out_layer_path,day_seg,out_fn_name);
  out_fn_dir = fileparts(out_fn);
  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  load(fn);
  load(in_fn);
  
  save(out_fn,'layerData','GPS_time','Latitude','Longitude','Elevation');
  
end

return;
