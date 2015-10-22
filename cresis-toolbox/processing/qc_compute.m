function qc_metric = qc_compute(img_fname)
% qc_metric = qc_compute(img_fname)
%
% Inputs:
% img_fname
% E.g.
% '/cresis/scratch2/mdce/mcords2/2011_Greenland_P3/CSARP_standard/20110507_01/Data_20110507_01_028.mat'
%
% landmarks is a file containing the parameter fields for
%   signal region (region for extracting the peak signal power)
%   clutter region (region for extracting the average clutter power)
%   noise region (region for extracting the average noise power)
%
% Outputs:
% qc_metric is a structure containing
%   .peak_snr which is the peak signal to noise ratio
%   .peak_scr which is the peak signal to clutter ratio
%   .peak_cnr which is the average clutter to noise ratio
%   .clutter_power which is the average clutter power
%   .noise_power which is the average noise power
%   
% Author: Isaac Tan
%
% Example at end of file

% ====================================================================
% Process and check arguments
% =====================================================================

if ~isempty(figure(201))
  close figure 201
end

clear landmarks
clc

% Load landmarks.mat file
% qc_landmark_fname = '/users/petan/scripts/matlab/landmarks.mat';
qc_landmark_fname = '/cresis/scratch2/mdce/landmarks.mat';

if exist(qc_landmark_fname,'file')
  load(qc_landmark_fname);
else
  error('The corresponding landmark data is missing');
end

if isempty(img_fname)
  
  array_season_name = {};
  array_img_id = {};
  array_landmark_type = {};
  
  landmark_count = 1;
  size_landmark = size(landmarks);

  season_name = landmarks(landmark_count).season_name;
  img_id = landmarks(landmark_count).frm_id;
  landmark_type = landmarks(landmark_count).type;
  
  array_season_name{end+1} = season_name;
  array_img_id{end+1} = img_id;
  array_landmark_type{end+1} = landmark_type;
  array_length = 1;
  
  qc_landmarks = {};
  qc_landmarks{end+1} = 'bottom';
  qc_landmarks{end+1} = 'clutter';
  qc_landmarks{end+1} = 'noise';
  

  for landmark_count = 2:size_landmark(2)

    season_name = landmarks(landmark_count).season_name;
    img_id = landmarks(landmark_count).frm_id;
    landmark_type = landmarks(landmark_count).type;

    is_qc_landmark = ismember(landmark_type,qc_landmarks);
    
    if is_qc_landmark
    
      member_season = find(~ismember(season_name,array_season_name{end}));
  %     member_img_id = find(~ismember(img_id,array_img_id{end}));
      member_img_id = strcmp(img_id,array_img_id{end});
      member_landmark_type = find(~ismember(landmark_type,array_landmark_type{end}));

      if ~isempty(member_season) || ~member_img_id
        array_season_name{end+1} = season_name;
        array_img_id{end+1} = img_id;
        array_landmark_type{end+1} = landmark_type;
        array_length = array_length + 1;
      elseif ~isempty(member_landmark_type)
        array_landmark_type{array_length} = [array_landmark_type{array_length},' ',landmark_type];
      end
      
    end

  end
    
  fprintf('Number of QC frames is %d\n\n',array_length);
  
  for m = 1:array_length
    fprintf('Season : %s, Frm_id : %s, Type: %s \n', array_season_name{m}, array_img_id{m}, array_landmark_type{m});
  end
  
  qc_metric.peak_snr = '';
  qc_metric.peak_scr = '';
  qc_metric.peak_cnr = '';
  
  return
end

test_frame = load (img_fname);
img_fname_short = img_fname(end-23:end-4);

img_fname_long = img_fname(end-36:end-4);

find_CSARP = strfind(img_fname,'CSARP');

layer_dir = img_fname(1:find_CSARP+4);
layer_dir = strcat(layer_dir,'_layerData');
layer_fname = strcat(layer_dir,img_fname_long,'.mat');

if exist(layer_fname,'file')
  load(layer_fname);
else
  error('The corresponding layer data is missing');
end

% This section extracts the peak signal power from the selected signal region
find_season = strfind(img_fname(1:find_CSARP-2),'/');
input_season_name = img_fname(find_season(end)+1:find_CSARP-2);

size_landmark = size(landmarks);

landmark_count = 1;
landmark_sel = [];

while landmark_count <= size_landmark(2)
  
  img_id = landmarks(landmark_count).frm_id;
  season_name = landmarks(landmark_count).season_name;
  
  if strcmp(img_id,img_fname_short(6:end)) && strcmp(season_name,input_season_name)
    landmark_sel = [landmark_sel landmark_count];
  end
  landmark_count = landmark_count+1;
end

if isempty(landmark_sel)
  error('No qc landmark of current filename is present in the landmarks.mat');
end

length_num_landmark = length(landmark_sel);

signal_landmark = 0;
clutter_landmark = [];
noise_landmark = 0;
for i=1:length_num_landmark
  if strcmp(landmarks(landmark_sel(i)).type,'bottom')
    signal_landmark = landmark_sel(i);
  end
  
  if strcmp(landmarks(landmark_sel(i)).type,'clutter')
    clutter_landmark = [clutter_landmark landmark_sel(i)];
  end  

  if strcmp(landmarks(landmark_sel(i)).type,'noise')
    noise_landmark = landmark_sel(i);
  end
end

% This section computes the SNR
if signal_landmark == 0
  fprintf('No signal landmark for current filename and thus no SNR computed\n');
elseif noise_landmark == 0
  fprintf('No noise landmark for current filename and thus no SNR computed\n');
else
  % This section extracts the average signal power from the selected
  % signal region

%   sig_x1 = round(landmarks(signal_landmark).rline_start);
%   sig_x2 = round(landmarks(signal_landmark).rline_stop);

  sig_x1 = find(abs(test_frame.GPS_time - landmarks(signal_landmark).gpstime_start) < 1.0,1);
  sig_x2 = find(abs(test_frame.GPS_time - landmarks(signal_landmark).gpstime_stop) < 1.0,1);
  
  sig_y1 = find(abs(test_frame.Time - landmarks(signal_landmark).rbin_start*1.0e-6) < 2.5e-8,1);
  sig_y2 = find(abs(test_frame.Time - landmarks(signal_landmark).rbin_stop*1.0e-6) < 2.5e-8,1);

  bottom_layer_rbins = layerData{2}.value{2};
  good_idxs = find(isfinite(bottom_layer_rbins.data));

  interpolated_bottom_layer = interp1(1:length(bottom_layer_rbins.data(good_idxs)),...
        bottom_layer_rbins.data(good_idxs),1:length(test_frame.Data),'linear','extrap');

  num_sig_rlines = sig_x2 - sig_x1 + 1;

  sig_power = 0;
  empty_count = 0;
  for i = 1:num_sig_rlines
%     tmp1 = find(abs(test_frame.Time - bottom_layer_rbins.data(sig_x1+i)) < 2.5e-8,1);
    tmp1 = find(abs(test_frame.Time - interpolated_bottom_layer(sig_x1+i)) < 2.5e-8,1);
    tmp2 = max(test_frame.Data([tmp1-5:tmp1+5],sig_x1+i));
    if ~isempty(tmp2)
      sig_power = sig_power + tmp2;
    else
      empty_count = empty_count+1;
    end
  end

  sig_power = sig_power/(num_sig_rlines-empty_count);
  
  % This section extracts the average noise power from the selected noise region
%   noise_x1 = round(landmarks(noise_landmark).rline_start);
%   noise_x2 = round(landmarks(noise_landmark).rline_stop);
  
  noise_x1 = find(abs(test_frame.GPS_time - landmarks(noise_landmark).gpstime_start) < 1.0,1);
  noise_x2 = find(abs(test_frame.GPS_time - landmarks(noise_landmark).gpstime_stop) < 1.0,1);
  
  noise_y1 = find(abs(test_frame.Time - landmarks(noise_landmark).rbin_start*1.0e-6) < 2.5e-8,1);
  noise_y2 = find(abs(test_frame.Time - landmarks(noise_landmark).rbin_stop*1.0e-6) < 2.5e-8,1);

  noise_power = mean(mean(test_frame.Data([noise_y1:noise_y2],[noise_x1:noise_x2])));
  
  qc_metric.peak_snr = 10*log10(sig_power/noise_power);

end

% This section computes the SCR
if signal_landmark == 0
  fprintf('No signal landmark for current filename and thus no SCR computed\n');
elseif isempty(clutter_landmark)
  fprintf('No clutter landmark for current filename and thus no SCR computed\n');
else
  
 if noise_landmark == 0 %This means that the signal power has not been computed yet  
  
    % This section extracts the average signal power from the selected
    % signal region

%     sig_x1 = round(landmarks(signal_landmark).rline_start);
%     sig_x2 = round(landmarks(signal_landmark).rline_stop);

    sig_x1 = find(abs(test_frame.GPS_time - landmarks(signal_landmark).gpstime_start) < 1.0,1);
    sig_x2 = find(abs(test_frame.GPS_time - landmarks(signal_landmark).gpstime_stop) < 1.0,1);

    sig_y1 = find(abs(test_frame.Time - landmarks(signal_landmark).rbin_start*1.0e-6) < 2.5e-8,1);
    sig_y2 = find(abs(test_frame.Time - landmarks(signal_landmark).rbin_stop*1.0e-6) < 2.5e-8,1);

  bottom_layer_rbins = layerData{2}.value{2};
  good_idxs = find(isfinite(bottom_layer_rbins.data));

  interpolated_bottom_layer = interp1(1:length(bottom_layer_rbins.data(good_idxs)),...
        bottom_layer_rbins.data(good_idxs),1:length(test_frame.Data),'linear','extrap');

  num_sig_rlines = sig_x2 - sig_x1 + 1;

  sig_power = 0;
  empty_count = 0;
  for i = 1:num_sig_rlines
%     tmp1 = find(abs(test_frame.Time - bottom_layer_rbins.data(sig_x1+i)) < 2.5e-8,1);
    tmp1 = find(abs(test_frame.Time - interpolated_bottom_layer(sig_x1+i)) < 2.5e-8,1);
    tmp2 = max(test_frame.Data([tmp1-5:tmp1+5],sig_x1+i));
    if ~isempty(tmp2)
      sig_power = sig_power + tmp2;
    else
      empty_count = empty_count+1;
    end
  end

  sig_power = sig_power/(num_sig_rlines-empty_count);
 
 end
  
  % This section extracts the average clutter power from the various selected
  % clutter regions
  num_clutter = length(clutter_landmark);
  clu_power = 0;

  for i = 1:num_clutter
%     clu_x1 = round(landmarks(clutter_landmark(i)).rline_start);
%     clu_x2 = round(landmarks(clutter_landmark(i)).rline_stop);
    clu_x1 = find(abs(test_frame.GPS_time - landmarks(clutter_landmark(i)).gpstime_start) < 1.0,1);
    clu_x2 = find(abs(test_frame.GPS_time - landmarks(clutter_landmark(i)).gpstime_stop) < 1.0,1);

    clu_y1 = find(abs(test_frame.Time - landmarks(clutter_landmark(i)).rbin_start*1.0e-6) < 2.5e-8,1);
    clu_y2 = find(abs(test_frame.Time - landmarks(clutter_landmark(i)).rbin_stop*1.0e-6) < 2.5e-8,1);

    clu_power = clu_power + mean(max(test_frame.Data([clu_y1:clu_y2],[clu_x1:clu_x2])));
  end
  
  clu_power = clu_power/num_clutter;
  
  qc_metric.peak_scr = 10*log10(sig_power/clu_power);

end

% This section computes the CNR
if isempty(clutter_landmark)
  fprintf('No clutter landmark for current filename and thus no CNR computed\n');
elseif noise_landmark == 0
  fprintf('No noise landmark for current filename and thus no CNR computed\n');
else
  
  if signal_landmark == 0 %This means that the clutter/noise power has not been computed yet
    
    % This section extracts the average clutter power from the various selected
    % clutter regions
    num_clutter = length(clutter_landmark);
    clu_power = 0;

    for i = 1:num_clutter
%     clu_x1 = round(landmarks(clutter_landmark(i)).rline_start);
%     clu_x2 = round(landmarks(clutter_landmark(i)).rline_stop);
      clu_x1 = find(abs(test_frame.GPS_time - landmarks(clutter_landmark(i)).gpstime_start) < 1.0,1);
      clu_x2 = find(abs(test_frame.GPS_time - landmarks(clutter_landmark(i)).gpstime_stop) < 1.0,1);

      clu_y1 = find(abs(test_frame.Time - landmarks(clutter_landmark(i)).rbin_start*1.0e-6) < 2.5e-8,1);
      clu_y2 = find(abs(test_frame.Time - landmarks(clutter_landmark(i)).rbin_stop*1.0e-6) < 2.5e-8,1);

      clu_power = clu_power + mean(max(test_frame.Data([clu_y1:clu_y2],[clu_x1:clu_x2])));
    end
    
    clu_power = clu_power/num_clutter;

    % This section extracts the average noise power from the selected noise region
%   noise_x1 = round(landmarks(noise_landmark).rline_start);
%   noise_x2 = round(landmarks(noise_landmark).rline_stop);
    noise_x1 = find(abs(test_frame.GPS_time - landmarks(noise_landmark).gpstime_start) < 1.0,1);
    noise_x2 = find(abs(test_frame.GPS_time - landmarks(noise_landmark).gpstime_stop) < 1.0,1);

    noise_y1 = find(abs(test_frame.Time - landmarks(noise_landmark).rbin_start*1.0e-6) < 2.5e-8,1);
    noise_y2 = find(abs(test_frame.Time - landmarks(noise_landmark).rbin_stop*1.0e-6) < 2.5e-8,1);

    noise_power = mean(mean(test_frame.Data([noise_y1:noise_y2],[noise_x1:noise_x2])));
    
  end
  
  qc_metric.peak_cnr = 10*log10(clu_power/noise_power);
end

% This section computes the Clutter Power
if signal_landmark == 0 && noise_landmark == 0 && ~isempty(clutter_landmark)

  % This section extracts the average clutter power from the various selected
  % clutter regions
  num_clutter = length(clutter_landmark);
  clu_power = 0;

  for i = 1:num_clutter
%     clu_x1 = round(landmarks(clutter_landmark(i)).rline_start);
%     clu_x2 = round(landmarks(clutter_landmark(i)).rline_stop);
    clu_x1 = find(abs(test_frame.GPS_time - landmarks(clutter_landmark(i)).gpstime_start) < 1.0,1);
    clu_x2 = find(abs(test_frame.GPS_time - landmarks(clutter_landmark(i)).gpstime_stop) < 1.0,1);

    clu_y1 = find(abs(test_frame.Time - landmarks(clutter_landmark(i)).rbin_start*1.0e-6) < 2.5e-8,1);
    clu_y2 = find(abs(test_frame.Time - landmarks(clutter_landmark(i)).rbin_stop*1.0e-6) < 2.5e-8,1);

    clu_power = clu_power + mean(max(test_frame.Data([clu_y1:clu_y2],[clu_x1:clu_x2])));
  end
  
  qc_metric.clu_power = 10*log10(clu_power/num_clutter);
  
end

% This section computes the Noise Power
if signal_landmark == 0 && isempty(clutter_landmark)

  % This section extracts the average noise power from the selected noise region
%   noise_x1 = round(landmarks(noise_landmark).rline_start);
%   noise_x2 = round(landmarks(noise_landmark).rline_stop);
    noise_x1 = find(abs(test_frame.GPS_time - landmarks(noise_landmark).gpstime_start) < 1.0,1);
    noise_x2 = find(abs(test_frame.GPS_time - landmarks(noise_landmark).gpstime_stop) < 1.0,1);

    noise_y1 = find(abs(test_frame.Time - landmarks(noise_landmark).rbin_start*1.0e-6) < 2.5e-8,1);
    noise_y2 = find(abs(test_frame.Time - landmarks(noise_landmark).rbin_stop*1.0e-6) < 2.5e-8,1);

    qc_metric.noise_power = 10*log10(mean(mean(test_frame.Data([noise_y1:noise_y2],[noise_x1:noise_x2]))));
    
end

% % This section is to compute the gradient of the test frame using spatial
% % filter
% 
% % Using Sobel kernel
% spatial_filter = fspecial('sobel');
% img_grad = abs(filter2(spatial_filter,test_frame.Data.*1e12,'same'));
% 
% % Using Laplacian of Gaussian kernel
% % sigma = 3;
% % img_grad = GradientMagnitude(test_frame.Data.*1e15,sigma);

figure(201)
imagesc(lp(test_frame.Data));
colormap(1-gray(256));
xlabel('Range line');
% ylabel('Two-way Propagation (us)');
ylabel('Range Bin Index');
title(sprintf('Test Data Frame ID: %s',img_fname_short(6:end)),'Interpreter','none');

% Drawing some annotations/overlays over the test frame image
if signal_landmark ~= 0
  rectangle('Position',[sig_x1 sig_y1 sig_x2-sig_x1+1 sig_y2-sig_y1+1],'EdgeColor',[1 0 0]); %Draw rectange over signal region
end

if ~isempty(clutter_landmark)
  for i = 1:num_clutter
%     clu_x1 = round(landmarks(clutter_landmark(i)).rline_start);
%     clu_x2 = round(landmarks(clutter_landmark(i)).rline_stop);
    clu_x1 = find(abs(test_frame.GPS_time - landmarks(clutter_landmark(i)).gpstime_start) < 1.0,1);
    clu_x2 = find(abs(test_frame.GPS_time - landmarks(clutter_landmark(i)).gpstime_stop) < 1.0,1);

    clu_y1 = find(abs(test_frame.Time - landmarks(clutter_landmark(i)).rbin_start*1.0e-6) < 2.5e-8,1);
    clu_y2 = find(abs(test_frame.Time - landmarks(clutter_landmark(i)).rbin_stop*1.0e-6) < 2.5e-8,1);

    rectangle('Position',[clu_x1 clu_y1 clu_x2-clu_x1+1 clu_y2-clu_y1+1],'EdgeColor',[0 1 0]); %Draw rectange over clutter region

  end
end

if noise_landmark ~= 0
  rectangle('Position',[noise_x1 noise_y1 noise_x2-noise_x1+1 noise_y2-noise_y1+1],'EdgeColor',[0 0 1]); %Draw rectange over noise region
end

return

% =====================================================================
% Example
% =====================================================================

fname = '/cresis/scratch2/mdce/mcords2/2011_Greenland_P3/CSARP_standard/20110507_01/Data_20110507_01_028.mat';
img_metric = qc_compute(fname);