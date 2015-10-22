% script rx_chan_equal_raw_load
%
% Script for loading the output of rx_chan_equal_raw.m

% params = read_param_xls(ct_filename_param('rds_param_2006_Greenland_TO.xls'),'20060530_01'); params(1).cmd.generic = 1;
% params = read_param_xls(ct_filename_param('rds_param_2009_Greenland_TO.xls'),'20090331_05'); params(1).cmd.generic = 1;
params = read_param_xls(ct_filename_param('rds_param_2008_Greenland_TO.xls'),[],'equal');
% params = read_param_xls(ct_filename_param('rds_param_2009_Greenland_TO.xls'),[],'equal');

feedthru_expected = [168.58	    -21.75	   -170.71	      0.00	   -167.48	     14.56]; % 2008
% feedthru_expected = [122.28	    -46.99	    167.64	      0.00	   -156.38	     58.29]; % 2009
% feedthru_expected = [-93.62	     35.22	     -0.00	    110.45	    176.14	    -20.50]; % 2009 antennas rotated

for param_idx = 1:length(params)
  param = params(param_idx);
  
  if ~isempty(regexpi(param.cmd.notes,'do not process')) || iscell(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  %   wf_adc_pairs = 12:15;
  %   ref_wf_adc_idx = 4;
  
  fprintf('%s\n', param.day_seg);
  
  base_dir = ct_filename_out(param,'','CSARP_equal');
  load(ct_filename_support(param,'','frames'))
  
  peak_offset = [];
  peak_val = [];
  gps_time = [];
  for frm = 1:length(frames.frame_idxs)
    tmp = load(sprintf(fullfile(base_dir,sprintf('Equal_%s_%03d_01.mat',param.day_seg,frm))));
    
    peak_offset = cat(2,peak_offset,tmp.peak_offset);
    peak_val = cat(2,peak_val,tmp.peak_val);
    gps_time = cat(2,gps_time,tmp.gps_time);
  end
  param_equal = tmp.param_equal;
  ref_wf_adc_idx = tmp.param_equal.equal.ref_wf_adc_idx;
  wf_adc_pairs = 1:size(tmp.param_equal.equal.imgs{1},1);
  peak_offset = peak_offset(wf_adc_pairs,:);
  peak_val = peak_val(wf_adc_pairs,:);
  
  if 0
    ft_peak_offset = zeros(size(peak_offset));
    ft_peak_val = ones(size(peak_val));
  else
    ft_peak_offset = [];
    ft_peak_val = [];
    
    base_dir = ct_filename_out(param,'','CSARP_feedthru');
    for frm = 1:length(frames.frame_idxs)
      tmp = load(sprintf(fullfile(base_dir,sprintf('Equal_%s_%03d_01.mat',param.day_seg,frm))));
      
      ft_peak_offset = cat(2,ft_peak_offset,tmp.peak_offset);
      ft_peak_val = cat(2,ft_peak_val,tmp.peak_val);
    end
    ft_peak_offset = ft_peak_offset(wf_adc_pairs,:);
    ft_peak_val = ft_peak_val(wf_adc_pairs,:);
  end
    
  [B,A] = butter(4,0.01);
  ft_pv = filtfilt(B,A,(ft_peak_val .* conj(repmat(ft_peak_val(ref_wf_adc_idx,:),[size(ft_peak_val,1) 1]))).');
  ft_pv = ft_pv ./ repmat(exp(j*feedthru_expected/180*pi), [size(ft_pv,1) 1]);
  % ft_pv = ft_pv./repmat(mean(ft_pv,1), [size(ft_pv,1) 1]);
  % ft_pv = exp(j*angle(ft_pv));
  
  if strcmpi('20090331_05',param.day_seg)
    % HACK for this segment because feed through was bad
    manual_track(1).('x') = reshape(double([61.403508771929864 8485.9649122807004 9100 12833.333333333332 ]),[4 1]);
    manual_track(1).('y') = reshape(double([2.2727272727271952 0.4545454545453822 -53.181818181818244 -76.81818181818187 ]),[4 1]);
    channels = 1:2;
    for idx=1:length(channels)
      ft_pv(:,channels(idx)) = exp(j*interp1(manual_track.x,manual_track.y,1:size(ft_pv,1),'linear','extrap')/180*pi);
    end
  end
  
  figure(2); clf;
  h = plot(angle(ft_pv)*180/pi,'.');
  legend(h,{'1','2','3','4','5','6'})
  xlabel('Measurement #');
  ylabel('Drift angle (deg)');
  grid on;
  
  if 1
    [B,A] = butter(4,0.05);
    pv = filtfilt(B,A,(peak_val .* conj(repmat(peak_val(ref_wf_adc_idx,:),[size(peak_val,1) 1]))).');
    
    figure(1); clf;
    h = plot(angle(pv./ft_pv)*180/pi,'.');
    %   h = plot(angle(pv)*180/pi,'.');
    legend(h,mat2cell(char(48+(1:length(h))),1,ones(size(h))))
    grid on;
    % hold on;
    % h = plot(angle(ft_pv)*180/pi,'.');
    % hold off;
    
    
    % Feedthru to calculate change in phase over time
    %   - If we take the mean of the feedthru and apply that to the surface
    %     equalization data in the same way we apply it during processing,
    %     then the phase values from surface equalization will be correct
    % Absolute phase correction
    chan_equal_deg_relative = angle(ft_pv)*180/pi;
    
%     param.equal.hack_mcrds_adc12 = 0;
%     param.equal.hack_mcrds_adc34 = 0;
%     param.equal.hack_mcrds_adc56 = 0;
    
    chan_equal_deg = zeros(size(pv));
    if param.equal.hack_mcrds_adc12 ~= 1
      chan_equal_deg(:,1:2) = pv(:,1:2);
      ft_pv(:,1:2) = 1;
    elseif param.equal.hack_mcrds_adc12 == 1
      chan_equal_deg(:,1:2) = pv(:,1:2) ./ ft_pv(:,1:2);
      chan_equal_deg(:,1) = pv(:,1) ./ ft_pv(:,2);
      chan_equal_deg(:,2) = pv(:,2) ./ ft_pv(:,2);
      ft_pv(:,1) = ft_pv(:,2);
    end
    
    if param.equal.hack_mcrds_adc34 ~= 1
      chan_equal_deg(:,3:4) = pv(:,3:4);
      ft_pv(:,3:4) = 1;
    elseif param.equal.hack_mcrds_adc34 == 1
      chan_equal_deg(:,3:4) = pv(:,3:4) ./ ft_pv(:,3:4);
    end
    
    if param.equal.hack_mcrds_adc56 ~= 1
      chan_equal_deg(:,5:6) = pv(:,5:6);
      ft_pv(:,5:6) = 1;
    elseif param.equal.hack_mcrds_adc56 == 1
      chan_equal_deg(:,5:6) = pv(:,5:6) ./ ft_pv(:,5:6);
      chan_equal_deg(:,5) = pv(:,5) ./ ft_pv(:,6);
      chan_equal_deg(:,6) = pv(:,6) ./ ft_pv(:,6);
      ft_pv(:,5) = ft_pv(:,6);
    end
    
    pv = chan_equal_deg;
    
    mean_pv = mean(pv);
    
    chan_equal_dB = 10*log10(mean(abs(pv).^2));
    
    fprintf('%-20s', 'Power Error (dB)');
    fprintf('%.2f ', round(chan_equal_dB*100)/100);
    fprintf('\n');
    
    fprintf('%-20s ', 'chan_equal_dB');
    for idx = 1:length(wf_adc_pairs)%1:size(param_equal.equal.imgs{1},1)
      wf_adc_idx = wf_adc_pairs(idx);
      wf = param_equal.equal.imgs{1}(wf_adc_idx,1);
      adc = param_equal.equal.imgs{1}(wf_adc_idx,2);
      rx_path = param_equal.radar.wfs(wf).rx_paths(adc);
      fprintf('%.2f ', round((param_equal.radar.wfs(wf).chan_equal_dB(rx_path) + chan_equal_dB(idx))*100)/100)
    end
    fprintf('\n');
    
    chan_equal_deg = angle(mean_pv)*180/pi + mean_without_outliers(angle(pv ./ repmat(mean_pv, [size(pv,1) 1]))*180/pi, 1);
    
    fprintf('%-20s', 'Angle Error (deg)');
    fprintf('%.2f ', round(chan_equal_deg*100)/100);
    fprintf('\n');
    
    fprintf('%-20s ', 'chan_equal_deg');
    for idx = 1:length(wf_adc_pairs)%1:size(param_equal.equal.imgs{1},1)
      wf_adc_idx = wf_adc_pairs(idx);
      wf = param_equal.equal.imgs{1}(wf_adc_idx,1);
      adc = param_equal.equal.imgs{1}(wf_adc_idx,2);
      rx_path = param_equal.radar.wfs(wf).rx_paths(adc);
      fprintf('%.2f ', round((param_equal.radar.wfs(wf).chan_equal_deg(rx_path) + chan_equal_deg(idx))*100)/100)
    end
    fprintf('\n');
    
    fprintf('%-20s', 'Tsys (ns)');
    dt = tmp.wfs(wf).time(2)-tmp.wfs(wf).time(1);
    Tsys = mean_without_outliers(peak_offset.',1) * dt;
    fprintf('%.2f ', round(Tsys*1e9*100)/100);
    fprintf('\n');
    
  end
  
  %   pause
  if 1
    records_fn = ct_filename_support(param,'','records');
    records = load(records_fn);
    [B,A] = butter(4,0.01);
    adc_phase_corr_deg = interp1(gps_time,filtfilt(B,A,angle(ft_pv)*180/pi),records.gps_time,'linear');
    plot(adc_phase_corr_deg)
    pause;
    for adc = 1:size(adc_phase_corr_deg,2)
      adc_phase_corr_deg(:,adc) = interp_finite(adc_phase_corr_deg(:,adc),NaN);
    end
    save(records_fn,'-APPEND','adc_phase_corr_deg');
    create_records_aux_files(records_fn);
  end

end


return




[B,A] = butter(4,0.05);
figure(1); clf;
h = plot(angle(filtfilt(B,A,(peak_val ./ repmat(peak_val(param_equal.equal.ref_wf_adc_idx,:),[size(peak_val,1) 1])).'))*180/pi,'.');
legend(h,{'1','2','3','4','5','6'})
grid on;

%   hold on;
%   h = plot(angle(filtfilt(B,A,(peak_val_corr ./ repmat(peak_val_corr(param_equal.equal.ref_wf_adc_idx,:),[size(peak_val,1) 1])).'))*180/pi,'.');
%   hold off;
%   set(h,'Marker','x');
%   legend(h,{'1','2','3','4','5','6'})
%   grid on;

figure(2); clf;
h = plot(medfilt1(peak_offset.',11));
legend(h,{'1','2','3','4','5','6'})
grid on;

wf = param_equal.equal.imgs{1}(1);
time = wfs(wf).time;
dt = time(2)-time(1);

fprintf('%.2f ',round((param_equal.radar.wfs(wf).Tsys*1e9 + mean(medfilt1(peak_offset.',11))*dt*1e9)*100)/100)
fprintf('\n');

fprintf('%.2f ', round((param_equal.radar.wfs(wf).chan_equal_deg + angle(mean(peak_val.',1))*180/pi)*100)/100)
fprintf('\n');

