% script rx_chan_equal_raw_analyze
%
% Reads outputs from rx_chan_equal_raw (CSARP_equal and CSARP_feedthru)
% and analyzes the results.
%
% Author: John Paden

params = read_param_xls(ct_filename_param('rds_param_2008_Greenland_TO.xls'),[],'equal');
% params = read_param_xls(ct_filename_param('rds_param_2009_Greenland_TO.xls'),[],'equal');


% feedthru_expected = [168.58	    -21.75	   -170.71	      0.00	   -167.48	     14.56]; % 2008 old and without correction
% feedthru_expected = [122.28	    -46.99	    167.64	      0.00	   -156.38	     58.29]; % 2009 old and without correction
% feedthru_expected = [-38.13	    146.93	     -0.00	   -172.38	     27.32	   -109.47]; % 2009 post-adc-slip correction
feedthru_expected = [-177.72	    -16.08	   -167.16	      0.00	   -166.31	      9.11]; % 2008 post-adc-slip correction
% feedthru_expected = [-93.62	     35.22	     -0.00	    110.45	    176.14	    -20.50]; % 2009 antennas rotated

load_equalization = true;
wf = 1;

auto_mode = false;

for param_idx = 1:length(params)
  param = params(param_idx);
  
  if ~isempty(regexpi(param.cmd.notes,'do not process')) || iscell(param.cmd.generic) || ~param.cmd.generic
    if auto_mode
      fprintf('%s\n', param.day_seg);
    end
    continue;
  end
  
  if ~auto_mode
    fprintf('%s %d %d %d\n', param.day_seg, param.equal.hack_mcrds_adc12, param.equal.hack_mcrds_adc34, param.equal.hack_mcrds_adc56);
  end
  
  % Load the frames file
  load(ct_filename_support(param,'','frames'))
  
  if load_equalization
    %% Load the CSARP_equal files (surface tracked)
    base_dir = ct_filename_out(param,'','CSARP_equal');
    peak_offset = [];
    peak_val = [];
    for frm = 1:length(frames.frame_idxs)
      tmp = load(sprintf(fullfile(base_dir,sprintf('Equal_%s_%03d_01.mat',param.day_seg,frm))));
       
      peak_offset = cat(2,peak_offset,tmp.peak_offset);
      peak_val = cat(2,peak_val,tmp.peak_val);
    end
  end
  param_equal = tmp.param_equal;
  if any(abs(tmp.param_equal.radar.wfs(1).Tsys - 1e-6*[0.541700000000000   0.536200000000000   0.539570000000000   0.540000000000000   0.539230000000000 0.539210000000000]) > 1e-12)
%     keyboard
  end
  
  %% Load the CSARP_feedthru files (feed through tracked)
  base_dir = ct_filename_out(param,'','CSARP_feedthru');
  ft_peak_offset = [];
  ft_peak_val = [];
  for frm = 1:length(frames.frame_idxs)
    tmp = load(sprintf(fullfile(base_dir,sprintf('Equal_%s_%03d_01.mat',param.day_seg,frm))));
    
    ft_peak_offset = cat(2,ft_peak_offset,tmp.peak_offset);
    ft_peak_val = cat(2,ft_peak_val,tmp.peak_val);
  end
  if any(abs(tmp.param_equal.radar.wfs(1).Tsys - 1e-6*[0.541700000000000   0.536200000000000   0.539570000000000   0.540000000000000   0.539230000000000 0.539210000000000]) > 1e-12)
%     keyboard
  end
  
  
  [B,A] = butter(4,0.01);
  ft_pv = filtfilt(B,A,(ft_peak_val .* conj(repmat(ft_peak_val(param_equal.equal.ref_wf_adc_idx,:),[size(ft_peak_val,1) 1]))).');
  ft_pv = ft_pv ./ repmat(exp(j*feedthru_expected/180*pi), [size(ft_pv,1) 1]);
  %ft_pv = ft_pv./repmat(mean(ft_pv,1), [size(ft_pv,1) 1]);
  if ~auto_mode
    fprintf('%20s\t', 'Feedthru Error:');
    fprintf('%10.2f\t', angle(mean(ft_pv,1))*180/pi)
    fprintf('\n');
    fprintf('%20s\t', 'Feedthru:');
    fprintf('%10.2f\t', angle(mean(ft_pv,1) .* exp(j*feedthru_expected/180*pi) )*180/pi)
    fprintf('\n');
  end
  
  if 0
    % HACK for manual tracking of phase
    [manual_track.x manual_track.y] = ginput; % Press return when done
    struct_to_matlab_cmds(manual_track)
  end
  if 0 && strcmpi('20090331_05',param.day_seg)
    % HACK for this segment because feed through was bad (manual tracking of phase required)
    manual_track(1).('x') = reshape(double([61.403508771929864 8485.9649122807004 9100 12833.333333333332 ]),[4 1]);
    manual_track(1).('y') = reshape(double([2.2727272727271952 0.4545454545453822 -53.181818181818244 -76.81818181818187 ]),[4 1]);
    channels = 1:2;
    for idx=1:length(channels)
      ft_pv(:,channels(idx)) = exp(j*interp1(manual_track.x,manual_track.y,1:size(ft_pv,1),'linear','extrap')/180*pi);
    end
  end
  
  figure(2); clf;
  h = plot(angle(ft_pv)*180/pi,'.');
  title(sprintf('%s', param.day_seg), 'interpreter', 'none');
  legend(h,{'1','2','3','4','5','6'})
  ylabel('Phase angle (deg)');
  xlabel('Estimate');
  grid on;
  
  hold on
    
  if load_equalization
    [B,A] = butter(4,0.05);
    
    pv = filtfilt(B,A,(peak_val .* conj(repmat(peak_val(param_equal.equal.ref_wf_adc_idx,:),[size(peak_val,1) 1]))).');
    pv = pv ./ repmat(exp(j*(param.radar.wfs(wf).chan_equal_deg - param_equal.radar.wfs(wf).chan_equal_deg)/180*pi), [size(pv,1) 1]);
    
    param.equal.hack_mcrds_adc12 = 0;
    param.equal.hack_mcrds_adc34 = 0;
    param.equal.hack_mcrds_adc56 = 0;
    
    chan_equal_deg = zeros(size(pv));
    if param.equal.hack_mcrds_adc12 == 0
      chan_equal_deg(:,1:2) = pv(:,1:2);
    elseif param.equal.hack_mcrds_adc12 == 1
      chan_equal_deg(:,1:2) = pv(:,1:2) ./ ft_pv(:,1:2);
      chan_equal_deg(:,1) = pv(:,1) ./ ft_pv(:,2);
      chan_equal_deg(:,2) = pv(:,2) ./ ft_pv(:,2);
    else
      chan_equal_deg(:,1:2) = pv(:,1:2) .* exp(-j*pi/2);
    end
    
    if param.equal.hack_mcrds_adc34 == 0
      chan_equal_deg(:,3:4) = pv(:,3:4);
    elseif param.equal.hack_mcrds_adc34 == 1
      chan_equal_deg(:,3:4) = pv(:,3:4) ./ ft_pv(:,3:4);
    else
      chan_equal_deg(:,3:4) = pv(:,3:4) .* exp(-j*pi/2);
    end
    
    if param.equal.hack_mcrds_adc56 == 0
      chan_equal_deg(:,5:6) = pv(:,5:6);
    elseif param.equal.hack_mcrds_adc56 == 1
      chan_equal_deg(:,5:6) = pv(:,5:6) ./ ft_pv(:,5:6);
      chan_equal_deg(:,5) = pv(:,5) ./ ft_pv(:,6);
      chan_equal_deg(:,6) = pv(:,6) ./ ft_pv(:,6);
    else
      chan_equal_deg(:,5:6) = pv(:,5:6) .* exp(-j*pi/2);
    end
    
    figure(1); clf;
    h = plot(angle(chan_equal_deg)*180/pi,'.');
    % h = plot(angle(pv./ft_pv)*180/pi,'.');
    % h = plot(angle(pv)*180/pi,'.');
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
    
%     pv = pv ./ ft_pv;
%     
%     mean_pv = mean(pv);
    
    %chan_equal_deg = angle(mean_pv)*180/pi + mean_without_outliers(angle(pv ./ repmat(mean_pv, [size(pv,1) 1]))*180/pi, 1);
    %chan_equal_deg = mean_without_outliers(angle(chan_equal_deg));
    chan_equal_deg = mean(chan_equal_deg);

    chan_equal_deg = chan_equal_deg / chan_equal_deg(param_equal.equal.ref_wf_adc_idx);
    chan_equal_deg = angle(chan_equal_deg) * 180/pi;
    if ~auto_mode
      fprintf('%20s\t', 'Error:');
      fprintf('%10.2f\t', chan_equal_deg)
      fprintf('\n');
    end
    fprintf('%-20s\t', param.day_seg);
    fprintf('%10.2f\t', angle(exp(j*(param.radar.wfs(wf).chan_equal_deg + chan_equal_deg)/180*pi))*180/pi)
    fprintf('\n');
    
  end
  
  if ~auto_mode
    pause
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

