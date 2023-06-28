% fns = get_filenames('E:\tmp\HFSounder\bagley_ground\','stacked_combined','','.mat');
% caxis_lims = [-110 0];
% fns = get_filenames('E:\tmp\HFSounder\malaspina\','stacked_combined','','.mat');
%   caxis_lims = [-80 30];
fns = get_filenames('E:\tmp\HFSounder\wrangells_clutter\','stacked_combined','','.mat');
  caxis_lims = [-80 30];

physical_constants;

if 0
  % Create reference waveform
  % E:\tmp\HFSounder\malaspina\stacked_combined_5.mat
  sig = data(300:700,227);
  save('malaspina_sig.mat','sig')
end
load('malaspina_sig.mat')

for fn_idx = 1:length(fns)
  fn = fns{fn_idx};
  [fn_dir,fn_name] = fileparts(fn);
  [~,fn_dir_name] = fileparts(fn_dir);
  load(fn);
  
  good_mask = rec.lat~=0;
  rec.lat = rec.lat(good_mask);
  rec.lon = rec.lon(good_mask);
  rec.ele = rec.ele(good_mask);
  rec.ch0 = rec.ch0(good_mask,:);
  rec.ch1 = rec.ch1(good_mask,:);
  
  data = rec.ch0.';
  data = data(353:end,:); % Estimate of where zero time occurs
  dt = rec.dt;
  Nt = size(data,1);
  time = dt * (0:Nt-1).';
  range = time * c/2 / sqrt(er_ice);
  T = dt*Nt;
  df = 1/T;
  freq = df * (0:Nt-1).';
  
  along_track = geodetic_to_along_track(rec.lat,rec.lon,rec.ele);
  
  % Coherent noise removal... not effective (detrending instead)
  %   if 1
  %     data = bsxfun(@minus,data,mean(data,2));
  %   end
  
  % Skip small images
  if size(data,2) < 20
    continue
  end
  
  % Detrend data
  if 0
    noise_mean = mean(mean(data(4000:6000,:)));
    data = data-noise_mean;
  elseif 1
    data = data - sgolayfilt(data,3,101);
  else
    for rline=1:size(data,2)
      data(:,rline) = data(:,rline) - detrend(data(1200:2000,rline));
      %       data(:,rline) = data(:,rline) - mean(data(4000:6000,rline));
    end
  end
  
  %   blank_mask = time < 3.93e-6;
  %   data(blank_mask,:) = 0;
  
  % Matched filter using cross correlation
  new_data=[];
  for rline = 1:size(data,2)
    new_data(:,rline) = xcorr(data(:,rline),sig);
  end
  data= new_data;
  data = data(end-Nt+1-80:end-80,:);
  
  % Down convert and then low pass filter data
  for rline = 1:size(data,2)
    data(:,rline) = data(:,rline).*exp(-j*2*pi*4e6*time);
  end
  [B,A] = butter(2,2/50);
  data = filtfilt(B,A,data);
  
  % 2D image filtering
  % data = medfilt2(abs(data).^2,[120,3]);
  data = medfilt2(abs(data).^2,[17,3]);
  
  figure(1); clf;
  imagesc([], freq/1e6, lp(fft(data)))
  
  figure(2); clf;
  %   imagesc(lp(data))
  %   imagesc([], time*1e6, lp(data))
  imagesc([], range, lp(data))
  ylim([0 2200]);
  colormap(1-gray(256));
  ylabel('Range, {\epsilon}_{ice}=3.15, (m)')
  xlabel('Range line')
  title([fn_dir_name '/' fn_name],'interpreter','none')
  caxis(caxis_lims)
  
  %   pause
  %     return
  
  last_underscore_idx = find(fn_name=='_',1,'last');
  new_fn = fullfile(fn_dir,sprintf('view%s.mat', fn_name(last_underscore_idx:end)));
  save(new_fn,'data','range');
  new_fn_fig = fullfile(fn_dir,sprintf('view%s.fig', fn_name(last_underscore_idx:end)));
  saveas(2,new_fn_fig);
  new_fn_fig = fullfile(fn_dir,sprintf('view%s.jpg', fn_name(last_underscore_idx:end)));
  saveas(2,new_fn_fig);
end


