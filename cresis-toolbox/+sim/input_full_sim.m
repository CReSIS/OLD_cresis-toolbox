hm;
param=[];
hdr=[];
data = [];
frames = [];
records = [];


%% INPUTS to sim.flightline_extract.m

param.sim.radar_name  = 'snow';
param.sim.season_name = '2012_Greenland_P3'; %'2018_Antarctica_DC8';
param.sim.day_seg     = '20120330_04'; % '20181010_01';

param.sim.radar_name  = 'snow';
param.sim.season_name = '2013_Greenland_P3';
param.sim.day_seg     = '20130327_02';

% param.radar_name  = 'snow';
% param.season_name = '2013_Greenland_P3';
% param.day_seg     = '20130327_02';

% param.sim.radar_name  = 'snow';
% param.sim.season_name = '2016_Greenland_Polar6';
% param.sim.day_seg     = '20160414_01';
%
% param.sim.radar_name  = 'snow';
% param.sim.season_name = '2018_Antarctica_DC8';
% param.sim.day_seg     = '20181010_01';
% %
% param.sim.radar_name  = 'rds';
% param.sim.season_name = '2014_Greenland_P3';
% param.sim.day_seg     = '20140325_05';

% if 1
%   param.sim.start_gps_time  = 1.333120368172158e+09;
%   param.sim.stop_gps_time   = 1.333120370976195e+09;
% else
%   param.sim.frame_idx = 2;
% end
if 0
  dd=load('/cresis/snfs1/dataproducts/ct_data/snow/2013_Greenland_P3/CSARP_post/CSARP_qlook/20130327_02/Data_20130327_02_002.mat');
  imagesc([],dd.Time,lp(dd.Data));
  % 4.68e-6
end

param.sim.monostatic = [1];

param.sim.presum_combine = 1;

param.sim.rx_combine = 0;
param.sim.rx_combine_weights = [1];

param.sim.wfs_combine = 0;
param.sim.wfs = 1;

param.sim.rxpath = 1;

param.sim.imgs = {[1 1]};% ,[2,1] {[ones([7 1]),(1:7).'], [2*ones([7 1]),(1:7).'], [3*ones([7 1]),(1:7).']};

% param.radar.lever_arm_fh = ''; % disabled for rds

c = physical_constants('c');

[ param, hdr, frames, records, exec_good ] = sim.flightline_extract(param);

% flightline_extract;
if ~exec_good; fprintf('flightline_extract executed incompletely\n'); return; end; %clear exec_good;

%% Load the following from the outputs of sim.flightline_extract
debug_plot_en  = 0;
data = [];
%     Nx = min(1000,length(param.gps.gps_time)); % temporary use
Nx = length(param.gps.gps_time); % rec_len

%point target for all rlines
for img = 1:length(param.sim.imgs)
  for wf_adc = 1:size(param.sim.imgs{img},1)
    
    wf = param.sim.imgs{img}(wf_adc,1);
    adc = param.sim.imgs{img}(wf_adc,2);
    wfs = param.signal.wfs(wf);
    %     phase_center = lever_arm(param, param.signal.wfs(wf).tx_weights, adc); %adc is rx_path(adc)?
    time        = wfs(wf).time;
    Tpd         = wfs(wf).Tpd;
    f0          = wfs(wf).f0;
    fc          = wfs(wf).fc;
    fs          = param.radar.fs;
    chirp_rate  = wfs(wf).chirp_rate;
    rx_tukey    = wfs(wf).tukey;
    t_NCO       = wfs(wf).DDC_NCO_delay;
    time_NCO    = time + t_NCO;
    
%     range = 660 * ones(Nx,1);
    range = sqrt( (param.gps.x - param.target.x).^2 + (param.gps.y - param.target.y).^2 + (param.gps.z - param.target.z).^2 );
    
    twtt = 2*range/c;
    
    param.sim.range = range;
    param.sim.twtt = twtt;
    
    if sum(hdr.DDC_dec{img}) == Nx; DDC_ON = 0; else DDC_ON = 1; end
    b = fir1(100,0.45,'low');
    
    Nt = hdr.Nt{img}(1);
    
    for rec = 1:Nx      %<#<#<#<#<#<#<#<#<#<#<#<#<#<#<#<#<#<#<#<#<#<
      td        = twtt(rec);
      fb        = chirp_rate*td;
      tmp_time  = time - td;
      
      
      if param.signal.wfs(wf).deramp  %deramp
        raw_data{wf,adc}(:,rec) = tukeywin_cont(tmp_time/Tpd-0.5,rx_tukey) .* cos(2*pi*f0*td + pi*chirp_rate*(2*time*td - td^2)); % Nt x Nrx;
        
        if DDC_ON   % Digital down conversion
          %           if mod(rec,2)
          offset = 4.858519900497512e+03/4;
          offset = 0;
%           hdr.DDC_freq{img}(rec) = round(hdr.DDC_freq{img}(rec) / offset)*offset;
          % 1.365710206292821e+08
          hdr.DDC_freq{img}(rec) = hdr.DDC_freq{img}(rec) + mod(rec,2)* offset; %%###########################
          
          %           dt=time_NCO(2)-time_NCO(1);
          DDC_data{wf,adc}(:,rec) = raw_data{wf,adc}(:,rec) .* exp(1i*2*pi*(-hdr.DDC_freq{img}(rec))*(time_NCO));
          
%           STEP(rec) = 1.1+0.5*mod(rec,2);
%           
%           f_STEP(rec) = STEP(rec)/Nt;
%           omega(rec) = 2*pi*f_STEP(rec);
%           new(:,rec) = exp(1i* omega(rec) * (0:Nt-1).');
%           
%           DDC_data{wf,adc}(:,rec) = raw_data{wf,adc}(:,rec) .* new(:,rec);
          
          %
                    fdata = fft(DDC_data{wf,adc}(:,rec));
          % %           fdata(1:1e4) = 0;
                    [max_val(rec),max_idx(rec)] = max(fdata);
          %           max_idx
          %           angle(max_val)
          %           [max_val,max_idx] = max(fdata(1e4+1:end));
          %           max_idx = max_idx + 1e4
          %           angle(max_val)
          %           k = round(STEP);
          %           %w0 = 2*pi/Nt * STEP;
          %           angle( 1-exp(j*w0*Nt) )
          %           angle( 1-exp(j*(w0-2*pi/Nt*k)) )
          %           angle( (1-exp(j*w0*Nt)) / (1-exp(j*(w0-2*pi/Nt*k))) )
          % %           k = Nt-k;
          % %           angle( (1-exp(j*w0*Nt)) / (1-exp(j*(w0-2*pi/Nt*k))) )
          
          if debug_plot_en
            freq = (-floor(Nt/2):floor((Nt-1)/2)) * fs/Nt;
            figure(23);clf(23);plot(freq/1e6,lp(fftshift(fft(raw_data{wf,adc}(:,rec))))); hold on; plot(freq/1e6,lp(fftshift(fft(DDC_data{wf,adc}(:,rec)))));grid on;
            title('IF vs NCO mixer output');
          end
          
          Nfilt = log2(hdr.DDC_dec{img}(rec));
          
          if debug_plot_en
            figure(97);hold off;plot(freq/1e6, fftshift(lp(fft(DDC_data{wf,adc}(:,rec)))));hold on;
            figure(98);hold off;plot(unwrap(angle((DDC_data{wf,adc}(:,rec)))));hold on;
            figure(99);hold off;plot(lp(DDC_data{wf,adc}(:,rec)));hold on;
          end
          
          prev_data = DDC_data{wf,adc}(:,rec);
          
          for idx = 1:Nfilt
            filt_data = filter(b,1,prev_data);
            
            if debug_plot_en
              tmp = [filt_data; zeros(Nt - size(filt_data,1),1)];
              figure(97);plot(freq/1e6,fftshift(lp(fft(tmp))));
              figure(98);plot(unwrap(angle((tmp))));
              figure(99);plot(lp(tmp));
            end
            
            dec_data = filt_data(1:2:end);
            
            if debug_plot_en
              tmp = [dec_data; zeros(Nt - size(dec_data,1),1)];
              figure(97);plot(freq/1e6,fftshift(lp(fft(tmp))));
              figure(98);plot(unwrap(angle((tmp))));
              figure(99);plot(lp(tmp));
            end
            
            prev_data = dec_data;
            
          end
          
          data{wf,adc}(:,rec) = [prev_data];%; zeros(Nt - size(prev_data,1),1)];

%           freq_axis = ifftshift(-floor(hdr.Nt{img}(rec)/2):floor((hdr.Nt{img}(rec)-1)/2)).';
%           if hdr.DDC_dec{img}(rec) == 2
%             Scale output
%             data{img}(1:hdr.Nt{img}(rec),rec,wf_adc) = ( 1/ hdr.DDC_dec{img}(rec) )...
%               * data{img}(1:hdr.Nt{img}(rec),rec,wf_adc);
%           elseif hdr.DDC_dec{img}(rec) == 4
%             Scale output
%             data{img}(1:hdr.Nt{img}(rec),rec,wf_adc) = ( 1/ (exp(1i*-pi/2)*hdr.DDC_dec{img}(rec)) ) ...
%               * data{img}(1:hdr.Nt{img}(rec),rec,wf_adc);
%             Delay of 100/4 = 25 relative to DDC_dec==2
%             extra_delay(rec) = 25;
%             data{img}(1:hdr.Nt{img}(rec),rec,wf_adc) = ifft(fft(data{img}(1:hdr.Nt{img}(rec),rec,wf_adc)) ...
%               ./ exp(1i*2*pi*25*freq_axis/hdr.Nt{img}(rec)));
%           elseif hdr.DDC_dec{img}(rec) == 8
%             Scale output
%             data{img}(1:hdr.Nt{img}(rec),rec,wf_adc) = (1 / (exp(1i*-pi/2)*hdr.DDC_dec{img}(rec)) ) ...
%               * data{img}(1:hdr.Nt{img}(rec),rec,wf_adc);
%             Delay of (100+2*100)/8 = 37.5 relative to DDC_dec==2
%             extra_delay(rec) = 37.5;
%             data{img}(1:hdr.Nt{img}(rec),rec,wf_adc) = ifft(fft(data{img}(1:hdr.Nt{img}(rec),rec,wf_adc)) ...
%               ./ exp(1i*2*pi*37.5*freq_axis/hdr.Nt{img}(rec)));
%           elseif hdr.DDC_dec{img}(rec) == 16
%             Scale output
%             data{img}(1:hdr.Nt{img}(rec),rec,wf_adc) = ( 1/ hdr.DDC_dec{img}(rec) )...
%               * data{img}(1:hdr.Nt{img}(rec),rec,wf_adc);
%             Delay of (100+2*100+4*100)/8 = 43.75 relative to DDC_dec==2
%             data{img}(1:hdr.Nt{img}(rec),rec,wf_adc) = ifft(fft(data{img}(1:hdr.Nt{img}(rec),rec,wf_adc)) ...
%               ./ exp(1i*2*pi*43.75*freq_axis/hdr.Nt{img}(rec)));
%           end
        else % NO DDC
          data{wf,adc}(:,rec) = raw_data{wf,adc}(:,rec);
        end
        
      else %pulsed
        data{img}(:,rec,wf_adc) = tukeywin_cont(tmp_time/Tpd-0.5,rx_tukey) * 0.5 .* exp(-1i*2*pi*f0*td + 1i*pi*chirp_rate*tmp_time.^2); % Nt x Nrx;
      end
      
      hdr.Nt{1}(rec) = size(data{1},1);
    end                   %<#<#<#<#<#<#<#<#<#<#<#<#<#<#<#<#<#<#<#<#<#<
    
  end
end



% clear dt t0 IF_filter_idx NCO DDC_filter_idx dec;
% Nt = length(time);
%
% hdr             = hdr;
% dt              = 1/param.radar.fs;
% t0              = time(1);
% IF_filter_idx   = zeros(1,Nx);
% NCO             = zeros(1,Nx);
% DDC_filter_idx  = zeros(1,Nx);
% dec             = ones(1,Nx);

param.records.file.version = 1000;

%%
fig_h = figure(123);clf(123);
subplot(221)
% try
imagesc(1:Nx,param.signal.wfs(1).time/1e-6,lp(data{1}));
% catch
%   imagesc(1:Nx,param.signal.wfs(1).time/1e-6,lp(data));
% end
title(param.sim.season_name,'Interpreter', 'none');
% hold on;
% plot(param.target.elev*2/3e8/1e-6,'x','LineWidth',2);
xlabel('rlines');
ylabel('fast-time, us');
subplot(222)
try
  plot(param.signal.wfs(1).time/1e-6,lp(data{1,1}(:,ceil(Nx/2))));
catch
  plot(param.signal.wfs(1).time/1e-6,lp(data(:,ceil(Nx/2))));
end
xlabel('fast-time, us');
ylabel('center rline');
title(param.sim.day_seg,'Interpreter', 'none');
axis tight;
subplot(223)
plot(1:Nx,range);
xlabel('rlines');
ylabel('range, m');
axis tight;
subplot(224)
plot(1:Nx,twtt/1e-6);
xlabel('rlines');
ylabel('twtt, us');
axis tight;

% figure(369);
% plot(1:Nx,twtt/1e-6);
% xlabel('rlines');
% ylabel('twtt, us');
% axis tight;

% fig_title = sprintf('%s %s',param.sim.season_name, param.sim.day_seg);
% print(fig_h, '-dpng', fig_title, '-r300');
% saveas(fig_h,fig_title);

% figure(124);clf(124);
% imagesc(1:Nx,param.signal.wfs.time/1e-6,lp(data_test{1,1}));
% hold on;
% plot(param.target.elev*2/3e8/1e-6,'x','LineWidth',2);
% xlabel('rlines');
% ylabel('fast-time, us');

%% Output_Files
for compress_this=1
  
  false_save_en = 0; % if 1; does not save the files
  
  fprintf('=====================================================================\n');
  fprintf('%s: %s (%s)\n', mfilename, param.sim.day_seg, datestr(now));
  fprintf('=====================================================================\n');
  
  % RAW DATA
  param.sim.out_fn_dir_raw_data = fullfile(gRadar.ct_tmp_path,'sim3D', ...
    sprintf('%s',param.sim.radar_name), ...
    sprintf('%s',param.season_name), ...
    sprintf('%s',param.day_seg(1:8)) );
  for img = 1:length(param.sim.imgs)
    for wf_adc = 1:size(param.sim.imgs{img},1)
      wf = param.sim.imgs{img}(wf_adc,1);
      adc = param.sim.imgs{img}(wf_adc,2);
      raw_data = data;
      param.sim.out_fn_raw_data = fullfile(param.sim.out_fn_dir_raw_data, ...
        sprintf('data_wfs_%02d_adc_%02d.mat',wf,adc) );
      fprintf('Saving raw_data %s (%s)\n', param.sim.out_fn_raw_data, datestr(now));
      if ~false_save_en; ct_save(param.sim.out_fn_raw_data, 'hdr', 'raw_data');end; %, 'dt', 't0', 'IF_filter_idx', 'NCO', 'DDC_filter_idx', 'dec'
    end
  end
  
  % FRAMES
  param.sim.out_fn_dir_frames = fullfile(gRadar.support_path,'frames', ...
    sprintf('%s',param.sim.radar_name), ...
    sprintf('%s',param.season_name) );
  param.sim.out_fn_frames = fullfile(param.sim.out_fn_dir_frames, ...
    sprintf('frames_%s.mat',param.day_seg) );
  fprintf('Saving frames %s (%s)\n', param.sim.out_fn_frames, datestr(now));
  if ~false_save_en; ct_save(param.sim.out_fn_frames, 'frames'); end;
  
  % RECORDS
  param.sim.out_fn_dir_records = fullfile(gRadar.support_path,'records', ...
    sprintf('%s',param.sim.radar_name), ...
    sprintf('%s',param.season_name) );
  param.sim.out_fn_records = fullfile(param.sim.out_fn_dir_records, ...
    sprintf('records_%s.mat',param.day_seg) );
  fprintf('Saving records %s (%s)\n', param.sim.out_fn_records, datestr(now));
  if ~false_save_en; ct_save(param.sim.out_fn_records,'-struct', 'records');
    if 0; fprintf('Creating auxiliary records files %s (%s)\n',param.sim.out_fn_records,datestr(now));
    create_records_aux_files(param.sim.out_fn_records);end; end;
  
  % PARAM
  param.sim.out_fn_dir_param = param.sim.out_fn_dir_raw_data;
  param.sim.out_fn_param = fullfile(param.sim.out_fn_dir_param, ...
    sprintf('param.mat') );
  param.fn = param.sim.out_fn_param;
  fprintf('Saving param %s (%s)\n', param.sim.out_fn_param, datestr(now));
  if ~false_save_en; ct_save(param.fn, 'param');end;
  
  if false_save_en; fprintf('===== FALSE_SAVE_EN: none of the above files are written =====\n');end;
  
end; clear compress_this;
%%

return;

% hdr.DDC_freq{img}(1:4)
% time_NCO(1:4)
% % x(t) = x(t) .* exp(j*2*pi*f*t)
% %
% % DDC_data{wf,adc}(:,rec) = raw_data{wf,adc}(:,rec) .* exp(1i*2*pi*(-hdr.DDC_freq{img}(rec))*(time_NCO));
%
% dt = time_NCO(2)-time_NCO(1);
% % idx = 15550
% idx = 0;
% f_data = bsxfun(@times,fft(data{1}(1:8175,:),[],1),exp(-1i*2*pi*(-hdr.DDC_freq{img})*dt*idx));
% [max_val, max_idx] = max(f_data,[],1);
%
% figure(678);
% clf;
% % plot(angle(max_val));
% diff(angle(max_val(1:2)))
% % hold on;

figure(1234);
% clf(1234);
subplot(211);plot(angle(max_val)*180/pi);ylabel('Degrees'); hold on;
subplot(212);plot(angle(max_val));ylabel('Radians'); hold on;

