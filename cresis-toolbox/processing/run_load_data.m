% Script run_load_data
%
% Several examples of using load_data
%
% Author: John Paden
%
% See also load_data, pulse_compress
% param = read_param_xls(ct_filename_param('snow_param_2012_Greenland_P3.xls'),'20120330_04');
param = read_param_xls(ct_filename_param('snow_param_2013_Greenland_P3.xls'),'20130420_01');
% param = read_param_xls(ct_filename_param('snow_param_2014_Greenland_P3.xls'),'20140428_01');

% Determine which records you want to load:
frames_fn = '';
frames_fn = ct_filename_support(param,frames_fn,'frames');
load(frames_fn);
frm = 271; % raw and DDC
% frm = 318; % raw
% frm = 319; % raw
% frm =56;
% frm = 317; % All raw
frm = 1;% All raw
% param.load_data.recs = [frames.frame_idxs(frm)-3000 frames.frame_idxs(frm+1)-1]; %+ 0 + [1 10000];
%   param.load_data.recs = [1 10000]; %+ 0 + [1 10000];
param.load_data.recs = frames.frame_idxs(frm) +[1 100];
% param.load_data.recs = [frames.frame_idxs(frm) frames.frame_idxs(frm+1)-1];
%   param.load_data.recs = 1.2e6 +[1 5000];

%   param.load_data.imgs = {[-1j 5]};
%   param.load_data.imgs = {[2 2; 2 3; 2 4; 2 5; 2 6; 2 7; 2 8; 2 9; 2 10; 2 11; 2 12; 2 13; 2 14; 2 15; 2 16]};
param.load_data.imgs                  = {[1 1]};
param.load_data.pulse_comp            = true;
param.load_data.raw_data              = false;
param.load_data.ft_wind               = @hanning;
param.load_data.combine_rx            = false;
param.radar.wfs(1).coh_noise_method = '';           %<#######
param.radar.wfs(1).coh_noise_method = 'analysis';           %<#######
%   param.radar.wfs(1).prepulse_H.type = 'inverse_filter';
%   param.radar.wfs(1).prepulse_H.fn = 'bottom';
%   param.radar.wfs(1).prepulse_H.smooth = @(H) fir_dec(H',hanning(123)'/sum(hanning(123)),1);
% Load data
[hdr,data] = load_data(param);

%%

plot_data = fir_dec(data{1},1);
% plot_data = fir_dec(data{1}(1:10:end,:),10);
[Nt,Nx] = size(plot_data);

figure(7); clf;
% imagesc(1:Nx,hdr.time{1}/1e-6,lp(plot_data));
imagesc(lp(plot_data));
figure(8); clf;
% imagesc(1:Nx,hdr.time{1}/1e-6,angle(plot_data));
imagesc(angle(plot_data));
link_figures([7 8]);

figure(5); clf(5);
subplot(211);plot(hdr.DDC_dec{1}); title('DDD\_dec');
subplot(212); plot(hdr.DDC_freq{1}/1e6); title('DDD\_freq, MHz');

figure(6); clf(6);
plot(hdr.nyquist_zone_hw{1},'g'); hold on;
plot(hdr.nyquist_zone_signal{1},'r');
legend({'hw','signal'});
title('nz');
link_axes([7 8 5 6],'x');

plot_data = fir_dec(data{1},10);
% plot_data = fir_dec(data{1}(1:10:end,:),10);
[Nt,Nx] = size(plot_data);

figure(3); clf;
% imagesc(1:Nx,hdr.time{1}/1e-6,lp(plot_data));
imagesc(lp(plot_data));
figure(4); clf;
% imagesc(1:Nx,hdr.time{1}/1e-6,angle(plot_data));
imagesc(angle(plot_data));
link_figures([3 4]);
return;
%%
if 0
  idx = 2580;
  real_data = data{1}(:,idx-1000:idx-1);
  idx = 1000;
  comp_data = data{1}(:,idx+1:idx+1000);
  
  real_mean =   nanmean(fir_dec(real_data,100),2);
  comp_mean = nanmean(fir_dec(comp_data,100),2) ;
  figure(994);clf(994);
  plot(180/pi*angle(real_mean));
  hold on;
  plot(180/pi*angle(comp_mean));
  plot(180/pi*angle(real_mean.*conj(comp_mean)),'.');
  legend({'real','complex','diff'});
  my_rbins = 3500:6000;
  180/pi*angle(nanmean(real_mean(my_rbins).*conj(comp_mean(my_rbins))))
  
  % idx = 2580;
  % real_data = data{1}(:,idx-1000:idx-1);
  % idx = 2600;
  % comp_data = data{1}(:,idx+1:idx+700);
  %
  % real_mean =   nanmean(fir_dec(real_data,100),2);
  % comp_mean = nanmean(fir_dec(comp_data,100),2) ;
  % figure(994);clf(994);
  % plot(180/pi*angle(real_mean));
  % hold on;
  % plot(180/pi*angle(comp_mean));
  % plot(180/pi*angle(real_mean.*conj(comp_mean)),'.');
  % legend({'real','complex','diff'});
  % my_rbins = 3500:6000;
  % 180/pi*angle(nanmean(real_mean(my_rbins).*conj(comp_mean(my_rbins))))
  
  
end

%%
if 0
  
  figure(1); clf;
  imagesc(1:round(Nx/100),hdr.time{1}/1e-6,lp(fir_dec(data{1},100)));
  figure(2); clf;
  imagesc(1:round(Nx/100),hdr.time{1}/1e-6,180/pi*angle(fir_dec(data{1},100)));
  link_figures([1 2]);
  
end

if 0
  data1 = plot_data;
  data1(isnan(data1)) = 0;
  figure(99);
  imagesc(lp(fft(data1)));
end
%   figure(100);clf(100);
%   line1 = data{1}(:,idx  );
%   line2 = data{1}(:,idx+1);
%   plot(lp(fft(line1))); hold on;
%   plot(lp(fft(line2)));
%   figure(101);clf(101);
%   plot(xcorr(line1,line2));

%%

plot_data = fir_dec(data{1},100);
[Nt,Nx] = size(plot_data);

figure(100);clf(100);
line_idx = 25;
line1 = nanmean(plot_data(:,line_idx+(-10:-1)),2);
line1(isnan(line1)) = 0;
% line1 = fft(line1);
line_idx = 58;
line2 = nanmean(plot_data(:,line_idx+(1:10)),2);
line2(isnan(line2)) = 0;
% line2 = fft(line2);

subplot(211)
plot(lp(line1)); hold on;
plot(lp(line2));
% plot(lp(conj(line1).*line2));
legend({'F(1)','F(2)'});
grid on;

subplot(212)
plot(lp(fft(conj(line1).*line2)));
grid on;
%   figure(3); clf;
%   aa=fir_dec(data{1},10);
%   plot(lp(aa([6114 6117 6118],:).'))
%

