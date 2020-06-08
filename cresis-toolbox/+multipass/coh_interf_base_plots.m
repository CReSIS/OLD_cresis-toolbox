clearvars -except gRadar tmp
close all
%File names for the different years
fn2011_2014 = '/cresis/snfs1/dataproducts/ct_data/rds/2011_Greenland_P3/CSARP_insar/rds_thule_2011_2014_wf2_insar.mat';
fn2011_2012 = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/rds_thule_2011_2012_wf2_insar3.mat';
fn2012_2014 = '/cresis/snfs1/dataproducts/ct_data/rds/2012_Greenland_P3/CSARP_insar/rds_thule_2012_2014_wf2_insar.mat';
fn2012_2013 = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/rds_thule_2012_2013_wf2_insar3.mat';
fn2013_2014 = '/cresis/snfs1/dataproducts/ct_data/rds/2013_Greenland_P3/CSARP_insar/rds_thule_2013_2014_wf2_insar.mat';
fn2014_2014 = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/rds_thule_2014_2Week_wf3_insar3.mat';
fn20142 = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/rds_thule_2014_1Day_wf3_insar3.mat';
fn2014allwf = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_multipass/rds_thule_2014_SameDay_allwf_multipass03.mat';

fns = {fn2014allwf};
figtitle = {'2014 Same Day'};
saven = {'2014_SameDay_allwf'};
masterid = 8;
%Settings 
set1 = {1:2:30};%{1:15,1:15, 1:15, 1:15, 1:7, 1:15,1:15,1:15};
set2 = {2:2:30};%{16:30,16:30, 16:30, 16:30, 8:22, 16:22, 16:30, 16:30};
if 0 %7 Element Setting
  set1 = {5:11, 5:11, 5:11, 5:11, 5:11};
  set2 = {20:26, 20:26, 16:22, 20:26, 20:26};
  saven  = cellfun(@(c)[c '_7Elements'],saven,'UniformOutput',false);
end
rbins = [240:420];
filt_x = 9; filt_y = 45;
%Convert phase to physical length
eps_ice = 3.15; wfuse = 2;
cls = physconst('LightSpeed');
savedir = '/cresis/snfs1/scratch/bmiller/More_Multipass/Check_Figures/';

for fn_id = 1:length(fns)
  if ~exist('tmp','var') || length(tmp)<fn_id
    tmp{fn_id} = load(fns{fn_id});
  end
  insar1 = mean(tmp{fn_id}.data(rbins,:,set1{fn_id}),3);
  insar2 = mean(tmp{fn_id}.data(rbins,:,set2{fn_id}),3);
  insar_data = insar2 .* conj(insar1) ./ (abs(insar1).*abs(insar2));
  %Block averaging (smooths image)
%   insar_data_filtrow = fir_dec(transpose(insar_data),ones(1,filt_y)/filt_y,1);
%   insar_data_filt = fir_dec(transpose(insar_data_filtrow),ones(1,filt_x)/filt_x);
  insar_data_filt = fir_dec(fir_dec(insar_data,ones(1,filt_y)/filt_y,1).',ones(1,filt_x)/filt_x).';
  %Normalization (makes phase begin at 0)
  if 1%~strcmp(saven{fn_id}(1:4),'2014')
    insar_data_filt = bsxfun(@times,insar_data_filt, ...
      exp(-1i*angle(insar_data_filt(round(size(insar_data_filt,1)/2),:))));
  end
  %Convert geodetic to along track
  alongx = geodetic_to_along_track(tmp{fn_id}.ref.lat,tmp{fn_id}.ref.lon,tmp{fn_id}.ref.elev)/1000;%km
  %Convert rbins to depth
  BW = tmp{fn_id}.ref.wfs(wfuse).f1-tmp{fn_id}.ref.wfs(wfuse).f0;
  propspeed = (cls/(2*sqrt(eps_ice)));
  crossy = rbins/BW*propspeed-mean(tmp{fn_id}.pass(1).surface)*cls/2;
  %Generate the title string
  title_str = figtitle{fn_id};
  
  %% Plot the interferogram
  figure(fn_id + 30); clf;
  imagesc(alongx,crossy,hsv_plot(insar_data_filt,-7));
  colormap(hsv(256))
  h_colorbar = colorbar;
  caxis([-pi pi])
  set(get(h_colorbar,'ylabel'),'string','Angle (radians)');
  grid on;
  xlabel('Along Track (km)');
  ylabel('Depth (m)');
  title(sprintf('%s Interferogram',title_str));
  grid on
  set(get(gcf,'Children'),'Fontsize',14)
  set(findall(gcf,'type','line'),'linewidth',2)
  savename = fullfile(savedir,sprintf('%s_Interferogram',saven{fn_id}));
  saveas(gcf,[savename '.fig'])
  saveas(gcf,[savename '.png'])
  %% Plot the coherence
  figure(fn_id+40); clf;
  imagesc(alongx,crossy,abs(insar_data_filt));
  h_colorbar = colorbar;
  set(get(h_colorbar,'ylabel'),'string','Coherence');
  xlabel('Along Track (km)');
  ylabel('Depth (m)');
  title(sprintf('%s Coherence',title_str));
  grid on
  set(get(gcf,'Children'),'Fontsize',14)
  set(findall(gcf,'type','line'),'linewidth',2)
  savename = fullfile(savedir,sprintf('%s_Coherence',saven{fn_id}));
  saveas(gcf,[savename '.fig'])
  saveas(gcf,[savename '.png'])
  
  %% Plot the absolute baseline
  absbase1 = []; absbase2 = [];
  for ab1 = set1{fn_id}
    absbase1(ab1,:) = sqrt(tmp{fn_id}.pass(ab1).ref_y.^2+tmp{fn_id}.pass(ab1).ref_z.^2);
  end
  for ab2 = set2{fn_id}
    absbase2(ab2,:) = sqrt(tmp{fn_id}.pass(ab2).ref_y.^2+tmp{fn_id}.pass(ab2).ref_z.^2);
  end
  %Plot the first two for the legend purposes
  figure(fn_id+50); clf;
  plot(alongx,absbase1(1,:),'r')
  hold on
  plot(alongx,absbase2,'b')
  plot(alongx,absbase1(2:end,:),'r')
  hold off
  leg = {sprintf('%s Elements',tmp{fn_id}.pass(set1{fn_id}(1)).param_multipass.day_seg(1:4))};
  leg{end+1} = sprintf('%s Elements',tmp{fn_id}.pass(set2{fn_id}(1)).param_multipass.day_seg(1:4));
  legend(leg)
  title(sprintf('%s Absolute Baselines',title_str))
  ylabel('Baseline to Reference Element (m)')
  ylim([0 60])
  xlabel('Along Track (km)')
  grid on
  set(get(gcf,'Children'),'Fontsize',14)
  set(findall(gcf,'type','line'),'linewidth',2)
  savename = fullfile(savedir,sprintf('%s_Baselines',saven{fn_id}));
  saveas(gcf,[savename '.fig'])
  saveas(gcf,[savename '.png'])
end

