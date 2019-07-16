if 0
  fn = '/cresis/snfs1/dataproducts/ct_data/rds/2011_Greenland_P3/CSARP_insar/rds_thule_2011_2014_insar.mat';
  set1 = 1:15;
  set2 = 16:30;
  title_str = '2011 to 2014';
  
elseif 1
  fn = '/cresis/snfs1/dataproducts/ct_data/rds/2012_Greenland_P3/CSARP_insar/rds_thule_2012_2014_wf1_insar.mat';
  set1 = 1:15;
  set2 = 16:30;
  title_str = '2012 to 2014';
elseif 1
  fn = '/cresis/snfs1/dataproducts/ct_data/rds/2013_Greenland_P3/CSARP_insar/rds_thule_2013_2014_wf1_insar.mat';
  set1 = 1:15;
  set2 = 16:22;
  title_str = '2013 to 2014';
end
rbins = [240:420];

tmp = load(fn);

insar1 = mean(tmp.data(rbins,:,set1),3);
insar2 = mean(tmp.data(rbins,:,set2),3);
insar_data = insar2 .* conj(insar1) ./ (abs(insar1).*abs(insar2));

insar_data_filt = fir_dec(fir_dec(insar_data,ones(1,31)/31,1).',ones(1,3)/3).';

figure(3001); clf;
imagesc(hsv_plot(insar_data_filt,-10));
colormap(hsv(256))
h_colorbar = colorbar;
caxis([-pi pi])
set(get(h_colorbar,'ylabel'),'string','Angle (radians)');
grid on;
xlabel('Range line');
ylabel('Range bin');
title(title_str);
set(3001,'Position',[56   470   867   426]);

figure(3000); clf;
imagesc(abs(insar_data_filt));
h_colorbar = colorbar;
set(get(h_colorbar,'ylabel'),'string','Coherence');
grid on;
xlabel('Range line');
ylabel('Range bin');
title(title_str);
set(3000,'Position',[56   470   867   426]);

[fn_dir,fn_name] = fileparts(fn);
fn_insar_phase = fullfile(fn_dir,[fn_name '_phase.fig']);
saveas(3001,fn_insar_phase);

[fn_dir,fn_name] = fileparts(fn);
fn_insar_coherence = fullfile(fn_dir,[fn_name '_coherence.fig']);
saveas(3000,fn_insar_coherence);


