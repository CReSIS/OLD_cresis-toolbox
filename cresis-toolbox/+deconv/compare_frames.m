% script deconv.compare_frames

%% User Settings

fn = '20120314_02/Data_20120314_02_061.mat';
first_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2012_Greenland_P3/CSARP_post/CSARP_deconv/';
second_dir = '/cresis/snfs1/dataproducts/ct_data/snow/2012_Greenland_P3/CSARP_post/CSARP_qlook/';

rlines = 407;
xlim_rng = [-360 50];

%% Automated Section

mdata = load(fullfile(first_dir,fn));
figure(1); clf;
imagesc(lp(mdata.Data));
colormap(1-gray(256));

fprintf('Deconvolution waveforms:\n');
unique(mdata.custom.deconv_filter_idx)
figure(2);clf;
plot(mdata.custom.deconv_filter_idx)

mdata2 = load(fullfile(second_dir,fn));
figure(3); clf;
imagesc(lp(mdata2.Data));
colormap(1-gray(256));


for rline = rlines
  figure(4); clf;
  [peak_val,peak_idx] = max(lp(mdata.Data(:,rline)));
  [peak_val,peak_idx2] = max(lp(mdata2.Data(:,rline)));
  h_plot = plot(lp(mdata.Data(:,rline)));
  xlim(peak_idx + xlim_rng);
  ylim(peak_val + [-50 3]);
  grid on;
  hold on;
  h_plot(2) = plot(peak_idx-peak_idx2+(1:size(mdata2.Data,1)), lp(mdata2.Data(:,rline)),'r')
  xlim(peak_idx + xlim_rng);
  ylim(peak_val + [-50 3]);
  xlabel('Range bins');
  ylabel('Relative power (dB)')
    
  legend(h_plot,'Deconvolution','No Deconvolution','location','best');
%  legend(h_plot,'first','second','location','best');
  
  [~,fn_name] = fileparts(fn);
  frm_str = fn_name(end-2:end)
  title(sprintf('%s_%s\n', mdata.param_get_heights.day_seg,frm_str),'Interpreter','none')
  if rline ~= rlines(end)
    keyboard
  end
end
