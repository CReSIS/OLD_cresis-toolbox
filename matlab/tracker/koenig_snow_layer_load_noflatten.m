%% Test Results
[B,A] = butter(2,0.03);
figure(1000); clf;
cdata.Data = [];
cdata.GPS_time = [];
cdata.lay_twtt = [];
for frm = 231:253
  mdata = load_L1B(sprintf('/cresis/snfs1/dataproducts/ct_data/snow/2011_Greenland_P3/CSARP_post/CSARP_qlook/20110329_01/Data_20110329_01_%03d.mat',frm));
  dt = mdata.Time(2)-mdata.Time(1);
  lay = load(sprintf('/cresis/snfs1/dataproducts/ct_data/snow/2011_Greenland_P3/CSARP_layerData/20110329_01/Data_20110329_01_%03d.mat',frm));
  figure(1000);
%   clf;
  surf_twtt_filtered = filtfilt(B,A,lay.layerData{1}.value{2}.data);
  surf_twtt_filtered_bin = round(interp1(mdata.Time,1:length(mdata.Time),surf_twtt_filtered,'linear','extrap'));
  surf_twtt_filtered_bin = surf_twtt_filtered_bin-max(surf_twtt_filtered_bin);
  filter_correction = surf_twtt_filtered - lay.layerData{1}.value{2}.data;
%   for rline=1:size(mdata.Data,2)
%     mdata.Data(:,rline) = circshift(mdata.Data(:,rline),[-surf_twtt_filtered_bin(rline) 0]);
%   end
%   mdata.Data = fir_dec(mdata.Data,hanning(33).',1);
  filter_correction = filter_correction - dt*surf_twtt_filtered_bin;
%   for rline=1:size(mdata.Data,2)
%     mdata.Data(:,rline) = circshift(mdata.Data(:,rline),[surf_twtt_filtered_bin(rline) 0]);
%   end

  twtt_norm = mean(lay.layerData{1}.value{2}.data+filter_correction);

%   imagesc([], mdata.Time - twtt_norm, lp(mdata.Data))
%   hold on;
%   min_twtt = inf;
%   max_twtt = -inf;
%   h_plot = [];
%   for idx = 1:length(lay.layerData)
%     if idx==1
%       lay_twtt = lay.layerData{idx}.value{2}.data+filter_correction - twtt_norm;
%     else
%       lay_twtt = lay.layerData{idx}.value{2}.data+filter_correction - twtt_norm + 0.8e-9;
%     end
%     h_plot(end+1) = plot(lay_twtt);
%     min_twtt = min(min_twtt,nanmin(lay_twtt));
%     max_twtt = max(max_twtt,nanmax(lay_twtt));
%   end
%   if isfinite(min_twtt)
%     ylim([min_twtt max_twtt])
%   end
%   caxis([-10 8])
%   fprintf('%d: 1\n',frm);
%   pause
%   set(h_plot,'Visible','off')
%   pause
%   set(h_plot,'Visible','on')
%   pause
%   set(h_plot,'Visible','off')
%   pause
%   set(h_plot,'Visible','on')

if isempty(cdata.Data)
  cdata.Time = mdata.Time;
end
cdata.Data = cat(2,cdata.Data,interp1(mdata.Time,mdata.Data,cdata.Time));
cdata.GPS_time = cat(2,cdata.GPS_time,mdata.GPS_time);
lay_twtt = [];
for idx = 1:length(lay.layerData)
  if idx==1
    lay_twtt(idx,:) = lay.layerData{idx}.value{2}.data;
  else
    lay_twtt(idx,:) = lay.layerData{idx}.value{2}.data;
  end
end
cdata.lay_twtt = cat(2,cdata.lay_twtt,lay_twtt);

%   pause
end


clf
imagesc([],cdata.Time,lp(cdata.Data))
hold on;
% for idx = 1:length(cdata.lay_twtt)
  plot(cdata.lay_twtt.');
% end




