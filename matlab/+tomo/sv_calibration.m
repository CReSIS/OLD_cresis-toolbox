
if 0
  fn = '/cresis/snfs1/dataproducts/ct_data/rds/2018_Greenland_P3/CSARP_music_imgs4_Nsig2/20180406_01/Data_20180406_01_004.mat';
  mdata=load(fn);
end
demv_threshold = 13;

Nt = size(mdata.Topography.img,1);
Nsv = size(mdata.Topography.img,2);
Nx = size(mdata.Topography.img,3);
nd = floor(Nsv/2)+1;
dem_bin = round(interp1(mdata.Time,1:length(mdata.Time),mdata.twtt(nd,:)));
dem_bin2 = min(Nt,round(interp1(mdata.Time,1:length(mdata.Time),mdata.twtt(nd,:)*2)));

tb = interp1(mdata.Time, 1:length(mdata.Time), mdata.twtt);
tb = interp_finite(tb);
demi = nan(size(tb));
demv = nan(size(tb));
for rline=1:Nx
  rline
  dtb = tb(:,rline) - tb(nd,rline);
  doa = (1:Nsv)-nd;
  rbins = dem_bin(rline) : dem_bin2(rline);
  
  % Right side
  doa_bins = nd:Nsv;
  [tb_sort,sort_idxs] = unique(tb(doa_bins,rline));
  surf_doa_bin = round(interp1(tb_sort, doa_bins(sort_idxs), rbins));
  mv = zeros(length(rbins),1);
  mi = zeros(length(rbins),1);
  for rbin_idx = 1:length(rbins)
    search_rng = round([-8-2*rbin_idx/length(rbins) 2-4*rbin_idx/length(rbins)]);
    rbin = rbins(rbin_idx);
    doa_bins = max(nd,surf_doa_bin(rbin_idx)+search_rng(1)) : min(Nsv,surf_doa_bin(rbin_idx)+search_rng(2));
    [mv(rbin_idx),mi(rbin_idx)] = max(mdata.Topography.img(rbin,doa_bins,rline),[],2);
    mi(rbin_idx) = mi(rbin_idx)+doa_bins(1);
  end
  demi(nd:end,rline) = interp1(rbins-rbins(1),mi-1,dtb(nd:end));
  demv(nd:end,rline) = interp1(rbins-rbins(1),mv,dtb(nd:end));
  
  % Left side
  doa_bins = 1:nd;
  [tb_sort,sort_idxs] = unique(tb(doa_bins,rline));
  surf_doa_bin = round(interp1(tb_sort, doa_bins(sort_idxs), rbins));
  nmv = zeros(length(rbins),1);
  nmi = zeros(length(rbins),1);
  for rbin_idx = 1:length(rbins)
    search_rng = round([0+2*rbin_idx/length(rbins)  8+2*rbin_idx/length(rbins)]);
    rbin = rbins(rbin_idx);
    doa_bins = max(1,surf_doa_bin(rbin_idx)+search_rng(1)) : min(nd,surf_doa_bin(rbin_idx)+search_rng(2));
    [nmv(rbin_idx),nmi(rbin_idx)] = max(mdata.Topography.img(rbin,doa_bins,rline),[],2);
    nmi(rbin_idx) = nmi(rbin_idx)+doa_bins(1);
  end
  demi(1:nd,rline) = interp1(rbins-rbins(1),nmi-1,dtb(1:nd));
  demv(1:nd,rline) = interp1(rbins-rbins(1),nmv,dtb(1:nd));
  if 1
    clf;
    plot(demi(:,rline))
    hold on;
    plot(1:Nsv,1:Nsv);
    grid on;
    figure(1); clf;
    imagesc(lp(mdata.Topography.img(:,:,rline)));
    ylim(rbins([1 end]));
    hold on;
    plot(tb(:,rline),'-x','LineWidth',2,'MarkerSize',10);
    demi_thresh = demi(:,rline);
    demi_thresh(lp(demv(:,rline))<demv_threshold) = NaN;
    demi_thresh(abs(demi_thresh-nd)<3) = NaN;
    plot(demi_thresh,tb(:,rline),'LineWidth',2);
    figure(2); clf;
    plot(lp(demv(:,rline)));
    grid on;
    keyboard
  end
end


demi_thresh = demi;
% demi_thresh(demv<10) = NaN;
% demi_thresh(abs(demi_thresh-nd)<3,:) = NaN;
% demi_thresh(1:9,:) = NaN;
demi_thresh(:,1000:end) = NaN;
% demi_thresh(:,1:25) = NaN;
figure(1); clf;
sv_error = bsxfun(@minus,demi_thresh,(1:Nsv).');
imagesc(sv_error)
caxis([-12 12]);
colorbar
% return

di = mean_without_outliers(demi_thresh.').';
% di(nd-5:nd+6) = NaN;
figure(2); clf;
plot(di);
hold on;
plot(1:Nsv,1:Nsv);
grid on;

figure(1); clf;
plot(sv_error,'.')
hold on;
plot(di - (1:Nsv).','x');
grid on;

sv_mean = di - (1:Nsv).';
sv_mean(nd + (-10:10)) = NaN;
begin_cut = 11;
end_cut = 10;
sv_mean(1:begin_cut) = NaN;
sv_mean(end-end_cut+2 : end) = NaN;
sv_mean = interp_finite(sv_mean);
plot(sv_mean,'k','LineWidth',2)

pp_order = 5;
doa_bins = 1+begin_cut:Nsv-end_cut;
pp = polyfit( doa_bins.', sv_mean(doa_bins),pp_order);
sv_mean_poly = nan(size(sv_mean));
sv_mean_poly(doa_bins) = polyval(pp,doa_bins);
sv_mean_poly = interp_finite(sv_mean_poly);
plot(sv_mean_poly,'m','LineWidth',2)

param = mdata.param_array;
if isfield(mdata,'theta')
  theta = mdata.theta;
else
  theta = mdata.param_array.array_param.theta;
end
theta_original = mdata.param_array.array_param.theta;

% theta_cal_fn = [ct_filename_ct_tmp(param,'','tomo_collate','theta_cal') '.mat'];
if ispc
  theta_cal_fn = 'E:\ct_tmp\tomo_collate\rds\2018_Greenland_P3\theta_cal_20180406_01.mat';
else
  theta_cal_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/tomo_collate/rds/2018_Greenland_P3/theta_cal_20180406_01.mat';
end
theta_cal.theta_original = theta_original;
theta_cal.theta = interp1(1:Nsv,theta,(1:Nsv).' + sv_mean_poly);
theta_cal_fn_dir = fileparts(theta_cal_fn);
if ~exist(theta_cal_fn_dir,'dir')
  mkdir(theta_cal_fn_dir);
end
save(theta_cal_fn,'-v7.3','-struct','theta_cal','theta','theta_original');

return;


