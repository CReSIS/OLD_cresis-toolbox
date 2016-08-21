function mdata_combined = surface_extractor(param,mdata)
% mdata_combined = surface_extractor(param,mdata)
%
% Description: Usually this function is called from tomo.collate. Tracks
%  surface on fused frame data.
%
% Inputs:
%   param = struct with processing parameters
%   mdata = contains frame data
%
% Outputs:
%   mdata_combined = contains fused frame data
%
% See also: tomo.collate
%
% Author: John Paden, Jordan Sprick, and Mingze Xu

nadir_img = 2;

% fn_dir: Directory where 3D image files are at
fn_dir = ct_filename_out(param,param.surf_extract.out_dir);

fprintf('Fusing Images...\n');
% Check for twtt DEM data
if ~isfield(mdata{1},'twtt')
  Nx = length(mdata{1}.GPS_time);
  Nsv = length(mdata{1}.param_combine.array_param.theta);
  twtt = NaN*zeros(Nsv,Nx);
end

% Interpolate images onto a common propagation time axis
for img  = [1 3]
  mdata{img}.Data = interp1(mdata{img}.Time,mdata{img}.Data,mdata{nadir_img}.Time);
  mdata{img}.Topography.img = interp1(mdata{img}.Time,mdata{img}.Topography.img,mdata{nadir_img}.Time);
  mdata{img}.Time = mdata{nadir_img}.Time;
end

% Convert Surface and Bottom variables from propagation time into image pixels
Surface_bin = interp1(mdata{1}.Time, 1:length(mdata{1}.Time), mdata{1}.Surface);
Surface_Mult_bin = interp1(mdata{1}.Time, 1:length(mdata{1}.Time), 2 * mdata{1}.Surface);
Bottom_bin = interp1(mdata{1}.Time, 1:length(mdata{1}.Time), mdata{1}.Bottom);
twtt_bin = interp1(mdata{1}.Time, 1:length(mdata{1}.Time), mdata{1}.twtt);
ice_mask = mdata{1}.ice_mask;
Bottom_bin(isnan(Bottom_bin)) = -1;
twtt_bin(isnan(twtt_bin)) = -1;

%% Automated labeling section
% =========================================================================
% Specify which range lines to browse
skip = 1;
rlines = 1:skip:size(mdata{1}.Topography.img,3);

% Template size
ms = 11;

% Training parameters and preparing fusion slices
mu = [];
sigma = [];
fusion_slices = [];
tic;
for rline = rlines
  fprintf('rline: %d (%f sec)\n', rline, toc);
  fusion = tomo.fuse(double(10*log10(mdata{1}.Topography.img(:,:,rline))), ...
    double(10*log10(mdata{2}.Topography.img(:,:,rline))), ...
    double(10*log10(mdata{3}.Topography.img(:,:,rline))));
  fusion = reshape(fusion, size(10*log10(mdata{1}.Topography.img(:,:,rline))));
  fusion_slices(:,:,rline) = fusion;
  fusion(fusion>13.5) = 13.5;
  [m, s] = tomo.train_params(fusion, double(twtt_bin(:,rline)));
  mu = [mu; m];
  sigma = [sigma; s];
end

mdata_combined = mdata{nadir_img};
mdata_combined.Topography.img = fusion_slices;
mdata_combined.Topography.mu = mu;
mdata_combined.Topography.sigma = sigma;
mdata_combined.ice_mask = ice_mask;

fn_combined = fullfile(fn_dir,sprintf('Data_%s_%03.0f.mat', ...
  param.day_seg,mdata{1}.frm));
save(fn_combined,'-struct','mdata_combined');

mdata_combined.frm = mdata{1}.frm;

%% Automated cycling section
% =========================================================================
% 3D surface
bottom_surface = [];

fprintf('Tracking Layer...\n');
for rline = rlines
  fprintf('rline: %d (%f sec)\n', rline, toc);
  detect_data = fusion_slices(:,:,rline);
  detect_data(detect_data>13.5) = 13.5;
  labels = tomo.detect(detect_data, double(twtt_bin(:,rline)), ...
    double(Bottom_bin(rline)), [], double(ice_mask(:,rline)), double(mean(mu)), double(mean(sigma)));
  
  % If empty, assign random values
  if isempty(labels)
    bottom_surface = [bottom_surface; randi([1 700], 1, size(bottom_surface, 2))];
  else
    bottom_surface = [bottom_surface; labels];
  end
  
  %     set(h_image(1),'CData',fusion_slices(:,:,rline));
  %     set(h_title(1),'String',sprintf('%d',rline));
  %     set(h_surf_plot(1),'YData',Surface_bin(rline))
  %     set(h_bot_plot(1),'YData',Bottom_bin(rline))
  %     set(h_surf_mult_plot(1),'YData',Surface_Mult_bin(rline))
  %     set(h_dem_plot(1),'YData',twtt_bin(:,rline))
  %     set(h_ice_bed_plot(1),'YData',labels)
  %     ylim(Surface_bin(rline)+[-25 500]);
  %     set(hplot,'XData',[rline rline]);
  
  % keyboard
end

save(fn_combined,'-append','bottom_surface');
mdata_combined.bottom_surface = bottom_surface;

return