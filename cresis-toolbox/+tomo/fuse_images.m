function mdata = fuse_images(param)
% mdata = tomo.fuse_images(param)
%
% Description: Usually this function is called from tomo.collate_task.
%   Combines img_II data files into a single combined data file.
%
% Inputs:
%   param: struct from parameter spreadsheet
%    .tomo_collate
%     .imgs
%     .master_img_idx
%     .in_dir
%
% Outputs:
%   NONE
%
% See also: tomo.run_collate, tomo.collate, tomo_collate_task,
%   tomo.fuse_images, tomo.add_icemask_surfacedem, tomo.create_surfData,
%
% Author: John Paden, Jordan Sprick, and Mingze Xu

% imgs: list of images to fuse
imgs = param.tomo_collate.imgs;
% master_img_idx: all images will be resampled onto this image
master_img_idx = param.tomo_collate.master_img_idx;

% in_dir: Directory where 3D image files are at
in_dir = ct_filename_out(param,param.tomo_collate.in_dir);

% Load images
fns = {};
for img = imgs(:).'
  fns{img} = fullfile(in_dir,sprintf('Data_img_%02.0f_%s_%03.0f.mat',img, ...
    param.day_seg,param.proc.frm));
  mdata{img} = load(fns{img});
end

% Interpolate images onto a common propagation time axis
slave_imgs = imgs;
slave_imgs(master_img_idx) = [];
for img = slave_imgs
  mdata{img}.Topography.img = interp1(mdata{img}.Time,mdata{img}.Topography.img,mdata{master_img_idx}.Time);
end

% Fuse images
if length(imgs) ~= 3
  error('Currently tomo.fuse only supports no less or more than 3 images.');
end
Topography = mdata{master_img_idx}.Topography;
for rline = 1:size(mdata{1}.Topography.img,3)
%   if ~mod(rline-1,100)
%     fprintf('  Rline %d of %d (%s)\n', rline, size(mdata{1}.Topography.img,3), datestr(now));
%   end
%   Topography.img(:,:,rline) = reshape(tomo.fuse(double(10*log10(mdata{imgs(1)}.Topography.img(:,:,rline))), ...
%     double(10*log10(mdata{imgs(2)}.Topography.img(:,:,rline))), ...
%     double(10*log10(mdata{imgs(3)}.Topography.img(:,:,rline)))), [size(mdata{1}.Topography.img,1) size(mdata{1}.Topography.img,2)]);
  Topography.img(:,:,rline) = reshape(tomo.fuse(double(mdata{imgs(1)}.Topography.img(:,:,rline)), ...
    double(mdata{imgs(2)}.Topography.img(:,:,rline)), ...
    double(mdata{imgs(3)}.Topography.img(:,:,rline))), [size(mdata{1}.Topography.img,1) size(mdata{1}.Topography.img,2)]);
end
% Topography.img = 10.^(Topography.img/10);

% Copy master image file to combined image filename
combined_fn = fullfile(in_dir,sprintf('Data_%s_%03.0f.mat',param.day_seg,param.proc.frm));
copyfile(fns{master_img_idx}, combined_fn);

% Append new fused image to the combined file
save(combined_fn,'-append','Topography');

mdata = mdata{master_img_idx};
mdata.Topography = Topography;

end
