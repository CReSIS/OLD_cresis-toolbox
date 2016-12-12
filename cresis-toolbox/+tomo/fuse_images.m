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

if ~param.tomo_collate.vertical_fuse;
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
else
  
  N_t_bins = zeros(1,length(mdata));
  for i=1:length(mdata)
    N_t_bins(i) = length(mdata{i}.Time);
  end
  [~,sort_idx] = sort(N_t_bins);
  dummy = mdata;
  for i=1:length(mdata)
    mdata(i) = dummy(sort_idx(i));
  end
  clear dummy;
  
  % mean adjustment
%   mean_intensity = zeros(1,length(mdata));
  min_intensity = zeros(1,length(mdata));
  for i = 1:length(mdata)
%     mean_intensity(i) = mean2(lp(mdata{i}.Topography.img));
      min_intensity(i) = mean(min(min(mdata{i}.Topography.img)));
  end
  for i = 1:length(mdata)
%     mdata{i}.Topography.img = mdata{i}.Topography.img*10^((max(mean_intensity)-mean_intensity(i))/10);
%     mdata{i}.Topography.img = mdata{i}.Topography.img*10^((max(min_intensity)-min_intensity(i))/10);
    mdata{i}.Topography.img(mdata{i}.Topography.img<max(min_intensity)) = max(min_intensity);
  end
  
  Topography = mdata{master_img_idx}.Topography;
  
  rlines = size(mdata{1}.Topography.img,3);
  
  Surface = mdata{1}.Surface;
  
  img_comb = param.tomo_collate.img_comb;
  
  slave_imgs = imgs(imgs~=master_img_idx);
  
  time_min = inf;
  time_max = -inf;
  for img_idx = 1:length(mdata)
    time_min = min([time_min;mdata{img_idx}.Time]);
    time_max = max([time_max;mdata{img_idx}.Time]);
  end
  
  dt = mode(diff(mdata{master_img_idx}.Time));
  Time_fused = [flipud(-(-mdata{master_img_idx}.Time(1):dt:-time_min).') ...
    ; mdata{master_img_idx}.Time(2:end-1) ; (mdata{master_img_idx}.Time(end):dt:time_max).'];
  
  for img_idx = slave_imgs
    if mdata{img_idx}.Time(1) > Time_fused(1)
      start_time_idx = interp1(Time_fused,1:length(Time_fused),mdata{img_idx}.Time(1),'nearest');
    else
      start_time_idx = 1;
    end
    if start_time_idx+length(mdata{img_idx}.Time)-1 > length(Time_fused)
      start_time_idx = start_time_idx-1;
    end
    mdata{img_idx}.Time = Time_fused(start_time_idx+(0:(length(mdata{img_idx}.Time)-1)));
  end
  
  img_combined = nan(length(Time_fused),size(mdata{master_img_idx}.Topography.img,2) ...
    ,size(mdata{master_img_idx}.Topography.img,3));
  
    first_img_bins = zeros(rlines,length(img_comb));
    first_comb_bins = zeros(rlines,length(img_comb)/2);
    last_img_bins = zeros(rlines,length(img_comb));
    for comb = 1:length(img_comb)
      img = floor(comb/2)+1;
      f = comb + mod(comb,2) - 1;
      l = comb + mod(comb,2);
      first_twtt = Surface+img_comb(f);
      if ~mod(comb,2)
        first_twtt(first_twtt>max(mdata{img-1}.Time)) = max(mdata{img-1}.Time);
      end
      first_img_bins(:,comb) = interp1(mdata{img}.Time,1:length(mdata{img}.Time),first_twtt,'nearest');
      if mod(comb,2)
        first_comb_bins(:,ceil(comb/2)) = interp1(Time_fused,1:length(Time_fused),Surface+img_comb(f),'nearest');
      elseif ~all(first_comb_bins(:,comb/2) == interp1(Time_fused,1:length(Time_fused),Surface+img_comb(f),'nearest').')
        error('');
      end
      if isfinite(img_comb(l));
        last_img_bins(:,comb) = interp1(mdata{img}.Time,1:length(mdata{img}.Time),Surface+img_comb(l),'nearest');
        last_img_bins((Surface+img_comb(l))>max(mdata{img}.Time),comb)=max(mdata{img}.Time);
      else
        last_img_bins(:,comb) = interp1(mdata{img}.Time,1:length(mdata{img}.Time),mdata{img+mod(comb,2)-1}.Time(end),'nearest');
      end
    end
  
  for rline = 1:rlines
      
    img_combined(1:first_comb_bins(rline,1),:,rline) = ...
      mdata{1}.Topography.img(1:first_img_bins(rline,1),:,rline);
%     for theta_idx = 1:size(mdata{1}.Topography.img,2);
      
      for comb = 1:length(img_comb)/2;

        N_bin = last_img_bins(rline,comb*2+[-1,0])-first_img_bins(rline,comb*2+[-1,0])+1;
        if diff(N_bin)
          error('');
        end
        N_bin = N_bin(1);
        
        weight = (.5*cos((0:N_bin-1)/N_bin*pi)+.5).';
        
        curr_imgs = comb+[0,1];
        
        img_combined(first_comb_bins(rline,comb)+(0:N_bin-1),:,rline) = ...
          bsxfun(@times,weight,mdata{curr_imgs(1)}.Topography.img( ...
          first_img_bins(rline,2*comb-1)+(0:N_bin-1),:,rline)) ...
          + bsxfun(@times,flipud(weight),mdata{curr_imgs(2)}.Topography.img( ...
          first_img_bins(rline,2*comb)+(0:N_bin-1),:,rline));
        
        img_combined(round(first_comb_bins(rline,comb)+N_bin-1)+(1:length(mdata{curr_imgs(2)}.Time)-last_img_bins(rline,2*comb)),:,rline) = ...
          mdata{curr_imgs(2)}.Topography.img((last_img_bins(rline,2*comb)+1):end,:,rline);

      end
%     end
  end
  Topography.img = img_combined;
  Time = Time_fused;
end

% Copy master image file to combined image filename
% combined_fn = fullfile(in_dir,sprintf('Data_%s_%03.0f.mat',param.day_seg,param.proc.frm));
combined_fn = fullfile(in_dir,sprintf('Data_%s_%03.0f.mat',param.day_seg,param.proc.frm));
copyfile(fns{master_img_idx}, combined_fn);

% Append new fused image to the combined file
save(combined_fn,'-append','Topography');

mdata = mdata{master_img_idx};
mdata.Topography = Topography;
if exist('Time_fused','var')
  mdata.Time = Time_fused;
  save(combined_fn,'-append','Time');
end

end
