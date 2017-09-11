function mdata = fuse_images(param)
% mdata = tomo.fuse_images(param)
%
% Description: Usually this function is called from tomo_collate_task.
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
  
  % normalized intensity of images
  min_intensity = zeros(1,length(mdata));
  for i = 1:length(mdata)
    min_intensity(i) = mean(min(min(mdata{i}.Topography.img)));
  end
  for i = 1:length(mdata)
    mdata{i}.Topography.img(mdata{i}.Topography.img<max(min_intensity)) = max(min_intensity);
  end
  
  Topography = mdata{master_img_idx}.Topography;
  
  rlines = size(mdata{1}.Topography.img,3);
  
  Surface = mdata{1}.Surface;
  
  img_comb = param.tomo_collate.img_comb;
  
  slave_imgs = imgs(imgs~=master_img_idx);
  
  % gets min and max twtts of imgs for sorting
  time_min = mdata{1}.Time(1);
  time_max = -inf;
  for img_idx = 1:length(mdata)
    time_max = max([time_max,mdata{img_idx}.Time(end)]);
  end
  
  % images universal time increment
  dt = mode(diff(mdata{master_img_idx}.Time));
  % time axis of combined image product
  Time_fused = [flipud(-(-mdata{master_img_idx}.Time(1):dt:-time_min).') ...
    ; mdata{master_img_idx}.Time(2:end-1) ; (mdata{master_img_idx}.Time(end):dt:time_max).'];
  
  % interpolates images time axis to combined image time axis
  for slave_img_idx = slave_imgs
    start_time_idx = interp1(Time_fused,1:length(Time_fused),mdata{slave_img_idx}.Time,'nearest');
    if isnan(start_time_idx(1))
      start_time_idx = find(~isnan(start_time_idx),1);
      range = 1:(length(mdata{slave_img_idx}.Time)-start_time_idx);
      mdata{slave_img_idx}.Time = Time_fused(range);
      mdata{slave_img_idx}.Topography.val = mdata{slave_img_idx}.Topography.val(range,:,:);
      mdata{slave_img_idx}.Topography.freq = mdata{slave_img_idx}.Topography.freq(range,:,:);
      mdata{slave_img_idx}.Topography.img = mdata{slave_img_idx}.Topography.img(range,:,:);
    else
      start_time_idx = start_time_idx(1);
      mdata{slave_img_idx}.Time = Time_fused(start_time_idx + 0:(length(mdata{slave_img_idx}.Time)-1));
    end
  end
  
  % combined image product
  img_combined = nan(length(Time_fused),size(mdata{master_img_idx}.Topography.img,2) ...
    ,size(mdata{master_img_idx}.Topography.img,3));
  
    % beginning index of fuse in images
    first_img_bins = zeros(rlines,length(img_comb));
    % last index of fuse in images
    last_img_bins = zeros(rlines,length(img_comb));
    % beginning index of fuse in combined image
    first_comb_bins = zeros(rlines,length(img_comb)/2);
    
    % loops thorough img_comb to determine fuse range of each fuse
    for comb = 1:length(img_comb)
      % target image of current img_comb components
      img = floor(comb/2)+1;
      % img_comb component of fuse start
      f = comb + mod(comb,2) - 1;
      % img_comb component of fuse end
      l = comb + mod(comb,2);
      
      % twtt associated with fuse start img_comb component for each slice
      first_twtt = Surface+img_comb(f);
      % twtt associated with fuse start img_comb component for each slice
      last_twtt = Surface+img_comb(l);
      
      % ???
      if ~mod(comb,2)
        first_twtt(first_twtt>max(mdata{img-1}.Time)) = max(mdata{img-1}.Time);
        if any(first_twtt>max(mdata{img-1}.Time))
          keyboard;
        end
      end
      
      % gets image index of fuse start
      first_img_bins(:,comb) = interp1(mdata{img}.Time,1:length(mdata{img}.Time),first_twtt,'nearest');
      
      % Gets fused image index of fuse start (if statement ensure it only
      % runs once for each combination (on first loop of each combo)
      if mod(comb,2)
        first_comb_bins(:,ceil(comb/2)) = interp1(Time_fused,1:length(Time_fused),first_twtt,'nearest');
      end
      
      % gets image index of fuse end
      if isfinite(img_comb(l));
        last_img_bins(:,comb) = interp1(mdata{img}.Time,1:length(mdata{img}.Time),last_twtt,'nearest');
        % ensures last index does not exceed length of image
        last_img_bins(last_twtt>max(mdata{img}.Time),comb) = max(mdata{img}.Time);
      else
        % uses last twtt of first image in combo
        last_img_bins(:,comb) = interp1(mdata{img}.Time,1:length(mdata{img}.Time),min(mdata{ceil(comb/2)}.Time(end),mdata{img}.Time(end)),'nearest');
      end
    end
  
  for rline = 1:rlines
    
    % inserts section of first image not in any fuse
    img_combined(1:first_comb_bins(rline,1),:,rline) = ...
      mdata{1}.Topography.img(1:first_img_bins(rline,1),:,rline);
      
      % loops through fuse combinations
      for comb = 1:length(img_comb)/2;

        N_bins = last_img_bins(rline,comb*2+[-1,0])-first_img_bins(rline,comb*2+[-1,0])+1;
        % Number of time bins included in the current image combination
        N_bin = min(N_bins);
        
        % weight used in image combination
        weight = (.5*cos((0:N_bin-1)/N_bin*pi)+.5).';
        
        % images included in the current image combination
        curr_imgs = comb+[0,1];
        
        % combines weighted images inside the fuse area
        img_combined(first_comb_bins(rline,comb)+(0:N_bin-1),:,rline) = ...
          bsxfun(@times,weight,mdata{curr_imgs(1)}.Topography.img( ...
          first_img_bins(rline,2*comb-1)+(0:N_bin-1),:,rline)) ...
          + bsxfun(@times,flipud(weight),mdata{curr_imgs(2)}.Topography.img( ...
          first_img_bins(rline,2*comb)+(0:N_bin-1),:,rline));
        
        % appends portion of image 2 in combination after the fuse area
        img_combined(first_comb_bins(rline,comb)+(N_bin:(length(mdata{curr_imgs(2)}.Time)-first_img_bins(rline,2*comb))),:,rline) = ...
          mdata{curr_imgs(2)}.Topography.img((first_img_bins(rline,2*comb)+N_bins):end,:,rline);
        
        % appends portion of image 1 that might be extend beyond above
        % portion of image 2
        if max(mdata{curr_imgs(2)}.Time)<max(mdata{curr_imgs(1)}.Time)
          img_combined(round(first_comb_bins(rline,comb))+((length(mdata{curr_imgs(2)}.Time)-first_img_bins(rline,2*comb)+1):(length(mdata{curr_imgs(1)}.Time)-first_img_bins(rline,2*comb-1))),:,rline) = ...
            mdata{curr_imgs(1)}.Topography.img((length(mdata{curr_imgs(2)}.Time)+1):end,:,rline);
        end

      end
      if any(any(isnan(img_combined(:,:,rline))))
        error('Image pixels were excluded from fuse.');
      end
  end
  Topography.img = img_combined;
  Time = Time_fused;
end

% Copy master image file to combined image filename
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
