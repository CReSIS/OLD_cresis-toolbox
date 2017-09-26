function mdata = fuse_images(param)
% mdata = tomo.fuse_images(param)
%
% Description: Usually this function is called from tomo_collate_task.
%   Combines Data_img_II_YYYYMMDD_SS_FFF.mat data files into a single combined
%   Data_YYYYMMDD_SS_FFF.mat file.
%
% Inputs:
%   param: struct from parameter spreadsheet
%
% Outputs:
%   NONE
%
% See also: tomo.run_collate, tomo.collate, tomo_collate_task,
%   tomo.fuse_images, tomo.add_icemask_surfacedem, tomo.create_surfData,
%
% Author: John Paden, Jordan Sprick, and Mingze Xu

% Nimgs: Number of images to fuse
Nimg = length(param.tomo_collate.imgs);

% in_dir: Directory where 3D image files are at
in_dir = ct_filename_out(param,param.tomo_collate.in_dir);

% Load images
fns = {};
for img_idx = 1:Nimg
  img = param.tomo_collate.imgs(img_idx);
  fns{img_idx} = fullfile(in_dir,sprintf('Data_img_%02.0f_%s_%03.0f.mat',img, ...
    param.day_seg,param.proc.frm));
  mdata{img_idx} = load(fns{img_idx});
end

% Copy first image file to combined "fused" image filename (this is a simple way
% to copy all the parameters and support variables over). We will update
% the 3D image and time fields in the file at the end of fusing.
combined_fn = fullfile(in_dir,sprintf('Data_%s_%03.0f.mat',param.day_seg,param.proc.frm));
copyfile(fns{1}, combined_fn);

Topography = [];
if strcmpi(param.tomo_collate.fuse_method,'horizontal')
  %% Horizontal Fuse
  % Interpolate images onto a common propagation time axis

  % dt: Fast time sample period 
  dt = mdata{1}.Time(2)-mdata{1}.Time(1);

  % Time: common time axis between all images
  start_time = -inf;
  stop_time = inf;
  for img_idx = 1:Nimg
    if mdata{img_idx}.Time(1) > start_time
      start_time = mdata{img_idx}.Time(1);
    end
    if mdata{img_idx}.Time(end) < stop_time
      stop_time = mdata{img_idx}.Time(end);
    end
  end
  Time = start_time:dt:stop_time;

  % Reinterpolate images onto common time axis
  Data = interp1(mdata{1}.Time,mdata{1}.Data,Time);
  for img_idx = 1:Nimg
    mdata{img_idx}.Topography.img = interp1(mdata{img_idx}.Time,mdata{img_idx}.Topography.img,Time);
  end
  
  % Nx: number of range lines
  Nx = size(mdata{1}.Topography.img,3);
  % Nsv: number of steering vector/directions of arrival
  Nsv = size(mdata{1}.Topography.img,2);

  % Create fuse weights
  fuse_weights = [];
  Topography.img = zeros(size(mdata{1}.Topography.img));
  for img_idx = 1:Nimg
    fuse_weights(img_idx,:) = exp( -((0:Nsv-1)-Nsv*((Nimg-img_idx)*2+1)/(2*Nimg)).^2 / (2*(Nsv/Nimg).^1.5) );
    Topography.img = Topography.img + bsxfun(@times,fuse_weights(img_idx,:), ...
      mdata{img_idx}.Topography.img);
  end
  Topography.img = bsxfun(@times,1./sum(fuse_weights),Topography.img);
  
  % Append new fused image to the combined file
  save(combined_fn,'-append','Time','Data','Topography');
  
  mdata = mdata{1};
  mdata.Time = Time;
  mdata.Data = Data;
  mdata.Topography.img = Topography.img;
  
elseif strcmpi(param.tomo_collate.fuse_method,'vertical')
  %% Vertical Fuse: Combine images into a single image (also trim time<0 values)
  Nx = size(mdata{1}.Topography.img,3);
  
  % Get the step size from the first image (all images should have the same step size)
  dt = mdata{1}.Time(2)-mdata{1}.Time(1);
  
  % Load each image and then combine with previous image
  for img = 1:length(mdata)
    if img == 1
      first_idx = find(mdata{img}.Time <= 0,1,'last');
      if ~isempty(first_idx)
        Time = mdata{img}.Time(first_idx:end);
        Data = mdata{img}.Data(first_idx:end,:,:);
        Topography.img = mdata{img}.Topography.img(first_idx:end,:,:);
      end
    else
      % Combine images
      % Data,Time => already loaded data
      % append.Data, append.Time => new data to append
      % New_Time, New_Data => Combined result
      
      % Interpolate image N onto already loaded data (assumption is that image
      % N-1 always comes before image N)
      New_Time = (Time(1) : dt : mdata{img}.Time(end)).';
      mdata{img}.Data = interp1(mdata{img}.Time,mdata{img}.Data,New_Time,'linear',0);
      mdata{img}.Topography.img = interp1(mdata{img}.Time,mdata{img}.Topography.img,New_Time,'linear',0);
      
      % Surface tracking image combine
      %  param.tomo_collate.img_comb(1): time after surface return where
      %    combine will happen
      %  param.tomo_collate.img_comb(2): minimum time that combine will occur
      %  param.tomo_collate.img_comb(3): guard time which specifies how
      %    many seconds at the end of img1 will not be used... this is
      %    important because the last samples of img1 will have low signal
      %    power and blurred because they will only have captured a portion
      %    of the chirp energy (typically this will be set to something
      %    close to the pulse duration for img1)
      %  param.tomo_collate.img_comb(4-6, 7-9, etc.): same fields as above
      %    except between images 2 and 3, 3 and 4, etc.
      
      Surface = interp_finite(mdata{img}.Surface,0);
      % First row of img_bins indicates the start of the blend-region
      img_bins = round(interp1(New_Time, 1:length(New_Time), ...
        max(Surface+param.tomo_collate.img_comb((img-2)*3+1),param.tomo_collate.img_comb((img-2)*3+2)), 'linear','extrap'));
      
      % Determine guard at end of image 1 that will not be used
      guard_bins = 1 + round(param.tomo_collate.img_comb((img-2)*3+3)/dt);
      
      % Check to make sure requested time is inside window and just
      % force the combination bin to occur at the second to last bin
      %   img_bins outside the img1 time window will be NaN due to interp1
      %   img_bins inside the img1 time window may still be larger than
      %     the guard allows
      max_good_time = length(Time)*ones(1,Nx);
      invalid_rlines = find(isnan(img_bins) ...
        | img_bins > max_good_time-guard_bins);
      img_bins(invalid_rlines) = max_good_time(invalid_rlines)-guard_bins;
      
      % Second row of img_bins indicates the end of the blend-region
      img_bins(2,:) = img_bins(1,:) + 10;
      
      difference = 10^(-0/10);
      
      % Combine images
      New_Data = zeros(size(mdata{img}.Data),'single');
      New_Topography = zeros(size(mdata{img}.Topography.img),'single');
      for rline = 1:Nx
        trans_bins = img_bins(1,rline)+1:img_bins(2,rline);
        weights = 0.5+0.5*cos(pi*linspace(0,1,length(trans_bins)).');
        if trans_bins(end) <= size(mdata{img}.Data,1)
          New_Data(:,rline) = [Data(1:img_bins(1,rline),rline); ...
            weights.*Data(trans_bins,rline) ...
            + difference*(1-weights).*mdata{img}.Data(trans_bins,rline); ...
            difference*mdata{img}.Data(img_bins(2,rline)+1:end,rline)];
          New_Topography(:,:,rline) = [Topography.img(1:img_bins(1,rline),:,rline); ...
            bsxfun(@times,weights,Topography.img(trans_bins,:,rline)) ...
            + bsxfun(@times,difference*(1-weights),mdata{img}.Topography.img(trans_bins,:,rline)); ...
            difference*mdata{img}.Topography.img(img_bins(2,rline)+1:end,:,rline)];
        else
          New_Data(:,rline) = Data(1:size(New_Data,1),rline);
          New_Topography(:,:,rline) = Topography.img(1:size(New_Data,1),:,rline);
        end
      end
      Time = New_Time;
      Data = New_Data;
      Topography.img = New_Topography;
    end
  end
  
  % Append new fused image to the combined file
  save(combined_fn,'-append','Time','Data','Topography');
  
  mdata = mdata{1};
  mdata.Time = Time;
  mdata.Data = Data;
  mdata.Topography.img = Topography.img;
end

end
