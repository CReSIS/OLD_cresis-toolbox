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

[~,radar_type,~] = ct_output_dir(param.radar_name);

% in_dir: Directory where 3D image files are at
in_dir = ct_filename_out(param,param.tomo_collate.in_path);

if ~isfield(param.tomo_collate, 'img_comb_trim') || isempty(param.tomo_collate.img_comb_trim)
  first_img = param.tomo_collate.imgs{1}(1);
  last_img = param.tomo_collate.imgs{end}(1);
  if strcmpi(radar_type,'deramp')
    param.tomo_collate.img_comb_trim = [0 0 0 inf];
  else
    % Set relative trim to be 50% support level or higher for pulse compression
    % Set absolute trim to be >= 0 time
    if iscell(param.array.imgs{first_img})
      % Multilook format (param_mode is 'array')
      wf_adc_list = param.array.imgs{first_img}{1};
    else
      wf_adc_list = param.array.imgs{first_img};
    end
    wf_first = wf_adc_list(1,1);
    if iscell(param.array.imgs{last_img})
      % Multilook format (param_mode is 'array')
      wf_adc_list = param.array.imgs{last_img}{1};
    else
      wf_adc_list = param.array.imgs{last_img};
    end
    wf_last = wf_adc_list(1,1);
    param.tomo_collate.img_comb_trim = [param.radar.wfs(wf_first).Tpd/2 -param.radar.wfs(wf_last).Tpd/2 0 inf];
  end
end

%% Horizontal fuse
% =========================================================================
mdata = {};
for v_img = 1:length(param.tomo_collate.imgs)
  Nimg = length(param.tomo_collate.imgs{v_img});
  fns = {};
  hdata = {};
  Topography = [];
  for h_img = 1:length(param.tomo_collate.imgs{v_img})
    img = param.tomo_collate.imgs{v_img}(h_img);
    wf = param.array.imgs{v_img}(1,1);
    fns{h_img} = fullfile(in_dir,sprintf('Data_img_%02.0f_%s_%03.0f.mat',img, ...
      param.day_seg,param.load.frm));
    hdata{h_img} = load(fns{h_img});
  end
  Nx = size(hdata{1}.Topography.img,3);
  
  if v_img == 1
    % Copy first image file to combined "fused" image filename (this is a simple way
    % to copy all the parameters and support variables over). We will update
    % the 3D image and time fields in the file at the end of fusing.
    combined_fn = fullfile(in_dir,sprintf('Data_%s_%03.0f.mat',param.day_seg,param.load.frm));
    fprintf('Creating %s (%s)\n', combined_fn, datestr(now));
    copyfile(fns{1}, combined_fn);
  end
  
  % Interpolate images onto a common propagation time axis
  % dt: Fast time sample period 
  dt = hdata{1}.Time(2)-hdata{1}.Time(1);
  % Time: common time axis between all images
  start_time = -inf;
  stop_time = inf;
  for img_idx = 1:Nimg
    if hdata{img_idx}.Time(1) > start_time
      start_time = hdata{img_idx}.Time(1);
    end
    if hdata{img_idx}.Time(end) < stop_time
      stop_time = hdata{img_idx}.Time(end);
    end
  end
  Time = start_time:dt:stop_time;

  % Reinterpolate images onto common time axis
  Data = interp1(hdata{1}.Time,hdata{1}.Data,Time);
  for img_idx = 1:Nimg
    hdata{img_idx}.Topography.img = interp1(hdata{img_idx}.Time,hdata{img_idx}.Topography.img,Time);
  end
  
  % Nx: number of range lines
  Nx = size(hdata{1}.Topography.img,3);
  % Nsv: number of steering vector/directions of arrival
  Nsv = size(hdata{1}.Topography.img,2);

  % Create fuse weights
  fuse_weights = [];
  Topography.img = zeros(size(hdata{1}.Topography.img));
  if Nimg > 1
    for img_idx = 1:Nimg
      fuse_weights(img_idx,:) = exp( -((0:Nsv-1)-Nsv*((Nimg-img_idx)*2+1)/(2*Nimg)).^2 / (2*(Nsv/Nimg).^1.5) );
      Topography.img = Topography.img + bsxfun(@times,fuse_weights(img_idx,:), ...
        hdata{img_idx}.Topography.img);
    end
    Topography.img = bsxfun(@times,1./sum(fuse_weights),Topography.img);
  else
    Topography.img = hdata{img_idx}.Topography.img;
  end
  
  mdata{v_img} = hdata{1};
  clear hdata;
  mdata{v_img}.Time = Time;
  mdata{v_img}.Data = Data;
  mdata{v_img}.Topography.img = Topography.img;
end

%% Vertical fuse
% =========================================================================
Topography = [];
for v_img = 1:length(param.tomo_collate.imgs)
  if v_img == 1
    Time = mdata{v_img}.Time;
    Data = mdata{v_img}.Data;
    Topography.img = mdata{v_img}.Topography.img;
    first_idx = find(Time >= Time(1)+param.tomo_collate.img_comb_trim(1) ...
      & Time >= param.tomo_collate.img_comb_trim(3),1,'first');
    if ~isempty(first_idx)
      Time = Time(first_idx:end);
      Data = Data(first_idx:end,:);
      Topography.img = Topography.img(first_idx:end,:,:);
    else
      error('Zero range bin length images not supported.');
    end
    if v_img == length(param.tomo_collate.imgs)
      last_idx = find(Time <= Time(end)+param.tomo_collate.img_comb_trim(2) ...
        & Time <= param.tomo_collate.img_comb_trim(4),1,'last');
      if ~isempty(last_idx)
        Time = Time(1:last_idx);
        Data = Data(1:last_idx,:);
        Topography.img = Topography.img(1:last_idx,:,:);
      else
        error('Zero range bin length images not supported.');
      end
    end
    
  else
    % Combine images
    % Data,Time => already loaded data
    % mdata{v_img}.Data, mdata{v_img}.Time => new data to mdata{v_img}
    % New_Time, New_Data => Combined result
    
    if v_img == length(param.tomo_collate.imgs)
      last_idx = find(mdata{v_img}.Time <= mdata{v_img}.Time(end)+param.tomo_collate.img_comb_trim(2) ...
        & mdata{v_img}.Time <= param.tomo_collate.img_comb_trim(4),1,'last');
      if ~isempty(last_idx)
        mdata{v_img}.Time = mdata{v_img}.Time(1:last_idx);
        mdata{v_img}.Data = mdata{v_img}.Data(1:last_idx,:);
        mdata{v_img}.Topography.img = mdata{v_img}.Topography.img(1:last_idx,:,:);
      else
        error('Zero range bin length images not supported.');
      end
    end
    
    % Interpolate image N onto already loaded data (assumption is that image
    % N-1 always comes before image N)
    New_Time = (Time(1) : dt : mdata{v_img}.Time(end)).';
    mdata{v_img}.Data = interp1(mdata{v_img}.Time,mdata{v_img}.Data,New_Time,'linear',0);
    mdata{v_img}.Topography.img = interp1(mdata{v_img}.Time,mdata{v_img}.Topography.img,New_Time,'linear',0);
    
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
    
    Surface = interp_finite(mdata{v_img}.Surface,0);
    % First row of img_bins indicates the start of the blend-region
    img_bins = round(interp1(New_Time, 1:length(New_Time), ...
      max(Surface+param.tomo_collate.img_comb((v_img-2)*3+1),param.tomo_collate.img_comb((v_img-2)*3+2)), 'linear','extrap'));
    
    % Determine guard at end of image 1 that will not be used
    guard_bins = 1 + round(param.tomo_collate.img_comb((v_img-2)*3+3)/dt);
    
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
    New_Data = zeros(size(mdata{v_img}.Data),'single');
    New_Topography = zeros(size(mdata{v_img}.Topography.img),'single');
    for rline = 1:Nx
      trans_bins = img_bins(1,rline)+1:img_bins(2,rline);
      weights = 0.5+0.5*cos(pi*linspace(0,1,length(trans_bins)).');
      if trans_bins(end) <= size(mdata{v_img}.Data,1)
        New_Data(:,rline) = [Data(1:img_bins(1,rline),rline); ...
          weights.*Data(trans_bins,rline) ...
          + difference*(1-weights).*mdata{v_img}.Data(trans_bins,rline); ...
          difference*mdata{v_img}.Data(img_bins(2,rline)+1:end,rline)];
        New_Topography(:,:,rline) = [Topography.img(1:img_bins(1,rline),:,rline); ...
          bsxfun(@times,weights,Topography.img(trans_bins,:,rline)) ...
          + bsxfun(@times,difference*(1-weights),mdata{v_img}.Topography.img(trans_bins,:,rline)); ...
          difference*mdata{v_img}.Topography.img(img_bins(2,rline)+1:end,:,rline)];
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
if param.ct_file_lock
  file_version = '1L';
else
  file_version = '1';
end
save(combined_fn,'-append','Time','Data','Topography','file_version');

mdata = mdata{1};
mdata.Time = Time;
mdata.Data = Data;
mdata.Topography.img = Topography.img;
