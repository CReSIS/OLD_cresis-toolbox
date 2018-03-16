function combine_wf_chan_ollie2(steady_param_file_name)
% combine_wf_chan_ollie2(steady_param_file_name)
%
% This script combines the receive channels and outputs the result
% for each waveform. It also combines the waveforms. It takes in
% f-k files one directory at a time and:
%  6. Combines the waveforms
%
% The assumption is that the directories in the input_path are named
% using the following convention:
%   PROC-TYPE-STRING_data_#{_SUBAPERTURE-STRING}
% where
%   PROC-TYPE-STRING can be 'fk','tdbp', or 'pc' for f-k migrated,time domain
%   back projected,and pulse compressed respectively ('fk' and tdbp supported)
%   _data_ is always present
%   #, \d+: one or more numbers
%   _SUBAPERTURE-STRING, {_[mp]\d\.\d}: optional subaperture string
% Examples:
%   fk_data_01_01: f-k migrated, frame 1, subaperture 1
%   fk_data_04_02: f-k migrated, frame 4, subaperture 2
%   fk_data_01_03: f-k migrated, frame 1, subaperture 3
%   pc_data_01: pulse compressed only, frame 1
%
% param = struct with processing parameters loaded from file
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_combine_wf_chan.m, combine_wf_chan.m,
%   combine_wf_chan_task.m

% =====================================================================
% General Setup
% =====================================================================

dbstack_info = dbstack;
  
load(steady_param_file_name,'steady_param');
param=steady_param;

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', dbstack_info(1).name, param.day_seg, datestr(now,'HH:MM:SS'));
fprintf('=====================================================================\n');

% =====================================================================
% Setup processing
% =====================================================================

% Get WGS84 ellipsoid parameters
physical_constants;

if ~isfield(param,'debug_level')
  param.debug_level = 1;
end

if ~isfield(param.sched,'rerun_only') || isempty(param.sched.rerun_only)
  param.sched.rerun_only = false;
end

% Handles multilooking syntax:
%  {{[1 1],[1 2],[1 3],[1 4],[1 5]},{[2 1],[2 2],[2 3],[2 4],[2 5]}}
%  If the image is a cell array it describes multilooking across apertures
if ~iscell(param.combine.imgs{1})
  % No special multilooking, reformat old syntax to new multilooking syntax
  for img = 1:length(param.combine.imgs)
    param.combine.imgs{img} = {param.combine.imgs{img}};
  end
end

for img = 1:length(param.combine.imgs)
  for ml_idx = 1:length(param.combine.imgs{img})
    % Imaginary image indices is for IQ combining during raw data load
    % which we do not need here.
    param.combine.imgs{img}{ml_idx} = abs(param.combine.imgs{img}{ml_idx});
  end
end

img_list = param.combine.imgs;
  
in_path = ct_filename_out(param, ...
  param.combine.in_path, 'CSARP_out');

array_path = ct_filename_out(param, ...
  param.combine.array_path, 'CSARP_out');

out_path = ct_filename_out(param, ...
  param.combine.out_path, sprintf('CSARP_%s', ...
  param.combine.method));

% Create the output directory
if ~exist(out_path,'dir')
  mkdir(out_path);
end

% Load frames file
load(ct_filename_support(param,'','frames'));

param.surf.manual = 0; % Turn manual pick off

%% Loop through all the frames
% =====================================================================
for frm_idx = 1:length(param.cmd.frms);
  frm = param.cmd.frms(frm_idx);
  
  if ct_proc_frame(frames.proc_mode(frm),param.csarp.frm_types)
    fprintf('%s combine %s_%03i (%i of %i) %s\n', param.radar_name, param.day_seg, frm, frm_idx, length(param.cmd.frms), datestr(now,'HH:MM:SS'));
  else
    fprintf('Skipping frame %s_%03i (no process frame)\n', param.day_seg, frm);
    continue;
  end
  
  %% Output directory
  param.combine.out_path = fullfile(array_path,sprintf('array_%03d', frm));
  
  %% Loop through all the images
  for img = 1:length(param.combine.imgs)
    
    %% Loop through all the chunks and combine
    Latitude = [];
    Longitude = [];
    Elevation = [];
    Roll = [];
    Pitch = [];
    Heading = [];
    GPS_time = [];
    Surface = [];
    Bottom = [];
    Data = [];
    Topography = [];
    chunk_fns = get_filenames(param.combine.out_path,'chk','',sprintf('img_%02d.mat',img));
    for chunk_idxs = 1:length(chunk_fns)
      tmp = load(chunk_fns{chunk_idxs});
      Time = tmp.Time;
      Latitude = [Latitude double(tmp.Latitude)];
      Longitude = [Longitude double(tmp.Longitude)];
      Elevation = [Elevation double(tmp.Elevation)];
      Roll = [Roll double(tmp.Roll)];
      Pitch = [Pitch double(tmp.Pitch)];
      Heading = [Heading double(tmp.Heading)];
      GPS_time = [GPS_time tmp.GPS_time];
      Surface = [Surface double(tmp.Surface)];
      Bottom = [Bottom double(tmp.Bottom)];
      Data = [Data tmp.Data];
      param_records = tmp.param_records;
      param_csarp = tmp.param_csarp;
      if chunk_idxs == 1
      param_combine = tmp.param_combine;
        param_combine.array_param.fcs{1}{1}.x = tmp.param_combine.array_param.fcs{1}{1}.x(:,tmp.param_combine.array_param.lines);
        param_combine.array_param.fcs{1}{1}.y = tmp.param_combine.array_param.fcs{1}{1}.y(:,tmp.param_combine.array_param.lines);
        param_combine.array_param.fcs{1}{1}.z = tmp.param_combine.array_param.fcs{1}{1}.z(:,tmp.param_combine.array_param.lines);
        param_combine.array_param.fcs{1}{1}.origin = tmp.param_combine.array_param.fcs{1}{1}.origin(:,tmp.param_combine.array_param.lines);
      else
        % Concatenate the fcs field
        param_combine.array_param.fcs{1}{1}.x = [param_combine.array_param.fcs{1}{1}.x tmp.param_combine.array_param.fcs{1}{1}.x(:,tmp.param_combine.array_param.lines)];
        param_combine.array_param.fcs{1}{1}.y = [param_combine.array_param.fcs{1}{1}.y tmp.param_combine.array_param.fcs{1}{1}.y(:,tmp.param_combine.array_param.lines)];
        param_combine.array_param.fcs{1}{1}.z = [param_combine.array_param.fcs{1}{1}.z tmp.param_combine.array_param.fcs{1}{1}.z(:,tmp.param_combine.array_param.lines)];
        param_combine.array_param.fcs{1}{1}.origin = [param_combine.array_param.fcs{1}{1}.origin tmp.param_combine.array_param.fcs{1}{1}.origin(:,tmp.param_combine.array_param.lines)];
      end
      if isfield(tmp,'Topography')
%         3D-surface is present so concatenate it too
%         Topography = cat(3,Topography,tmp.Topography);
%         Concatenate all the fields under struct Topography: valR, bins, val, freq
%         and img.
        fields = fieldnames(tmp.Topography);      
        if chunk_idxs == 1
          for field_idx = 1:length(fields)
            Topography.(fields{field_idx}) = tmp.Topography.(fields{field_idx});
          end       
        else        
          for field_idx = 1:length(fields)
            max_dim = length(size(tmp.Topography.(fields{field_idx})));
            Topography.(fields{field_idx}) = cat(max_dim,Topography.(fields{field_idx}),tmp.Topography.(fields{field_idx}));
          end       
        end
        
      end
    end
    
    % =====================================================================
    % Save output
    if length(param.combine.imgs) == 1
      out_fn = fullfile(out_path, sprintf('Data_%s_%03d.mat', ...
        param.day_seg, frm));
    else
      out_fn = fullfile(out_path, sprintf('Data_img_%02d_%s_%03d.mat', ...
        img, param.day_seg, frm));
    end
    fprintf('  Writing output to %s\n', out_fn);
    if isempty(Topography)
      % Do not save 3D surface
      save('-v7.3',out_fn,'Time','Latitude','Longitude', ...
        'Elevation','GPS_time','Data','Surface','Bottom', ...
        'param_combine','param_records','param_csarp', ...
        'Roll', 'Pitch', 'Heading');
    else
      % Save 3D surface
      save('-v7.3',out_fn,'Topography','Time','Latitude', ...
        'Longitude','Elevation','GPS_time','Data','Surface','Bottom', ...
        'param_combine','param_records','param_csarp', ...
        'Roll', 'Pitch', 'Heading');
    end
  end
  
  if isempty(param.combine.img_comb)
    % No image combining is required
    continue;
  end
  
  if length(param.combine.img_comb) ~= 3*(length(param.combine.imgs)-1)
    warning('param.combine.img_comb not the right length. There should be 3 entries for each image combination interface ([Tpd second image for surface saturation, -inf for second image blank, Tpd first image to avoid roll off] is typical). Set correctly here and update param spreadsheet before dbcont.');
    keyboard
  end
  
  %% Load each image and then combine with previous image (also trim time<0 values)
  for img = 1:length(param.combine.imgs)
    
    if length(param.combine.imgs) == 1
      out_fn = fullfile(out_path, sprintf('Data_%s_%03d.mat', ...
        param.day_seg, frm));
    else
      out_fn = fullfile(out_path, sprintf('Data_img_%02d_%s_%03d.mat', ...
        img, param.day_seg, frm));
    end
    if img == 1
      load(out_fn);
      first_idx = find(Time <= 0,1,'last');
      if ~isempty(first_idx)
        Time = Time(first_idx:end);
        Data = Data(first_idx:end,:);
      end
    else
      append = load(out_fn,'Time','Data');
      %% Combine images
      % Data,Time => already loaded data
      % append.Data, append.Time => new data to append
      % New_Time, New_Data => Combined result
      
      % Interpolate image N onto already loaded data (assumption is that image
      % N-1 always comes before image N)
      dt = Time(2)-Time(1);
      New_Time = (Time(1) : dt : append.Time(end)).';
      append.Data = interp1(append.Time,append.Data,New_Time,'linear',0);
      
      % Surface tracking image combine
      %  param.combine.img_comb(1): time after surface return where
      %    combine will happen
      %  param.combine.img_comb(2): minimum time that combine will occur
      %  param.combine.img_comb(3): guard time which specifies how
      %    many seconds at the end of img1 will not be used... this is
      %    important because the last samples of img1 will have low signal
      %    power and blurred because they will only have captured a portion
      %    of the chirp energy (typically this will be set to something
      %    close to the pulse duration for img1)
      %  param.combine.img_comb(4-6, 7-9, etc.): same fields as above
      %    except between images 2 and 3, 3 and 4, etc.
      
      Surface = interp_finite(Surface,0);
      % First row of img_bins indicates the start of the blend-region
      img_bins = round(interp1(New_Time, 1:length(New_Time), ...
        max(Surface+param.combine.img_comb((img-2)*3+1),param.combine.img_comb((img-2)*3+2)), 'linear','extrap'));
      
      % Determine guard at end of image 1 that will not be used
      guard_bins = 1 + round(param.combine.img_comb((img-2)*3+3)/dt);
      
      % Check to make sure requested time is inside window and just
      % force the combination bin to occur at the second to last bin
      %   img_bins outside the img1 time window will be NaN due to interp1
      %   img_bins inside the img1 time window may still be larger than
      %     the guard allows
      max_good_time = length(Time)*ones(1,size(Data,2));
      invalid_rlines = find(isnan(img_bins) ...
        | img_bins > max_good_time-guard_bins);
      img_bins(invalid_rlines) = max_good_time(invalid_rlines)-guard_bins;
      
      % Second row of img_bins indicates the end of the blend-region
      img_bins(2,:) = img_bins(1,:) + 1;
      
      difference = 10^(-0/10);
      
      % Combine images
      New_Data = zeros(size(append.Data),'single');
      for rline = 1:size(New_Data,2)
        trans_bins = img_bins(1,rline)+1:img_bins(2,rline);
        weights = 0.5+0.5*cos(pi*linspace(0,1,length(trans_bins)).');
        if trans_bins <= size(append.Data,1)
          New_Data(:,rline) = [Data(1:img_bins(1,rline),rline); ...
            weights.*Data(trans_bins,rline) ...
            + difference*(1-weights).*append.Data(trans_bins,rline); ...
            difference*append.Data(img_bins(2,rline)+1:end,rline)];
        else
          New_Data(:,rline) = Data(1:size(New_Data,1),rline);
        end        
      end
      Time = New_Time;
      Data = New_Data;
    end
  end
  
  %% Save output
  out_fn = fullfile(out_path, sprintf('Data_%s_%03d.mat', ...
    param.day_seg, frm));
  fprintf('  Writing output to %s\n', out_fn);
  save('-v7.3',out_fn,'Time','Latitude','Longitude', ...
    'Elevation','GPS_time','Data','Surface','Bottom', ...
    'param_combine','param_records','param_csarp', ...
    'Roll', 'Pitch', 'Heading');
end

return;