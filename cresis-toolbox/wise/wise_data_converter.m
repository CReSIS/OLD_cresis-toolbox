% convert WISE data into CReSIS data format to be loaded by picker
wise_data_path = '/cresis/data2/WISE/Helheim_2008/';
out_path = '/cresis/scratch2/mdce/wise/';
fns.bed = get_filenames(wise_data_path,'bed','','txt','recursive');
fns.gps = get_filenames(wise_data_path,'','datgps.dat','','recursive');
fns.data = get_filenames(wise_data_path,'','','dat','recursive');  % no echogram data in 20090409T090025
dt = 1/100e6;
c = 3e8;

seg1 = 0;
seg2 = 0;
seg3 = 0;
for ii = 1:length(fns.data)
  clear data gps layerData 
  % Load radar echogram
  fid = fopen(fns.data{ii},'r');
  Nx = fread(fid,1,'long');
  Ny = fread(fid,1,'long');
  xAxis = fread(fid,Nx,'float32');
  yAxis = fread(fid,Ny,'float32');
  Data = fread(fid,[Ny Nx],'float32');
  Time = [0:size(Data,1)-1]*dt;
  Depth = Time*c/sqrt(3.15)/2;
  fclose(fid);
  str_locater = strfind(fns.data{ii},'T');
  id_str = fns.data{ii}(str_locater(2) - 8 : str_locater(2)+6);
  date_str = id_str(1:8);
  if strcmp(date_str,'20080507')
    seg1 = seg1 + 1;
    seg = [date_str,sprintf('_%02d',seg1)];
    data_dir = fullfile(out_path,'2008_Greenland_TO_wise','CSARP_standard',seg);    % use wise to differ from mcrds of the same year and location
    layer_dir = fullfile(out_path,'2008_Greenland_TO_wise','CSARP_layerData',seg);
  elseif strcmp(date_str,'20090414')
    seg2 = seg2 + 1;
    seg = [date_str,sprintf('_%02d',seg2)];
    data_dir = fullfile(out_path,'2009_Greenland_TO_wise','CSARP_standard',seg);
    layer_dir = fullfile(out_path,'2009_Greenland_TO_wise','CSARP_layerData',seg);
  elseif strcmp(date_str,'20090409')
    seg3 = seg3 + 1;
    seg = [date_str,sprintf('_%02d',seg3)];
    data_dir = fullfile(out_path,'2009_Greenland_TO_wise','CSARP_standard',seg);
    layer_dir = fullfile(out_path,'2009_Greenland_TO_wise','CSARP_layerData',seg);
  end
  if ~exist(data_dir,'dir')
    mkdir(data_dir);
  end
  if ~exist(layer_dir,'dir')
    mkdir(layer_dir);
  end
  data_fn = fullfile(data_dir,['Data_',seg,'_001.mat']);
  layer_fn = fullfile(layer_dir,['Data_',seg,'_001.mat']);  

  % Load gps 
  for ij = 1:length(fns.gps)
    if isempty(strfind(fns.gps{ij},id_str))
      continue
    else
      fn_gps = fns.gps{ij};
      break
    end
  end
  gps = dlmread(fn_gps);
  Latitude = gps(:,2)';
  Longitude = gps(:,3)';
  Elevation = gps(:,4)';
  
  % Convert gps sod to seconds since Jan 1 1970
  year = str2num(date_str(1:4));
  month = str2num(date_str(5:6));
  day = str2num(date_str(7:8));
  GPS_time = gps(:,1)';
  if GPS_time(1)>86400
    GPS_time = GPS_time-86400;
    day = day +1;
  end
  if GPS_time(1) == 0
    hour = str2num(id_str(10:11));
    minute = str2num(id_str(12:13));
    second = str2num(id_str(14:15));;
    GPS_time0 = datenum_to_epoch(datenum(year,month,day,hour,minute,second));
    GPS_time = GPS_time + GPS_time0;
  else
    % remove same GPS_time and handle jumps
    [GPS_time,m,n] = unique(GPS_time); % gps time provided in first file
    Latitude = Latitude(m);
    Longitude = Longitude(m);
    Elevation = Elevation(m);
    dgps_time = diff(GPS_time);
    jump_idxs = find(dgps_time>1.5);
    if ~isempty(jump_idxs)
      gps_time_step = mean(dgps_time(1:jump_idxs(1)-1));
      for jump_idx = 1:length(jump_idxs)
        GPS_time(jump_idxs(jump_idx)+1:end) = GPS_time(jump_idxs(jump_idx)+1:end) - dgps_time(jump_idxs(jump_idx)) + gps_time_step;
      end
    end
  end
  GPS_time_itp = linspace(GPS_time(1),GPS_time(end),size(Data,2));
  Latitude = interp1(GPS_time,Latitude,GPS_time_itp);
  Longitude = interp1(GPS_time,Longitude,GPS_time_itp);
  Elevation = interp1(GPS_time,Elevation,GPS_time_itp);
  GPS_time = GPS_time_itp;
  
  % Load ice surface and bottom picks
  for ik = 1:length(fns.bed)
    if isempty(strfind(fns.bed{ik},id_str))
      continue
    else
      fn_bed = fns.bed{ik};
      break
    end
  end
  bed_data = load(fn_bed);
  bed_data(bed_data == -9999) = NaN;
  Surface = size(Data,1)-bed_data(:,2)';
  Surface = (Surface-1)*dt;
  Bottom = size(Data,1)- bed_data(:,3)';
  Bottom = (Bottom-1)*dt;
  
  % Create a layer struct
  for layer_idx = 1:2
    % Manually picked points
    %  inf/nan: no pick
    %  finite: propagation time to target
    layerData{layer_idx}.value{1}.data ...
      = inf*ones(1,Nx);
    % Automatically generated points
    %  inf/nan: no pick
    %  finite: propagation time to target
    if layer_idx == 1 % && isfield(lyr,'Surface')
      layerData{layer_idx}.value{2}.data = Surface;
    elseif layer_idx == 2 % && isfield(lyr,'Bottom')
      layerData{layer_idx}.value{2}.data = Bottom;
    else
      layerData{layer_idx}.value{2}.data ...
        = inf*ones(1,Nx);
    end
    % Quality control level
    %  1: good
    %  2: moderate
    %  3: derived
    layerData{layer_idx}.quality ...
      = ones(1,Nx);
  end
  
  % save converted data
  save(data_fn,'Latitude', 'Longitude', 'Elevation', 'GPS_time', 'Surface', 'Bottom', 'Data', 'Time', 'Depth');
  save(layer_fn, 'GPS_time', 'Latitude', 'Longitude', 'Elevation', 'layerData');
end
