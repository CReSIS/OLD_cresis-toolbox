% script wise_data_nc_converter.m
%
% convert WISE data into CReSIS data format to be loaded by
% picker
%
% 1. Run wise_data_converter or wise_data_nc_converter
%   - Creates CSARP_standard and CSARP_layerData
% 2. Run fix_old_data_to_cr1
%
% 20100321_12_001

wise_data_path = '/cresis/snfs1/data/WISE/'; 
out_path = '/cresis/scratch2/mdce/wise/';
fns = get_filenames(wise_data_path,'','','','recursive');
splitDir = regexp(fns,'/','split');
times = {};
years = [];
months = [];
days = [];
for ii = 1:length(splitDir)   % get time stamps of each file
  fn = splitDir{ii}(7);
  str_locater = cell2mat(strfind(fn,'T'));
  times{ii} = fn{1}(str_locater - 8 : str_locater+6);
end
[times,fn_idxs] = sort(times); % sort the files according to time stamps

c = 3e8; % speed of light in air
date_str0 = times{1}(1:8);
seg0 = 0;
for ii = 1:length(times) 
%   fns{ii} = '/cresis/snfs1/data/WISE/2010_Greenland_TO/chirp_2_5_0.8V_1200s_20100319T163441.nc';
  fn = fns{fn_idxs(ii)};
  fprintf('%s (%s)\n', fn, datestr(now,'HH:MM:SS'));
  try
    Data = flipud((ncread(fn,'real_part')+1j*ncread(fn,'imag_part')).');
  catch err
    warning('FAILED TO LOAD');
    continue;
  end
  
  %Create/Update Directory 
  date_str = times{ii}(1:8);
  if ~strcmp(date_str,date_str0)  % reset the segment number for a new date
    seg0 = 0;
    date_str0 = date_str;
  end
  seg0 = seg0 + 1;
  seg = [date_str,sprintf('_%02d',seg0)];
  year = times{ii}(1:4);
  data_dir = fullfile(out_path,[year,'_Greenland_TO_wise'],'CSARP_standard',seg);    % use wise to differ from mcrds of the same year and location
  layer_dir = fullfile(out_path,[year,'_Greenland_TO_wise'],'CSARP_layerData',seg); 
  if ~exist(data_dir,'dir')
    mkdir(data_dir)
  end
  if ~exist(layer_dir,'dir')
    mkdir(layer_dir)
  end
  data_fn = fullfile(data_dir,['Data_',seg,'_001.mat']);
  layer_fn = fullfile(layer_dir,['Data_',seg,'_001.mat']);
  fprintf('  %s\n', data_fn);
 
  % read y axis
  y_axis = ncread(fn,'yaxis'); % yaxis is meter in air
  y_axis1 = double(flipud(y_axis));


  % read geo location coordinates
  Latitude = double(ncread(fn,'latitude'));
  Longitude = double(ncread(fn,'longitude'));
  Elevation = double(ncread(fn,'aircraft_altitude'));

  % Convert gps sod to seconds since Jan 1 1970
  year = str2num(times{ii}(1:4));
  month = str2num(times{ii}(5:6));
  day = str2num(times{ii}(7:8));
  hour = str2num(times{ii}(10:11));
  minute = str2num(times{ii}(12:13));
  second = str2num(times{ii}(14:15));;
  GPS_time0 = datenum_to_epoch(datenum(year,month,day,hour,minute,second));
  GPS_time = linspace(GPS_time0,GPS_time0+1200,size(Data,2)); % assuming 2 min for each file
 
  % Load ice surface and bottom picks
  Bottom = double(ncread(fn,'bedrock'));
  Surface = double(ncread(fn,'surface'));

  
%   figure(1); clf;
%   imagesc([],y_axis1,lp(Data))
%   set(gca,'YDir','normal')
%   hold on
%   plot(Surface,'k')
%   plot(Bottom,'k')
%   hold off;


  %% Multilook and decimate in along-track and make sure vectors are correct shape also
  Data = fir_dec(Data,10);
  GPS_time = fir_dec(reshape(GPS_time,[1 numel(GPS_time)]),10);
  Latitude = fir_dec(reshape(Latitude,[1 numel(Latitude)]),10);
  Longitude = fir_dec(reshape(Longitude,[1 numel(Longitude)]),10);
  Elevation = fir_dec(reshape(Elevation,[1 numel(Elevation)]),10);
  Surface = fir_dec(reshape(Surface,[1 numel(Surface)]),10);
  Bottom = fir_dec(reshape(Bottom,[1 numel(Bottom)]),10);
  
  %% Convert data back to twtt from range
  First_Time = (Elevation - y_axis1(1))*2/c;
  Last_Time = (Elevation - y_axis1(end))*2/c;
  t0 = max(0,max(First_Time));
  dt = abs(y_axis1(2)-y_axis1(1)) * 2/c;
  Nt = floor((min(Last_Time)-t0)/dt);
  Time = t0 + dt*(0:Nt-1);
  er_ice = 3.15;
  for rline = 1:size(Data,2)
    % Find in air indices
    air_idxs = y_axis1 >= Surface(rline);
    rline_time = (Elevation(rline) - y_axis1(air_idxs))*2/c;
    ice_idxs = y_axis1 < Surface(rline);
    rline_time = cat(1,rline_time,(Elevation(rline)-Surface(rline))*2/c - (y_axis1(ice_idxs) - Surface(rline))*2/c);
    Data(1:Nt,rline) = interp1(rline_time,Data(:,rline),Time);
  end
  Data = Data(1:Nt,:);
  %convert from [m]to [s]
  Surface1 = (Elevation - Surface)*2/c;
  Bottom =  Surface1 + (Surface-Bottom)*2/(c/sqrt(er_ice));
  Surface = Surface1;

%   figure(2); clf;
%   imagesc([],Time,lp(Data))
%   hold on
%   plot(Surface,'k')
%   plot(Bottom,'k')
%   hold off;
  
  Time = reshape(Time,[numel(Time) 1]); 
  Depth = Time * c/2;
  
  % Create a layer struct
  for layer_idx = 1:2
    % Manually picked points
    %  inf/nan: no pick
    %  finite: propagation time to target
    Nx=length(Data);
    
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