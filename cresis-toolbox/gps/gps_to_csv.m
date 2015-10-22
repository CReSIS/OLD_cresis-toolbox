function gps_to_csv()
%% gps_to_csv.m
%
% Take GPS flight data and:
% (1) converts the GPS file to a structure
% (2) converts the structure to a CSV file that can be plotted or
%     manipulated with MATLAB or GIS software.
% (3) plots the data on a GeoTIFF
% There are already separate functions that read the GPS files into a
% structure; this script outputs them to a convenient, usable format.
%
% This function supports the following types of GPS files:
%   param.input_type = 1   --> .gps (NMEA) files
%   param.input_type = 2   --> ppp.txt (Novatel) files
%   param.input_type = 3   --> .out (Applanix) files
%   param.input_type = 4   --> .traj (ATM) files
%   param.input_type = 5   --> ln100g.asc (Litton 100) files
%
% EXAMPLE
% -------------------------------------------------------------------------
% param.input_type = 1;
% param.region = 'Antarctica';
% param.fn = ('P:\metadata\2012_Greenland_P3\BD960_10Apr12_PPPK_P13Jun12.out');
% param.fn = ('/cresis/projects/metadata/2012_Greenland_P3/BD960_10Apr12_PPPK_P13Jun12.out');
% param.dec = 11;
% param.year = 2011;
% param.month = 27;
% param.day = 10;
% param.time_reference = 'utc';
% param.csv_merge = false;
% csv_out_path = 'Y:\sfoga\2012_Greenland_P3_GPS.csv';
% csv_out_path = '/cresis/scratch1/sfoga/2012_Greenland_P3_GPS.csv';
% gps_to_csv(csv_out_path,param);
% -------------------------------------------------------------------------
%
% See 'edit gps_to_csv.m' to see the full list of param options.
%
% Date: 6 August 2012
% Author: Steven Foga
%
% See also read_gps_applanix.m, read_gps_atm.m, read_gps_litton.m,
% read_gps_novatel.m, read_gps_traj.m

% Specify GPS files to plot

% fns = get_filenames('/mnt/products/csarp_support/gps/2013_Antarctica_Basler/','gps','','.mat');
fns = get_filenames('/mnt/products/csarp_support/gps/2013_Antarctica_Ground/','gps','','.mat');

% Specify output directory
% out_fn_dir = '/home/cresis1/gps_csv/';
out_fn_dir = '/home/cresis1/gps_csv_accum/';

% Output along track spacing
spacing = 1000;

for fn_idx = 1:length(fns)
  fn = fns{fn_idx};
  [~,fn_name] = fileparts(fn);
  out_fn = fullfile(out_fn_dir,[fn_name '.csv']);
  fprintf('Converting %s\n to %s\n', fn, out_fn);
  
  gps = load(fn);
  along_track = geodetic_to_along_track(gps.lat,gps.lon,gps.elev);
  decim_idxs = get_equal_alongtrack_spacing_idxs(along_track,spacing);

  if ~exist(out_fn_dir,'dir')
    mkdir(out_fn_dir);
  end
  
  [fid,msg] = fopen(out_fn,'w');
  if fid < 0
    error('Failed to open %s: %s\n', out_fn, msg);
  end
  
  fprintf(fid, 'LAT,LON\n');
  fprintf(fid, '%f,%f\n',[gps.lat(decim_idxs); gps.lon(decim_idxs)]);
  
  fclose(fid);
  
end

return;






















return;


%% Example User Input

% Select type of input files that need to be converted to CSV files
%   param.input_type = 1   --> .gps (NMEA) files
%   param.input_type = 2   --> ppp.txt (Novatel) files
%   param.input_type = 3   --> .out (Applanix) files
%   param.input_type = 4   --> .traj (ATM) files
%   param.input_type = 5   --> ln100g.asc (Litton 100) files
%
%   param.input_type = 3;

% Specify region for GeoTIFF plot.
%   param.region = 'Greenland';
%   param.region = 'Antarctica';
%   param.region = 'Canada';

%   param.region = 'Greenland';

% Find all the input files from a directory
%  Example for one file:
%   param.fn = '/cresis/projects/metadata/2009_Antarctica_TO/rover_ppp_antarctica_12122009_LC.txt';
%
%  Example for multiple files:
%   param.fn = get_filenames('P:\metadata\2010_Greenland_P3_INS\','ln100g','','.asc');
%   param.fn = get_filenames('P:\metadata\2011_Antarctica_DC8\','','PPPK','.out');
%  For more help on finding multiple files, see 'get_filenames' Help Doc.
%   param.fn = get_filenames('/cresis/projects/metadata/2012_Greenland_P3/','','','PPPK.out');
%   param.fn = ('/cresis/projects/metadata/2012_Greenland_P3/BD960_10Apr12_PPPK_P13Jun12.out');

% Output file (use .csv extension)
%   csv_out_path = '/cresis/scratch2/sfoga/2012_Greenland_P3_GPS.csv';
%   csv_out_path = 'Z:\sfoga\test_gps_files\2010_Greenland_DC8_out.csv';

% Decimation - to keep every nth point
% Recommended: set .gps to 25 (Default), Novatel and Applanix (.out) to
% 1000 (Default).
% No recommended decimation on Littion or ATM; Default is 1000.
%   param.dec = 1000;

% Year/Month/Day of first day of data -- this helps to determine in what
% GPS week that data starts.
%   param.year = 2012;
%   param.month = 03;
%   param.day = 14;

% Time reference: should be 'utc' (Default), but try 'gps' if that fails.
%   param.time_reference = 'utc';

% If true, there will be only one CSV output file. If false, there will
% be one CSV file for each GPS file that is converted.
%   param.csv_merge = false ;


%% Automated section
fprintf('Start ... %s\n',datestr(now,'HH:MM:SS'));
clf; global gRadar;

% clear gps lat lon geotiff proj RGB R;

if strcmpi(param.region,'Greenland') | strcmpi(param.region,'Canada')
  geotiff = ct_filename_gis(gRadar,'greenland/Landsat-7/Greenland_natural.tif');
  proj = geotiffinfo(geotiff);
elseif strcmpi(param.region,'Antarctica')
  geotiff = ct_filename_gis(gRadar,'antarctica/Landsat-7/Antarctica_LIMA_480m.tif');
  proj = geotiffinfo(geotiff);
else
  fprintf('Region not specified or incorrectly spelled. GeoTIFF cannot be displayed.\n\n');
end
% Geotiff
figure(1);
tic; fprintf('Loading %s GeoTIFF...\n',param.region);
mapshow(geotiff); hold on; toc;

switch param.input_type
  case 1
    % if param.input_type == 1
    %% NMEA (.gps) to CSV
    
    % PARAMETERS
    
    header_lat = 'LAT';
    header_lon = 'LON';
    idx_size = size(param.fn);
    
    param.combine = 0;
    
    
    if isempty(param.time_reference)
      param.time_reference = 'utc';
    else
    end
    
    if isempty(param.dec)
      param.dec = 25;
    else
    end
    
    out_folder = fileparts(csv_out_path);
    if ~ispc
      out_folder = strcat(out_folder,'/');
    else
      out_folder = strcat(out_folder,'\');
    end
    
    % Build structure and write file(s)
    if param.csv_merge == false % One CSV file for every input
      fprintf('Processing %d file(s)...\n',idx_size(1));
      for fn_idx = 1:idx_size(1)
        if idx_size(1) > 1 % If more than one param.fn entry exists
          tic; gps = read_gps_nmea(param.fn{fn_idx},param); toc;
        else % If only one param.fn entry exists
          tic; gps = read_gps_nmea(param.fn,param); toc;
        end
        
        lat = gps.lat(1:param.dec:end);
        lon = gps.lon(1:param.dec:end);
        
        % Write out file
        if idx_size(1) > 1 % If more than one param.fn entry exists
          [a filename b] = fileparts(param.fn{fn_idx});
        else % If only one param.fn entry exists
          [a filename b] = fileparts(param.fn);
        end
        
        if ~ispc
          out_path = strcat('/',out_folder,filename,'.csv');
        else
          out_path = strcat(out_folder,filename,'.csv');
        end
        
        % Open and output CSV file
        output_file = strcat(filename,'.csv');
        fprintf('Writing out file %s and plotting map...\n',output_file);
        fid = fopen(out_path,'w');
        
        fprintf(fid,'%s,%s\n',header_lat,header_lon);
        for data_idx = 1:length(lat)
          fprintf(fid,'%f,%f\n',lat(data_idx),lon(data_idx));
          [x,y] = projfwd(proj,lat(data_idx),lon(data_idx));
          plot(x,y,'b-'); hold on;
        end
        fclose(fid);
      end
    else % One CSV file for any number of inputs
      fprintf('Processing %d file(s)...\n',idx_size(1));
      lat = {};
      lon = {};
      gps = {};
      for fn_idx = 1:idx_size(1)
        fprintf('File %d of %d...\n',fn_idx,idx_size(1));
        if idx_size(1) > 1 % If more than one param.fn entry exists
          tic; gps = read_gps_nmea(param.fn{fn_idx},param); toc;
          lat{fn_idx} = gps.lat(1:param.dec:end);
          lon{fn_idx} = gps.lon(1:param.dec:end);
        else % If only one param.fn entry exists
          tic; gps = read_gps_nmea(param.fn,param); toc;
          lat = gps.lat(1:param.dec:end);
          lon = gps.lon(1:param.dec:end);
        end
        
      end
      [a b c] = fileparts(csv_out_path);
      output_file = strcat(b,c);
      fprintf('Writing out file %s and plotting map...\n',output_file);
      fid = fopen(csv_out_path,'w');
      fprintf(fid,'%s,%s\n',header_lat,header_lon);
      for data_idx = 1:idx_size(1)
        for coord_idx = 1:length(lat{data_idx})
          fprintf(fid,'%f,%f\n',lat{data_idx}(coord_idx),lon{data_idx}(coord_idx));
          [x,y] = projfwd(proj,lat{data_idx}(coord_idx),lon{data_idx}(coord_idx));
          plot(x,y,'b-'); hold on;
        end
      end
      fclose(fid);
      
    end
    
  case 2
    %   elseif param.input_type == 2
    %% Novatel to CSV
    header_lat = 'LAT';
    header_lon = 'LON';
    idx_size = size(param.fn);
    
    % PARAMETERS
    if isempty(param.time_reference)
      param.time_reference = 'utc';
    else
    end
    
    if isempty(param.dec)
      param.dec = 1000;
    else
    end
    
    out_folder = fileparts(csv_out_path);
    if ~ispc
      out_folder = strcat(out_folder,'/');
    else
      out_folder = strcat(out_folder,'\');
    end
    
    % Build structure and write file(s)
    if param.csv_merge == false % One CSV file for every input
      fprintf('Processing %d file(s)...\n',idx_size(1));
      for fn_idx = 1:length(param.fn)
        clear gps a filename b out_path;
        
        fprintf('File %d of %d...\n',fn_idx,idx_size(1));
        idx_size = size(param.fn);
        if idx_size(1) > 1 % If more than one param.fn entry exists
          tic; gps = read_gps_novatel(param.fn{fn_idx},param); toc;
        else % If only one param.fn entry exists
          tic; gps = read_gps_novatel(param.fn,param); toc;
        end
        
        lat = gps.lat(1:param.dec:end);
        lon = gps.lon(1:param.dec:end);
        
        % Write out file
        if idx_size(1) > 1 % If more than one param.fn entry exists
          [a filename b] = fileparts(param.fn{fn_idx});
        else % If only one param.fn entry exists
          [a filename b] = fileparts(param.fn);
        end
        
        if ~ispc
          out_path = strcat('/',out_folder,filename,'.csv');
        else
          out_path = strcat(out_folder,filename,'.csv');
        end
        
        
        % Open and output CSV file
        output_file = strcat(filename,'.csv');
        fprintf('Writing out file %s and plotting map...\n',output_file);
        fid = fopen(csv_out_path,'w');
        fprintf(fid,'%s,%s\n',header_lat,header_lon);
        for data_idx = 1:length(lat)
          fprintf(fid,'%f,%f\n',lat(data_idx),lon(data_idx));
          [x,y] = projfwd(proj,lat{data_idx},lon{data_idx});
          plot(x,y,'b-'); hold on;
        end
        fclose(fid);
      end
    else % One CSV file for any number of inputs
      fprintf('Processing %d file(s)...\n',idx_size(1));
      lat = {};
      lon = {};
      gps = {};
      for fn_idx = 1:length(param.fn)
        fprintf('File %d of %d...\n',fn_idx,idx_size(1));
        if idx_size(1) > 1 % If more than one param.fn entry exists
          tic; gps = read_gps_novatel(param.fn{fn_idx},param); toc;
        else % If only one param.fn entry exists
          tic; gps = read_gps_novatel(param.fn,param); toc;
        end
        lat{fn_idx} = gps{fn_idx}.lat(1:param.dec:end);
        lon{fn_idx} = gps{fn_idx}.lon(1:param.dec:end);
      end
      [a b c] = fileparts(csv_out_path);
      output_file = strcat(b,c);
      fprintf('Writing out file %s and plotting map...\n',output_file);
      fid = fopen(csv_out_path,'w');
      fprintf(fid,'%s,%s\n',header_lat,header_lon);
      for data_idx = 1:length(lat)
        for coord_idx = 1:length(lat{data_idx})
          fprintf(fid,'%f,%f\n',lat{data_idx}(coord_idx),lon{data_idx}(coord_idx));
          [x,y] = projfwd(proj,lat{data_idx}(coord_idx),lon{data_idx}(coord_idx));
          plot(x,y,'b-'); hold on;
        end
      end
      fclose(fid);
    end
    
  case 3
    %   elseif param.input_type == 3
    %% Applanix (.out) to CSV
    header_lat = 'LAT';
    header_lon = 'LON';
    idx_size = size(param.fn);
    
    % PARAMETERS
    if isempty(param.time_reference)
      param.time_reference = 'utc';
    else
    end
    
    if isempty(param.dec)
      param.dec = 2;
    else
    end
    
    out_folder = fileparts(csv_out_path);
    if ~ispc
      out_folder = strcat(out_folder,'/');
    else
      out_folder = strcat(out_folder,'\');
    end
    
    if param.csv_merge == false % One CSV file for every input
      fprintf('Processing %d file(s)...\n',idx_size(1));
      for fn_idx = 1:idx_size(1)
        clear gps a filename b out_path;
        
        fprintf('File %d of %d...\n',fn_idx,idx_size(1));
        if idx_size(1) > 1 % If more than one param.fn entry exists
          tic; gps = read_gps_applanix(param.fn{fn_idx},param); toc;
        else % If only one param.fn entry exists
          tic; gps = read_gps_applanix(param.fn,param); toc;
        end
        
        lat = gps.lat(1:param.dec:end);
        lon = gps.lon(1:param.dec:end);
        % write out file
        if idx_size(1) > 1 % If more than one param.fn entry exists
          [a filename b] = fileparts(param.fn{fn_idx});
        else % If only one param.fn entry exists
          [a filename b] = fileparts(param.fn);
        end
        
        if ~ispc
          out_path = strcat('/',out_folder,filename,'.csv');
        else
          out_path = strcat(out_folder,filename,'.csv');
        end
        
        
        % Open and output CSV file
        output_file = strcat(filename,'.csv');
        fprintf('Writing out file %s and plotting map...\n',output_file);
        fid = fopen(out_path,'w');
        fprintf(fid,'%s,%s\n',header_lat,header_lon);
        for data_idx = 1:length(lat)
          fprintf(fid,'%f,%f\n',lat(data_idx),lon(data_idx));
          [x,y] = projfwd(proj,lat(data_idx),lon(data_idx));
          plot(x,y,'b-'); hold on;
        end
        fclose(fid);
      end
      
    else % One CSV file for any number of inputs
      fprintf('Processing %d file(s)...\n',idx_size(1));
      lat = {};
      lon = {};
      gps = {};
      for fn_idx = 1:length(param.fn)
        fprintf('Processing file %d of %d ...\n',fn_idx,idx_size(1));
        
        idx_size = size(param.fn);
        if idx_size(1) > 1 % If more than one param.fn entry exists
          tic; gps = read_gps_applanix(param.fn{fn_idx},param); toc;
        else % If only one param.fn entry exists
          tic; gps = read_gps_applanix(param.fn,param); toc;
        end
        lat{fn_idx} = gps{fn_idx}.lat(1:param.dec:end);
        lon{fn_idx} = gps{fn_idx}.lon(1:param.dec:end);
      end
      [a b c] = fileparts(csv_out_path);
      output_file = strcat(b,c);
      fprintf('Writing out file %s and plotting map...\n',output_file);
      fid = fopen(csv_out_path,'w');
      fprintf(fid,'%s,%s\n',header_lat,header_lon);
      for data_idx = 1:length(lat)
        for coord_idx = 1:length(lat{data_idx})
          fprintf(fid,'%f,%f\n',lat{data_idx}(coord_idx),lon{data_idx}(coord_idx));
          [x,y] = projfwd(proj,lat{data_idx}(coord_idx),lon{data_idx}(coord_idx));
          plot(x,y,'b-'); hold on;
        end
      end
      fclose(fid);
    end
    
  case 4
    %  elseif param.input_type == 4
    %% ATM (.traj) to CSV
    header_lat = 'LAT';
    header_lon = 'LON';
    idx_size = size(param.fn);
    
    % PARAMETERS
    % param.headerlines = 0;
    % param.delimiter_type = ',';
    % param.input_format = '%f%f%f%f%f%f%f%f';
    if isempty(param.time_reference)
      param.time_reference = 'utc';
    else
    end
    
    if isempty(param.dec)
      param.dec = 1000;
    else
    end
    
    out_folder = fileparts(csv_out_path);
    if ~ispc
      out_folder = strcat(out_folder,'/');
    else
      out_folder = strcat(out_folder,'\');
    end
    
    % Build structure and write file(s)
    if param.csv_merge == false % One CSV file for every input
      fprintf('Processing %d file(s)...\n',idx_size(1));
      for fn_idx = 1:length(param.fn)
        clear gps a filename b out_path;
        fprintf('Processing file %d of %d ...\n',fn_idx,length(param.fn));
        
        idx_size = size(param.fn);
        if idx_size(1) > 1 % If more than one param.fn entry exists
          tic; gps = read_gps_traj(param.fn{fn_idx},param); toc;
        else % If only one param.fn entry exists
          tic; gps = read_gps_traj(param.fn,param); toc;
        end
        
        lat = gps.lat(1:param.dec:end);
        lon = gps.lon(1:param.dec:end);
        
        % write out file
        if idx_size(1) > 1 % If more than one param.fn entry exists
          [a filename b] = fileparts(param.fn{fn_idx});
        else % If only one param.fn entry exists
          [a filename b] = fileparts(param.fn);
        end
        
        if ~ispc
          out_path = strcat('/',out_folder,filename,'.csv');
        else
          out_path = strcat(out_folder,filename,'.csv');
        end
        
        
        % Open and output CSV file
        output_file = strcat(filename,'.csv');
        fprintf('Writing out file %s and plotting map...\n',output_file);
        fid = fopen(csv_out_path,'w');
        fprintf(fid,'%s,%s\n',header_lat,header_lon);
        for data_idx = 1:length(lat)
          fprintf(fid,'%f,%f\n',lat(data_idx),lon(data_idx));
          [x,y] = projfwd(proj,lat{data_idx},lon{data_idx});
          plot(x,y,'b-'); hold on;
        end
        fclose(fid);
      end
    else % One CSV file for any number of inputs
      fprintf('Processing %d file(s)...\n',idx_size(1));
      lat = {};
      lon = {};
      gps = {};
      for fn_idx = 1:length(param.fn)
        fprintf('Processing file %d of %d ...\n',fn_idx,length(param.fn));
        idx_size = size(param.fn);
        if idx_size(1) > 1 % If more than one param.fn entry exists
          tic; gps = read_gps_traj(param.fn{fn_idx},param); toc;
        else % If only one param.fn entry exists
          tic; gps = read_gps_traj(param.fn,param); toc;
        end
        lat{fn_idx} = gps{fn_idx}.lat(1:param.dec:end);
        lon{fn_idx} = gps{fn_idx}.lon(1:param.dec:end);
      end
      [a b c] = fileparts(csv_out_path);
      output_file = strcat(b,c);
      fprintf('Writing out file %s and plotting map...\n',output_file);
      fid = fopen(csv_out_path,'w');
      fprintf(fid,'%s,%s\n',header_lat,header_lon);
      for data_idx = 1:length(lat)
        for coord_idx = 1:length(lat{data_idx})
          fprintf(fid,'%f,%f\n',lat{data_idx}(coord_idx),lon{data_idx}(coord_idx));
          [x,y] = projfwd(proj,lat{data_idx}(coord_idx),lon{data_idx}(coord_idx));
          plot(x,y,'b-'); hold on;
        end
      end
      fclose(fid);
    end
    
    
  case 5
    %   elseif param.input_type == 5
    %% Litton 100 (ln100.asc) to CSV
    header_lat = 'LAT';
    header_lon = 'LON';
    idx_size = size(param.fn);
    
    % PARAMETERS
    if isempty(param.time_reference)
      param.time_reference = 'utc';
    else
    end
    
    if isempty(param.dec)
      param.dec = 250;
    else
    end
    
    out_folder = fileparts(csv_out_path);
    if ~ispc
      out_folder = strcat(out_folder,'/');
    else
      out_folder = strcat(out_folder,'\');
    end
    
    if param.csv_merge == false % One CSV file for every input
      fprintf('Processing %d file(s)...\n',idx_size(1));
      % Build structure and write file(s)
      for fn_idx = 1:length(param.fn)
        clear gps a filename b out_path;
        fprintf('Processing file %d of %d ...\n',fn_idx,length(param.fn));
        idx_size = size(param.fn);
        if idx_size(1) > 1 % If more than one param.fn entry exists
          tic; gps = read_gps_litton(param.fn{fn_idx},param); toc;
        else % If only one param.fn entry exists
          tic; gps = read_gps_litton(param.fn,param); toc;
        end
        lat = gps.lat(1:param.dec:end);
        lon = gps.lon(1:param.dec:end);
        
        % write out file
        if idx_size(1) > 1 % If more than one param.fn entry exists
          [a filename b] = fileparts(param.fn{fn_idx});
        else % If only one param.fn entry exists
          [a filename b] = fileparts(param.fn);
        end
        
        if ~ispc
          out_path = strcat('/',out_folder,filename,'.csv');
        else
          out_path = strcat(out_folder,filename,'.csv');
        end
        
        % Open and output CSV file
        output_file = strcat(filename,'.csv');
        fprintf('Writing out file %s and plotting map...\n',output_file);
        fid = fopen(csv_out_path,'w');
        fprintf(fid,'%s,%s\n',header_lat,header_lon);
        for data_idx = 1:length(lat)
          fprintf(fid,'%f,%f\n',lat(data_idx),lon(data_idx));
          [x,y] = projfwd(proj,lat(data_idx),lon(data_idx));
          plot(x,y,'b-'); hold on;
        end
        fclose(fid);
      end
    else % One CSV file for any number of inputs
      fprintf('Processing %d file(s)...\n',idx_size(1));
      lat = {};
      lon = {};
      gps = {};
      for fn_idx = 1:length(param.fn)
        fprintf('Processing file %d of %d ...\n',fn_idx,length(param.fn));
        idx_size = size(param.fn);
        if idx_size(1) > 1 % If more than one param.fn entry exists
          tic; gps = read_gps_litton(param.fn{fn_idx},param); toc;
        else % If only one param.fn entry exists
          tic; gps = read_gps_litton(param.fn,param); toc;
        end
        lat{fn_idx} = gps{fn_idx}.lat(1:param.dec:end);
        lon{fn_idx} = gps{fn_idx}.lon(1:param.dec:end);
      end
      [a b c] = fileparts(csv_out_path);
      output_file = strcat(b,c);
      fprintf('Writing out file %s and plotting map...\n',output_file);
      fid = fopen(csv_out_path,'w');
      fprintf(fid,'%s,%s\n',header_lat,header_lon);
      for data_idx = 1:length(lat)
        for coord_idx = 1:length(lat{data_idx})
          fprintf(fid,'%f,%f\n',lat{data_idx}(coord_idx),lon{data_idx}(coord_idx));
          [x,y] = projfwd(proj,lat{data_idx}(coord_idx),lon{data_idx}(coord_idx));
          plot(x,y,'b-'); hold on;
        end
      end
      fclose(fid);
    end
  otherwise
    %   else
    fprintf('Invalid param.input_type selected.\n\n');
end
fprintf('Done ... %s\n',datestr(now,'HH:MM:SS'));