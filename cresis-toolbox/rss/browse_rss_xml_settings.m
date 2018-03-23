clc
clearvars
%% Constants
version_num = '4.0';
radar_name = 'hfrds2';
season = '2016_Greenland_TO';
season_name = '2016_Greenland_TOdtu';
data_dir = 'Z:/';
if strcmpi(radar_name,'hfrds2')
    rad_dir = 'HF_Sounder/';
end
geotiff_fn = 'Greenland/Landsat-7/Greenland_natural_150m';
file_version = 412;

all_dir_bool = false; %Turn to true to search all folders in the directory
    %Used when all_dir_bool is false
flt_dates = {'20161101', '20161102'};%, '20161107', '20161108', '20161110', '20161111', '20161112'}; 

chk_dir = fullfile(data_dir,rad_dir,season);%Directory to check for flights
if all_dir_bool
    %Check for the dirs
    chk_me = dir(chk_dir);
    %Extract the flight dates (format: YYYYMMDD)
    load_id = 1; %Initialized for loading the flight dates
    flt_dates = {};
    for chk_id = 3:length(chk_me) %3 is used to skip the . and .. directories
        flt_dates{load_id} = chk_me(chk_id).name(1:8);
        load_id = load_id+1;
    end
    %Extract the unique dates
    flt_dates = unique(flt_dates);
end

if ~exist('flt_dates','var')
    error('No flight dates have been specified. Either insert them manually or search for them by setting the all_dir_bool variable to true')
elseif isempty(flt_dates)
    error('No flight dates have been specified or none were found in %s',chk_dir)
end
%% Build the vectors worksheet for an rss based radar
for flt_id = 1:length(flt_dates)
%Input constant values
    %Date (YYYYMMDD format) This is just the day of interest
        flt_date = flt_dates{flt_id}; %flight date
                    
    %base_dir contains a string which should go to the base of the 
        %dataset; this string should be identical for every segment
        base_dir = fullfile(data_dir,rad_dir,season);
        base_dir_load = fullfile(rad_dir,season);
        
    %Given a base_dir and a date, find all of the folders that match
    full_dir = fullfile(base_dir,flt_date);
    dir_names = dir(sprintf('%s*',full_dir));
    
    if isempty(dir_names)
        error('This radar does not have records corresponding to that date.')
    end
    
    data_fns_chk = get_filenames(fullfile(base_dir,dir_names(1).name),flt_date,'','.dat');
      if isempty(data_fns_chk) %Mainly for 20161101
          %Check for segment folders
          dir_names= dir(sprintf('%s/seg*',full_dir));
          base_dir = full_dir;
          base_dir_load = fullfile(rad_dir,season,flt_date);%May need to be changed
      end
    %Generate the file prefixes
    for dir_id = 1:length(dir_names)
        %adc_folder_name contains a string which gives the rest of the path
            adc_folder_name = dir_names(dir_id).name;        

        %Using base_dir, adc_folder_name, and prefix in the get_filenames
            %function, the file ids can be found

          % Get raw data files associated with this directory
          data_fns = get_filenames(fullfile(base_dir,adc_folder_name),flt_date,'','.dat');

          fn_datenums = [];

          % Get the date information out of the filename
          for data_fn_idx = 1:length(data_fns)
            fn = data_fns{data_fn_idx}; %Define current filename
            [path, fn, ext] = fileparts(fn); %Separate the path and filename
            fn_parts = textscan(fn,'%s%s%s','Delimiter','_');%Get the time string for each individual test
            time_str = fn_parts{2}{1};
            if data_fn_idx ==1
                if ~exist('idx','var')
                    idx =1;%initialize loader for start and stop. Also the number of segments
                else
                    idx = idx+1;
                end
                curr_time_str{idx} = time_str;
                flt_segs{idx} = flt_date;
                file_prefix{idx} = sprintf('%s_%s_',flt_date,time_str);
                adc_fold_vec{idx} = adc_folder_name; %Load this for later
                %Determine the starting and ending file ids
                file_fns = get_filenames(fullfile(base_dir,adc_folder_name),...
                    file_prefix{idx},'','.dat');
                start_idx(idx) = 1; %Should be one
                stop_idx(idx) = length(file_fns);
            else
                if ~strcmpi(time_str,curr_time_str{idx})
                    idx = idx+1;
                    flt_segs{idx} = flt_date;
                    file_prefix{idx} = sprintf('%s_%s_',flt_date,time_str);
                    curr_time_str{idx} = sprintf('%s_%s_',flt_date,time_str);       
                    adc_fold_vec{idx} = adc_folder_name; %Load this for later
                    curr_time_str{idx} = time_str;
                    %Determine the starting and ending file ids
                    file_fns = get_filenames(fullfile(base_dir,adc_folder_name),...
                        file_prefix{idx},'','.dat');
                    start_idx(idx) = 1;%Should be one
                    stop_idx(idx) = length(file_fns);
                end
            end
          end
    end
end

%% Load segment data structures
for seg_id = 1:idx
      %Put in the date
      seg_dat{seg_id}.date = flt_segs{seg_id};
      %Put in the segment number
      seg_dat{seg_id}.seg_id = seg_id;
      %Put in the start id
      seg_dat{seg_id}.start_idx = start_idx(seg_id);
      %Put in the stop id
      seg_dat{seg_id}.stop_idx = stop_idx(seg_id);
      %Put in the base_dir
      seg_dat{seg_id}.base_dir = base_dir_load;
      %Put in the completed adc_folder_name
      seg_dat{seg_id}.adc_folder_name = adc_fold_vec{seg_id};
      %Load the file_prefix
      seg_dat{seg_id}.file_prefix = file_prefix{seg_id};
      %Load the geotiff_fn
      seg_dat{seg_id}.geotiff_fn = geotiff_fn;
      %Load the file version
      seg_dat{seg_id}.file_version = file_version;
end

%% Print out the cmd worksheet parameters
%Print the cmd sheet headers
cmd_headers
%Print each line 
for seg_id = 1:idx
      %Print the line
      cmd_line(seg_dat{seg_id});
end

%% Print out the vector worksheet parameters
%Print the vectors sheet headers
vectors_headers
%Print each line 
for seg_id = 1:idx
      %Print the line
      vectors_line(seg_dat{seg_id});
end

%% Print out the records worksheet parameters
%Print the records sheet headers
records_headers
%Print each line 
for seg_id = 1:idx
      %Print the line
      records_line(seg_dat{seg_id});
end
