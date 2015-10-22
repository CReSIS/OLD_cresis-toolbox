% Script Information/Details
% see also run_cross_over_analysis

function cross_over_analysis(param)
% cross_over_analysis(param)
%
% This function loads fin, finds the cross-overs, and writes the output
% to a CSV and MAT file.
%
% param = struct controlling cross-over analysis
%  .fin = input file name or path
%  .data_type.type = is an integer (1 or 2)
%    1: CSV arbitrary format (fin is a single CSV file)
%       Also have to specify param.datatype struct with these fields:
%        .format = Default is '%f%f%f%f%f%s%f%f%f' (textscan argument)
%        .delim = Default is ',' (textscan argument)
%        .headerlines = Default is 1 (textscan argument)
%        .col_lat = Default is 1
%        .col_lon = Default is 2
%        .col_time = Default is 3
%        .col_thickness = Default is 4
%        .col_frame = Default is 6
%    2: MAT layer file format (fin is a path or a directory tree
%       containing the layer files)
%       Also have to specify param.datatype struct with these fields:
%        .segs = Default is {} (empty cell array does all segments)
%                e.g. segs = {'20091118_01','20091102_01'}
%  .fout = output file name (a .csv and .mat file will be created)
%  .min_sample_spacing = user-determined along-track "delta-x" between
%     points. Data is down-sampled to this spacing. Down-sampling is
%     performed to make cross over analysis run faster.
%     Default is 100 m. Leave undefined or empty for default.
%  .scale_check = distance check. It is a scale factor multiplied by
%     min_sample_spacing, so if min_sample_spacing is 100 m and
%     scale_check is 5, the distance check for cross overs is 500 m.
%     Default is 5. Leave undefined or empty for default.
%  .min_idx_sep = index check
%     Default is 20 indices. Leave undefined or empty for default.
%  .UTC_check = time check
%     Default is 100 sec. Leave undefined or empty for default.
%  .run_sort = Boolean to run sort of points based on geographic and
%     time information for each point.
%     Run points through a sorter which sorts them by flight line and
%     in the order that the data was collected. This is useful is the
%     data points you are loading are randomly sorted (e.g. an export
%     from ARC will typically produce unsorted points). This is NOT
%     needed when loading .mat layer files.
%     Default is to not run. Leave undefined or empty for default.
%  .debug_level = 1 is default
%
% Examples: See run_cross_over_analysis
%
% Author: Brady Maasen, John Paden
%
% See also plot_cross_overs, run_cross_over_analysis

fprintf('==========================================================\n');
fprintf('==========================================================\n');

%% Check Input Arguments

if ~isfield(param,'min_sample_spacing') || isempty(param.min_sample_spacing)
  param.min_sample_spacing = 100;
end

if ~isfield(param,'scale_check') || isempty(param.scale_check)
  param.scale_check = 5;
end

if ~isfield(param,'min_idx_sep') || isempty(param.min_idx_sep)
  param.min_idx_sep = 20;
end

if ~isfield(param,'UTC_check') || isempty(param.UTC_check)
  param.UTC_check = 100;
end

if ~isfield(param,'run_sort') || isempty(param.run_sort)
  param.run_sort = false;
end

if ~isfield(param,'debug_level') || isempty(param.debug_level)
  param.debug_level = 1;
end

if ~isfield(param,'data_type') || isempty(param.data_type)
  param.data_type.type = 2;
  param.data_type.format = '%f%f%f%f%f%s%f%f%f';
  param.data_type.delim = ',';
  param.data_type.headerlines = 1;
  param.data_type.col_lat = 1;
  param.data_type.col_lon = 2;
  param.data_type.col_time = 3;
  param.data_type.col_thickness = 4;
  param.data_type.col_frame = 6;
else
  if ~isfield(param.data_type,'type') || isempty(param.data_type.type)
    param.data_type.type = 2;
  end
  if ~isfield(param.data_type,'segs') || isempty(param.data_type.segs)
    param.data_type.segs = {};
  end
  if ~isfield(param.data_type,'format') || isempty(param.data_type.format)
    param.data_type.format = '%f%f%f%f%f%s%f%f%f';
  end
  if ~isfield(param.data_type,'delim') || isempty(param.data_type.delim)
    param.data_type.delim = ',';
  end
  if ~isfield(param.data_type,'headerlines') || isempty(param.data_type.headerlines)
    param.data_type.headerlines = 1;
  end
  if ~isfield(param.data_type,'col_lat') || isempty(param.data_type.col_lat)
    param.data_type.col_lat = 1;
  end
  if ~isfield(param.data_type,'col_lon') || isempty(param.data_type.col_lon)
    param.data_type.col_lon = 2;
  end
  if ~isfield(param.data_type,'col_time') || isempty(param.data_type.col_time)
    param.data_type.col_time = 3;
  end
  if ~isfield(param.data_type,'col_thickness') || isempty(param.data_type.col_thickness)
    param.data_type.col_thickness = 4;
  end
  if ~isfield(param.data_type,'col_frame') || isempty(param.data_type.col_frame)
    param.data_type.col_frame = 6;
  end
end

tic
physical_constants;

%% Get the Data

fprintf('Loading data (%.1f sec)\n',toc);

if param.data_type.type == 2
  fns = get_filenames(param.fin,'Data','','*.mat',struct('recursive',1));
  if ~isempty(param.data_type.segs)
    keep_mask = boolean(zeros(size(fns)));
    for file_idx = 1:length(fns)
      [fn_dir fn_name] = fileparts(fns{file_idx});
      day_seg = fn_name(6:16);
      if ~isempty(strmatch(day_seg,param.data_type.segs))
        keep_mask(file_idx) = 1;
      end
    end
    fns = fns(keep_mask);
  end
  
  lat = [];
  lon = [];
  UTC = [];
  thickness = [];
  seg_frame = {};
  for file_idx = 1:length(fns)
    fprintf('  Loading file %d of %d (%.1f sec)\n', file_idx, length(fns), toc);
    lay = load(fns{file_idx});
    lat = cat(1,lat,lay.Latitude.');
    lon = cat(1,lon,lay.Longitude.');
    UTC = cat(1,UTC,lay.GPS_time.');
    thickness = cat(1,thickness, ...
      ((lay.layerData{2}.value{2}.data - lay.layerData{1}.value{2}.data)*(c/2)/sqrt(er_ice)).');
    [fn_dir fn_name] = fileparts(fns{file_idx});
    frm_id = {fn_name(6:20)};
    seg_frame = cat(1,seg_frame,repmat(frm_id,[length(lay.Latitude) 1]));
  end
elseif param.data_type.type == 1
  fid = fopen(param.fin,'r');
  X = textscan(fid,param.data_type.format,'delimiter',param.data_type.delim, ...
    'headerlines',param.data_type.headerlines); 
  lat = X{param.data_type.col_lat};
  lon = X{param.data_type.col_lon};
  UTC = X{param.data_type.col_time};
  thickness = X{param.data_type.col_thickness};
  seg_frame = X{param.data_type.col_frame};

elseif param.data_type.type == 3 || param.data_type.type == 4
  fid = fopen(param.fin,'r');
  X = textscan(fid,param.data_type.format,'delimiter',param.data_type.delim, ...
    'headerlines',param.data_type.headerlines);
  lat = X{param.data_type.col_lat};
  lon = X{param.data_type.col_lon};
  UTC = X{param.data_type.col_time};
  thickness = X{param.data_type.col_thickness};
  seg_frame = repair_frame(X{param.data_type.col_YYYYMMDD},X{param.data_type.col_segment},X{param.data_type.col_frame});
else
  fprintf('Data format %d not supported\n', param.data_type.type);
  return;
end


%% Process the Data
% remove NaNs in the data
lat_nan = find(isnan(lat) == 1);
lon_nan = find(isnan(lon) == 1);
nan_idxs = union(lat_nan,lon_nan);
if(~isempty(nan_idxs))
    fprintf('Warning: NaN values found in data and are removed.\n');
    lat(nan_idxs) =[];
    lon(nan_idxs) =[];
    UTC(nan_idxs) =[];
    thickness(nan_idxs) = [];
    seg_frame(nan_idxs) = [];
end

% Convert lat/lon to meters (Use UTM or STEREOGRAPHIC)
if  ~isempty(find(lat<-80)) || ~isempty(find(lat>84))
	% Outside of UTM range.
	[x,y,mstruct] = geodetic_to_stereographic(lat,lon);
else
	% Within UTM range.
	[x,y,mstruct] = geodetic_to_utm(lat,lon);
end

% ------------------------------------------------------------------------%
% Sorting original data by source or frame
% Used when data is not sorted by segment
% Necessary when processing data from ArcGIS
if param.run_sort
  fprintf('Sorting data (%.1f sec)\n',toc);
  
  
  % ------------------------------------------------------------------------%
  sources = unique(seg_frame);
  for cur_source = 1:length(sources)
    isource = strmatch(sources{cur_source},seg_frame);%indices of rows w/ this source
    fprintf('  Sorting source %s (%.1f sec)\n', sources{cur_source}, toc);
    inc = 0;
    sort_mat = zeros(size(isource));
    done_mask = zeros(size(isource));
    keep_going = 0;
    if param.debug_level >= 2
      figure(1); clf;
      plot(x(isource),y(isource),'b.');
      plot_idx = 1;
    end
    while ~all(done_mask) % ends when sort_mat is full
      if ~keep_going %lines are executed at the beginning of a flight line
        inc = inc+1;
        cur_idx = find(done_mask == 0,1,'first');
        sort_mat(inc) = isource(cur_idx);
        done_mask(cur_idx) = 1;
        keep_going = 1;
        if param.debug_level >= 2
          hold on;
          plot(x(sort_mat(plot_idx:inc-1)),y(sort_mat(plot_idx:inc-1)),'r.-');
          keyboard
          plot(x(sort_mat(plot_idx:inc-1)),y(sort_mat(plot_idx:inc-1)),'g.-');
          fprintf('============================================\n');
          hold off;
          plot_idx = inc;
        end
      else
        cur_idx = next_idx;
      end
      dist = sqrt((x(isource(cur_idx))-x(isource)).^2 + (y(isource(cur_idx))-y(isource)).^2);
      bool_dist = dist < 2.5*param.min_sample_spacing*param.scale_check; %idxs that fit dist criteria
      time = UTC(isource) - UTC(isource(cur_idx));
      bool_time = time >= 0 & time < 30; %idxs that fit time criteria
      good_idxs = find(bool_dist & bool_time & ~done_mask); %idxs that fit criteria and haven't been picked yet
      if isempty(good_idxs)
        keep_going = 0;
      else
        [mtime next_idx] = min(time(good_idxs));
        if mtime < 0
          keyboard
        end
        next_idx = good_idxs(next_idx);
        if param.debug_level >= 2
          fprintf('%10d %10.2f\n', cur_idx, dist(next_idx));
        end
        %save the info
        inc = inc+1;
        sort_mat(inc) = isource(next_idx);
        done_mask(next_idx) = 1;
      end
    end
    seg_frame(isource) = seg_frame(sort_mat);
    x(isource) = x(sort_mat);
    y(isource) = y(sort_mat);
    UTC(isource) = UTC(sort_mat);
    thickness(isource) = thickness(sort_mat);
  end
  [seg_frame sort_idxs] = sort(seg_frame);
  x = x(sort_idxs);
  y = y(sort_idxs);
  UTC = UTC(sort_idxs);
  thickness = thickness(sort_idxs);
  
  %% Debugging Code
  
  if 0
    % % debugging code
    % idxs = strmatch('Line90',seg_frame)
    % mask = zeros(size(idxs));
    % mask to know which ones you have done, while all(mask)
    valid_mask = zeros(size(seg_frame));
    for idx = 1:length(seg_frame)
      tmp = strfind(seg_frame{idx},'valid');
      if ~isempty(tmp)
        valid_mask(idx) = 1;
      end
    end
    ids = strmatch('1997_3_valid',seg_frame);
    ids = strmatch('1997_f_valid',seg_frame);
    ids = strmatch('1998_3_valid',seg_frame);
    ids = strmatch('2001_3_valid',seg_frame);
    ids = strmatch('2005_3_valid',seg_frame);
    ids = strmatch('2005_f_valid',seg_frame);
    ids = strmatch('2006_3_valid',seg_frame);
    seg_frame1 = seg_frame(ids);
    id_col = X{7};
    id1 = id_col(sort_mat);
    figure;
    plot(id1)
    hold on;
    plot(id1,'r.')
    hold off;
    
    plot(lon(ids),lat(ids))
    
    sort_string = seg_frame;
    for idx = 1:length(sort_string)
      sort_string{idx} = [sort_string{idx} UTC_string(idx,:)];
    end
  end
  
  %UTC_string = vec2mat(sprintf('%09.2f',UTC),9);
  %   UTC_string = vec2mat(sprintf('%07f',X{7}),7);
  %   sort_string = seg_frame;
  %   for idx = 1:length(sort_string)
  %     sort_string{idx} = [sort_string{idx} UTC_string(idx,:)];
  %   end
  %   [sort_string sort_idxs] = sort(sort_string);
  %   lat = lat(sort_idxs);
  %   lon = lon(sort_idxs);
  %   UTC = UTC(sort_idxs);
  %   thickness = thickness(sort_idxs);
end

%% Down-Sample Data

fprintf('Down-sample data (%.1f sec)\n',toc);

% Find the along-track distance
along_track = [0; cumsum(sqrt(diff(x).^2 + diff(y).^2))];
good_idxs = zeros(size(along_track));
good_idxs(1) = 1;
last_idx = 1;
for jj = 2:length(along_track)
  if along_track(jj)-along_track(last_idx) >= param.min_sample_spacing
    good_idxs(jj) = 1;
    last_idx = jj;
  end
end

good_idxs = find(good_idxs);

x = x(good_idxs);
y = y(good_idxs);
% along_track = along_track(good_idxs); % Debugging
UTC = UTC(good_idxs);
thickness = thickness(good_idxs);
seg_frame = seg_frame(good_idxs);


%% Find the indexes where cross-overs happen

fprintf('Find all potential cross-overs (%.1f sec)\n',toc);

poi(1) = round(length(x)/10);
poi(2) = poi(1)*2;
poi(3) = poi(1)*4;
poi(4) = poi(1)*6;
poi(5) = poi(1)*8;
poi(6) = poi(1)*9;
% Memory allocation block_size (forces Matlab to grow matrices block_size
% at a time so not too much time is wasted in memory allocation routine)
block_size = 1e6;
cross_idxs = zeros(block_size,2);
cross_dist = zeros(block_size,1);
count = 0;
prev_time = 0;
for jj = 1:length(x)
  if jj == poi(1)
    fprintf('  About 10%% complete (%.1f sec)\n', toc);
  elseif jj == poi(2)
    fprintf('  About 20%% complete (%.1f sec)\n', toc);
  elseif jj == poi(3)
    fprintf('  About 40%% complete (%.1f sec)\n', toc);
  elseif jj == poi(4)
    fprintf('  About 60%% complete (%.1f sec)\n', toc);
  elseif jj == poi(5)
    fprintf('  About 80%% complete (%.1f sec)\n', toc);
  elseif jj == poi(6)
    fprintf('  About 90%% complete (%.1f sec)\n', toc);
  end
  cur_dist = (x(jj)-x(jj+1:end)).^2 + (y(jj)-y(jj+1:end)).^2; %vector of current dist's
  cur_dist_bool = cur_dist < (param.min_sample_spacing*param.scale_check)^2;
  cur_dist_idxs = jj + find(cur_dist_bool); %indexes that fit criteria
  for kk = 1:length(cur_dist_idxs) %scan all close enough in dist
    t1 = toc;
    if t1 >= prev_time + 10
      fprintf('    index is %d/%d (%.1f sec)\n',jj,length(x), t1);
      prev_time = t1;
    end
    if abs(UTC(cur_dist_idxs(kk))-UTC(jj))>param.UTC_check
      count = count + 1;
      if count > size(cross_idxs,1)
        cross_idxs(size(cross_idxs,1)+block_size,2) = 0;
        cross_dist(size(cross_dist,1)+block_size,1) = 0;
      end
      cross_idxs(count,:) = [jj,cur_dist_idxs(kk)]; %save the two indexes
      cross_dist(count) = cur_dist(cur_dist_idxs(kk)-jj); %save the distance
    end
  end
end

%taking off the zeros in cross_idxs and cross_dist
cross_idxs = cross_idxs(1:count,:);
cross_dist = cross_dist(1:count);

% get a plot of what we've found
if param.debug_level >= 2
  figure;plot(x,y,'.k');
  hold on;
  plot(x(cross_idxs(1:end)),y(cross_idxs(1:end)),'r+');
end


%% Find the indexes where cross-overs happen

fprintf('Culling out bad cross overs (%.1f sec)\n',toc);

cross_over = cross_idxs;
dist = cross_dist;
% min_tuple = zeros(1,2);
done_mask = zeros(size(cross_over,1),1);
for jj = 1:size(cross_over,1)-1
  if (~done_mask(jj))
    min_dist = inf;
    kk = jj;
    while kk <= size(cross_over,1) && (~done_mask(kk) && cross_over(kk,1) <= cross_over(jj,1) + param.min_idx_sep || done_mask(kk,1))%abs(UTC(cross_over(kk,1)) - UTC(cross_over(jj,1))) <= 30)%
      %if ~done_mask(kk,1) && abs(UTC(cross_over(kk,2)) - UTC(cross_over(jj,2))) <= 30
      if ~done_mask(kk,1) && abs(cross_over(kk,2) - cross_over(jj,2)) <= param.min_idx_sep
        if dist(kk) < min_dist
          min_dist = dist(kk);
          min_idx = kk;
        end
        done_mask(kk) = 1;
      end
      kk = kk +1;
      %[cross_over(jj:kk+10,:) done_mask(jj:kk+10)] % For debugging
    end
    done_mask(min_idx) = 2;
  end
end
good_ones = find(done_mask == 2);
dist = dist(good_ones);
cross_over = cross_over(good_ones,:);

if param.debug_level >= 2
  figure;plot(x,y,'.k');
  hold on;
  plot(x(cross_over),y(cross_over),'r+')
end


%% Preparing outputs and writing to output files

fprintf('Creating outputs and output files (%.1f sec)\n',toc);

% Finding angle between crossing lines
cross_angle = zeros(size(cross_over,1),1);
for row = 1:size(cross_over,1)
  if cross_over(row,1) < length(x)
    x_coord1 = x(cross_over(row,1)+1) - x(cross_over(row,1));
    y_coord1 = y(cross_over(row,1)+1) - y(cross_over(row,1));
  else
    x_coord1 = x(cross_over(row,1)) - x(cross_over(row,1)-1);
    y_coord1 = y(cross_over(row,1)) - y(cross_over(row,1)-1);
  end
  vec1 = [x_coord1,y_coord1];
  
  if cross_over(row,2) < length(x)
    x_coord2 = x(cross_over(row,2)+1) - x(cross_over(row,2));
    y_coord2 = y(cross_over(row,2)+1) - y(cross_over(row,2));
  else
    x_coord2 = x(cross_over(row,2)) - x(cross_over(row,2)-1);
    y_coord2 = y(cross_over(row,2)) - y(cross_over(row,2)-1);
  end
  vec2 = [x_coord2,y_coord2];
  cross_angle(row) = acosd(dot(vec1,vec2)/sqrt(dot(vec1,vec1)*dot(vec2,vec2)));
end
cross_angle = real(cross_angle);
dist = sqrt(dist);

% =========================================================================
% Get our lat/lon back
% =========================================================================
[lat lon] = minvtran(mstruct,x,y);
% outfile = fopen(fout,'a');
Latitude = lat(cross_over);
Longitude = lon(cross_over);
Frame_ID = seg_frame(cross_over);
Thickness = thickness(cross_over);
UTC_time = UTC(cross_over);
Error = Thickness(:,1) - Thickness(:,2);
absError = abs(Error);

%% Saving and Printing Data

csv_header = cell(1,14);
csv_header{1} = 'LATA';
csv_header{2} = 'LONA';
csv_header{3} = 'THICKA';
csv_header{4} = 'FRAMEA';
csv_header{5} = 'TIMEA';
csv_header{6} = 'LATB';
csv_header{7} = 'LONB';
csv_header{8} = 'THICKB';
csv_header{9} = 'FRAMEB';
csv_header{10} = 'TIMEB';
csv_header{11} = 'ABSDIFF';
csv_header{12} = 'ANGLE';
csv_header{13} = 'DISTANCE';

fn_csv = strcat(param.fout,'.csv');
fn_csv_dir = fileparts(fn_csv);
% Check for the existence of the output directory
if ~exist(fn_csv_dir,'dir')
  warning('Output directory %s does not exist. Creating...',fn_csv_dir);
  mkdir(fn_csv_dir);
end
fprintf('  Saving %s\n', fn_csv);

fid_csv = fopen(fn_csv,'w');

fprintf(fid_csv,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',csv_header{1},csv_header{2},csv_header{3},csv_header{4},...
  csv_header{5},csv_header{6},csv_header{7},csv_header{8},csv_header{9},csv_header{10},csv_header{11},csv_header{12},csv_header{13});

csv_combine = zeros(length(Longitude),13);

csv_combine(:,1) = Latitude(:,1);
csv_combine(:,2) = Longitude(:,1);
csv_combine(:,3) = Thickness(:,1);
csv_combine(:,4) = str2double(Frame_ID(:,1));
csv_combine(:,5) = UTC_time(:,1);
csv_combine(:,6) = Latitude(:,2);
csv_combine(:,7) = Longitude(:,2);
csv_combine(:,8) = Thickness(:,2);
csv_combine(:,9) = str2double(Frame_ID(:,2));
csv_combine(:,10) = UTC_time(:,2);
csv_combine(:,11) = absError;
csv_combine(:,12) = cross_angle;
csv_combine(:,13) = dist;

for jj = 1:length(Longitude)
  fprintf(fid_csv,'%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n',csv_combine(jj,1),csv_combine(jj,2),...
    csv_combine(jj,3),csv_combine(jj,4),csv_combine(jj,5),csv_combine(jj,6),csv_combine(jj,7),csv_combine(jj,8),...
    csv_combine(jj,9),csv_combine(jj,10),csv_combine(jj,11),csv_combine(jj,12),csv_combine(jj,13));
end

fclose(fid_csv);

mat_fout = strcat(param.fout,'.mat');
fprintf('  Saving %s\n', mat_fout);
save(mat_fout,'Latitude','Longitude','Frame_ID','Thickness','UTC_time','cross_angle','dist','absError');


%Calculate and save RMS and RMeS
RMS = sqrt(mean(Error.^2));
RMeS = sqrt(median(Error.^2));

txt_fout = strcat(param.fout,'.txt');
fprintf('  Saving %s\n', txt_fout);

errorFout = strcat(param.fout,'_Errors','.txt');
fid_error = fopen(errorFout,'wt');
fprintf(fid_error,'RMS,RMeS\n');
fprintf(fid_error,'%2.4f,%2.4f\n',RMS,RMeS);
fclose(fid_error);

%End Script
fprintf('Done (%.1f sec)\n', toc);

return;

