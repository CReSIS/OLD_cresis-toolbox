function post_gis_qc(zip_fn)
% Runs a post Quality Check for a GIS post zipped folder.
%
% post_gis_qc(zip_fn)
%
% zip_fn: Absolute path to zipped GIS post folder.
%
% Author: Kyle Purdon

try
  
  fprintf('================================\n');
  fprintf('          POST GIS QC           \n');
  fprintf('================================\n');
  
  %% Unpackage and gather list of files
  
  % Get Base Folder and Filename
  [basepath,base_fn,tmp] = fileparts(zip_fn);
  basepath = strcat(basepath,filesep);
  
  fprintf('Unpackaging Project Folder ... ');
  project_fns = unzip(zip_fn,basepath);
  fprintf('%s\n',datestr(now,13));
  
  % Create Setup Files
  fprintf('Parsing Filenames ... ');
  error_messages = {};
  zip_folder_names = cell(1,length(project_fns));
  zip_file_names = cell(1,length(project_fns));
  for file_idx = 1: length(project_fns)
    sep_idxs = strfind(project_fns{file_idx},filesep);
    zip_folder_names{1,file_idx} = project_fns{file_idx}(sep_idxs(end-1)+1:sep_idxs(end)-1);
    zip_file_names{1,file_idx} = project_fns{file_idx}(sep_idxs(end)+1:length(project_fns{file_idx}));
  end
  fprintf('%s\n',datestr(now,13));
  
  %% Check for correct folder structure
  
  folders = unique(zip_folder_names);
  
  fprintf('Checking Folder structure ... ');
  if ~strcmp('boundaries',folders)
    error_messages{end+1} = 'Boundaries folder does not exist.';
  end
  if ~strcmp('errors',folders)
    error_messages{end+1} = 'Errors folder does not exist.';
  end
  if ~strcmp('flightlines',folders)
    error_messages{end+1} = 'Flightlines folder does not exist.';
  end
  if ~strcmp('grids',folders)
    error_messages{end+1} = 'Grids folder does not exist.';
  end
  if ~strcmp('preview_images',folders)
    error_messages{end+1} = 'Preview_Images folder does not exist.';
  end
  if ~strcmp('readme',folders)
    error_messages{end+1} = 'Readme folder does not exist.';
  end
  fprintf('%s\n',datestr(now,13));
  
  % Kill Block
  if ~isempty(error_messages)
    fprintf('ERROR MESSAGES:\n');
    for error_idx = 1:length(error_messages)
      fprintf('(%d) %s \n',error_idx,error_messages{error_idx});
    end
    error('Errors Found. See Messages Above.');
  end
  
  %% Check for correct file in folder structure
  fprintf('Seperating Files by Folder ... ');
  files_by_folders = {};
  for file_idx = 1:length(folders)
    files_by_folders{1,end+1} = strcmp(folders{file_idx},zip_folder_names);
  end
  % Note: Follow alphabetical Order.
  bound_files     = zip_file_names(files_by_folders{1});
  error_files     = zip_file_names(files_by_folders{2});
  fline_files     = zip_file_names(files_by_folders{3});
  grid_files      = zip_file_names(files_by_folders{4});
  preview_files   = zip_file_names(files_by_folders{5});
  readme_files    = zip_file_names(files_by_folders{6});
  fprintf('%s\n',datestr(now,13));
  
  fprintf('Checking # of Files by Folder ... ');
  if length(bound_files)  ~= 6
    error_messages{end+1} = 'Incorrect Number of Boundary Files.';
  end
  if length(error_files)  ~= 4
    error_messages{end+1} = 'Incorrect Number of Error Files.';
  end
  if length(fline_files)  ~= 7
    error_messages{end+1} = 'Incorrect Number of Flightline Files.';
  end
  if length(grid_files)   ~= 7
    error_messages{end+1} = 'Incorrect Number of Grid Files.';
  end
  if length(preview_files)~= 4
    error_messages{end+1} = 'Incorrect Number of Preview_Images Files.';
  end
  if length(readme_files) ~= 1
    error_messages{end+1} = 'Incorrect Number of README Files.';
  end
  fprintf('%s\n',datestr(now,13));
  
  % Kill Block
  if ~isempty(error_messages)
    fprintf('ERROR MESSAGES:\n');
    for error_idx = 1:length(error_messages)
      fprintf('(%d) %s \n',error_idx,error_messages{error_idx});
    end
    error('Errors Found. See Messages Above.');
  end
  
  %% Check for corrected file naming
  fprintf('Checking for Correct Files ... ');
  % Boundary Files
  good_bound_fns = {'_StudyArea.shp','_StudyArea.dbf','_StudyArea.prj','_StudyArea.sbn','_StudyArea.sbx','_StudyArea.shx'};
  for bfn_idx = 1:length(good_bound_fns)
    full_bound_fns = strcat(base_fn,good_bound_fns{bfn_idx});
    if ~strcmp(full_bound_fns,bound_files)
      error_messages{end+1} = sprintf('Filename Error On: [%s]',full_bound_fns);
    end
  end
  % Error Files
  good_error_fns = {'_Crossovers.csv','_Crossovers_Preview.png','_Crossovers.mat','_Crossovers_Errors.txt'};
  for efn_idx = 1:length(good_error_fns)
    full_error_fns = strcat(base_fn,good_error_fns{efn_idx});
    if ~strcmp(full_error_fns,error_files)
      error_messages{end+1} = sprintf('Filename Error On: [%s]',full_error_fns);
    end
  end
  % Flightline Files
  good_fline_fns = {'_Flightlines.txt','_Flightlines.shp','_Flightlines.dbf','_Flightlines.prj','_Flightlines.sbn','_Flightlines.sbx','_Flightlines.shx'};
  for flfn_idx = 1:length(good_fline_fns)
    full_fline_fns = strcat(base_fn,good_fline_fns{flfn_idx});
    if ~strcmp(full_fline_fns,fline_files)
      error_messages{end+1} = sprintf('Filename Error On: [%s]',full_fline_fns);
    end
  end
  % Grid Files
  good_grid_fns = {'_XYZGrid.txt','_bottom.txt','_surface.txt','_thickness.txt','_Bottom.prj','_Surface.prj','_Thickness.prj',};
  for gfn_idx = 1:length(good_grid_fns)
    full_grid_fns = strcat(base_fn,good_grid_fns{gfn_idx});
    if ~strcmpi(full_grid_fns,grid_files) % Ignore case, tool created.
      error_messages{end+1} = sprintf('Filename Error On: [%s]',full_grid_fns);
    end
  end
  % Preview_Image Files
  good_preview_fns = {'_Bottom_Preview.png','_Flightline_Preview.png','_Surface_Preview.png','_Thickness_Preview.png'};
  for pfn_idx = 1:length(good_preview_fns)
    full_preview_fns = strcat(base_fn,good_preview_fns{pfn_idx});
    if ~strcmp(full_preview_fns,preview_files)
      error_messages{end+1} = sprintf('Filename Error On: [%s]',full_preview_fns);
    end
  end
  % README Files
  good_readme_fns = {'_README.txt'};
  for rfn_idx = 1:length(good_readme_fns)
    full_readme_fns = strcat(base_fn,good_readme_fns{rfn_idx});
    if ~strcmp(full_readme_fns,readme_files)
      error_messages{end+1} = sprintf('Filename Error On: [%s]',full_readme_fns);
    end
  end
  fprintf('%s\n',datestr(now,13));
  
  % Kill Block
  if ~isempty(error_messages)
    fprintf('ERROR MESSAGES:\n');
    for error_idx = 1:length(error_messages)
      fprintf('(%d) %s \n',error_idx,error_messages{error_idx});
    end
    error('Errors Found. See Messages Above.');
  end
  
  %% File Content and Format validation
  fprintf('Checking File Formats:\n');
  fprintf('------------------------------\n')
  % Boundary Files
  % ======================================================
  fprintf('\tBoundary Files ... ');
  try
    bound_shp_fn = strcat(basepath,base_fn,'\boundaries\',full_bound_fns(1:length(full_bound_fns)-4),'.shp');
    bound_shp = shapeinfo(bound_shp_fn);
    if bound_shp.NumFeatures ~= 1
      error_messages{end+1} = 'Boundary cannot be a multi-part feature.';
    end
  catch ME
    error_messages{end+1} = 'Cannot read boundary shapefile.';
  end
  fprintf('%s\n',datestr(now,13));
  
  % Error Files
  % ======================================================
  fprintf('\tError Files ... ');
  try
    error_csv_fn = strcat(basepath,base_fn,'\Errors\',base_fn,'_Crossovers.csv');
    error_mat_fn = strcat(basepath,base_fn,'\Errors\',base_fn,'_Crossovers.mat');
    error_txt_fn = strcat(basepath,base_fn,'\Errors\',base_fn,'_Crossovers_Errors.txt');
    
    good_error_mat_fields = {'dist','cross_angle','Latitude','Longitude','Frame_ID',...
      'Thickness','UTC_time','absError'};
    good_error_csv_fields = {'LATA','LONA','THICKA','FRAMEA','TIMEA','LATB','LONB',...
      'THICKB','FRAMEB','TIMEB','ABSDIFF','ANGLE','DISTANCE'};
    good_error_txt_fields = {'RMS','RMeS'};
    
    % Check MAT File for correct fields
    error_mat_data = load(error_mat_fn); error_mat_fields = fieldnames(error_mat_data);
    if length(error_mat_fields) == length(good_error_mat_fields)
      for error_idx = 1:length(good_error_mat_fields)
        if ~strcmp(good_error_mat_fields{error_idx},error_mat_fields)
          error_messages{end+1} = sprintf('Field [%s] missing from error MAT file.',good_error_mat_fields{error_idx});
        end
      end
    else
      error_messages{end+1} = 'Errors MAT file has incorrect # of fields.';
    end
    % Check CSV File for correct fields
    error_csv_id = fopen(error_csv_fn);
    error_csv_header = textscan(error_csv_id,'%s%s%s%s%s%s%s%s%s%s%s%s%s',1);
    error_csv_header = char(error_csv_header{1});
    error_csv_fields = regexp(error_csv_header,',','split');
    if length(error_csv_fields) == length(good_error_csv_fields)
      for error_idx = 1:length(good_error_csv_fields)
        if ~strcmp(good_error_csv_fields{error_idx},error_csv_fields)
          error_messages{end+1} = sprintf('field [%s] missing from error csv file.',good_error_csv_fields{error_idx});
        end
      end
    else
      error_messages{end+1} = 'Errors csv file has incorrect # of fields.';
    end
    
    % Check that errors CSV and MAT have same number of values.
    error_mat_num_vals = length(error_mat_data.Latitude);
    error_csv_num_vals = 0;
    while(fgets(error_csv_id) ~= -1)
      error_csv_num_vals = error_csv_num_vals+1;
    end
    fclose(error_csv_id);
    if error_mat_num_vals < error_csv_num_vals - 1 % Account for CSV header
      error_messages{end+1} = 'Errors csv file has more rows than MAT file.';
    elseif error_mat_num_vals > error_csv_num_vals - 1
      error_messages{end+1} = 'Errors MAT file has more rows than csv file.';
    end
    
    % Check the errors TXT file
    error_txt_id = fopen(error_txt_fn);
    error_txt_header = textscan(error_txt_id,'%s%s',1);
    error_txt_header = char(error_txt_header{1});
    error_txt_fields = regexp(error_txt_header,',','split');
    if length(error_txt_fields) == length(good_error_txt_fields)
      for error_idx = 1:length(good_error_txt_fields)
        if ~strcmp(good_error_txt_fields{error_idx},error_txt_fields)
          error_messages{end+1} = sprintf('field [%s] missing from error txt file.',good_error_txt_fields{error_idx});
        end
      end
    else
      error_messages{end+1} = 'Errors txt file has incorrect # of fields.';
    end
    fclose(error_txt_id);
  catch ME
    error_messages{end+1} = 'Error Reading Errors Data.';
  end
  fprintf('%s\n',datestr(now,13));
  
  % Flightline Files
  % ======================================================
  fprintf('\tFlightlines Files ... ');
  try
    fline_shp_fn = strcat(basepath,base_fn,'\flightlines\',full_fline_fns(1:length(full_fline_fns)-4),'.shp');
    fline_txt_fn = strcat(basepath,base_fn,'\flightlines\',full_fline_fns(1:length(full_fline_fns)-4),'.txt');
    fline_shp = shapeinfo(fline_shp_fn);
    fl_txt_id = fopen(fline_txt_fn);
    fl_txt_header = textscan(fl_txt_id,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s',1);
    good_fline_fields = {'A_SURF','A_BED','THICK','DATATYPE','SEASON','LAT','LON',...
      'UTCTime','ELEVATION','YYYYMMDD','SEGMENT','FRAME','SURFACE','BOTTOM','QUALITY',};
    
    % Check that shapefile has correct fields
    if length(fline_shp.Attributes) == length(good_fline_fields)
      fline_shp_fields = {};
      for flshp_idx = 1:length(fline_shp.Attributes)
        fline_shp_fields{end+1} = fline_shp.Attributes(flshp_idx).Name;
      end
      for fl_idx = 1:length(good_fline_fields)
        if ~strcmp(good_fline_fields{fl_idx},fline_shp_fields);
          error_messages{end+1} = sprintf('Field [%s] missing from flightline shapefile.',good_fline_fields{fl_idx});
        end
      end
    else
      error_messages{end+1} = 'Flightline shapefile has incorrect # of fields.';
    end
    
    % Check that text file has correct fields
    fl_txt_header = char(fl_txt_header{1});
    fline_txt_fields = regexp(fl_txt_header,',','split');
    if length(fline_txt_fields) == length(good_fline_fields)
      for fl_idx = 1:length(fline_txt_fields)
        if ~strcmp(good_fline_fields{fl_idx},fline_txt_fields);
          error_messages{end+1} = sprintf('Field [%s] missing from flightline text file.',good_fline_fields{fl_idx});
        end
      end
    else
      error_messages{end+1} = 'Flightline text file has incorrect # of fields.';
    end
    
    % Check that the shapefile and text file have the same number of values
    fl_shp_num_vals = fline_shp.NumFeatures;
    fl_txt_num_vals = 0;
    while (fgets(fl_txt_id) ~= -1),
      fl_txt_num_vals = fl_txt_num_vals+1;
    end
    fclose(fl_txt_id);
    if fl_shp_num_vals < fl_txt_num_vals - 1 % Account for TXT header
      error_messages{end+1} = 'Flightline text file has more rows than shapefile.';
    elseif fl_shp_num_vals > fl_txt_num_vals - 1
      error_messages{end+1} = 'Flightline shapefile has more rows than text file.';
    end
  catch ME
    error_messages{end+1} = 'Error Reading Flightline Data.';
  end
  fprintf('%s\n',datestr(now,13));
  
  % Grid Files
  % ======================================================
  fprintf('\tGrid Files ... ');
  try
    grid_surface_fn = strcat(basepath,base_fn,'\grids\',lower(base_fn),'_surface.txt');
    grid_thickness_fn = strcat(basepath,base_fn,'\grids\',lower(base_fn),'_thickness.txt');
    grid_bottom_fn = strcat(basepath,base_fn,'\grids\',lower(base_fn),'_bottom.txt');
    grid_xyz_fn = strcat(basepath,base_fn,'\grids\',base_fn,'_XYZGrid.txt');
    % Try to read ASCII data
    try
      grid_surface_data = arcgridread(grid_surface_fn);
    catch ME
      error_messages{end+1} = 'ASCII Grid Surface did not read correctly.';
    end
    try
      grid_thickness_data = arcgridread(grid_thickness_fn);
    catch ME
      error_messages{end+1} = 'ASCII Grid Thickness did not read correctly.';
    end
    try
      grid_bottom_data = arcgridread(grid_bottom_fn);
    catch ME
      error_messages{end+1} = 'ASCII Grid Bottom did not read correctly.';
    end
    
    % Check grids (Surface - Bottom = Thickness)
    grid_thick_diff = grid_thickness_data - (grid_surface_data - grid_bottom_data);
    if max(max(abs(grid_thick_diff))) > 0.1
      error_messages{end+1} = 'Grids: Surface - Bottom does not equal Thickness.';
    end
    
    % Try to read XYZGrid TXT file
    grid_xyz_id = fopen(grid_xyz_fn);
    grid_xyz_header = textscan(grid_xyz_id,'%s%s%s%s%s',1);
    good_grid_xyz_fields = {'LAT','LON','THICK','SURFACE','BOTTOM'};
    grid_xyz_header = char(grid_xyz_header{1});
    grid_xyz_fields = regexp(grid_xyz_header,',','split');
    if length(grid_xyz_fields) == length(good_grid_xyz_fields)
      for grid_idx = 1:length(grid_xyz_fields)
        if ~strcmp(good_grid_xyz_fields{grid_idx},grid_xyz_fields);
          error_messages{end+1} = sprintf('Field [%s] missing from grids XYZ file.',good_grid_xyz_fields{grid_idx});
        end
      end
    else
      error_messages{end+1} = 'Grid XYZ file has incorrect # of fields.';
    end
    
    % Confirm XYZ Grid has same number of values as grids.
    [r c] = size(grid_surface_data); grid_ascii_num_vals = r*c;
    grid_xyz_num_vals = 0;
    while(fgets(grid_xyz_id) ~= -1)
      grid_xyz_num_vals = grid_xyz_num_vals+1;
    end
    fclose(grid_xyz_id);
    if grid_ascii_num_vals < grid_xyz_num_vals -1 % Account for TXT header.
      error_messages{end+1} = 'Grid XYZ file has more rows than ASCII.';
    elseif grid_ascii_num_vals > grid_xyz_num_vals -1
      error_messages{end+1} = 'Grid ASCII file has more rows than XYZ file.';
    end
  catch ME
    error_messages{end+1} = 'Error Reading Grids Data.';
  end
  fprintf('%s\n',datestr(now,13));
  
  % README Files
  % ======================================================
  fprintf('\tREADME ... ');
  readme_fn = strcat(basepath,base_fn,'\readme\',base_fn,'_README.txt');
  try
    readme_id = fopen(readme_fn);
    if readme_id == -1
      error_messages{end+1} = 'Failed to open README file.';
    end
    fclose(readme_id);
  catch ME
    error_messages{end+1} = 'Error Reading README File.';
  end
  fprintf('%s\n',datestr(now,13));
  
  % Kill Block
  if ~isempty(error_messages)
    fprintf('ERROR MESSAGES:\n');
    for error_idx = 1:length(error_messages)
      fprintf('(%d) %s \n',error_idx,error_messages{error_idx});
    end
    error('Errors Found. See Messages Above.');
  end
  
  %% Visual Check of Preview_Images
  fprintf('Checking Preview Images:\n');
  fprintf('------------------------------\n')
  
  bottom_preview_fn = strcat(basepath,base_fn,'\preview_images\',base_fn,'_Bottom_Preview.png');
  surface_preview_fn = strcat(basepath,base_fn,'\preview_images\',base_fn,'_Surface_Preview.png');
  thickness_preview_fn = strcat(basepath,base_fn,'\preview_images\',base_fn,'_Thickness_Preview.png');
  flightline_preview_fn = strcat(basepath,base_fn,'\preview_images\',base_fn,'_Flightline_Preview.png');
  errors_preview_fn = strcat(basepath,base_fn,'\errors\',base_fn,'_Crossovers_Preview.png');
  
  fprintf('\tThe preview images will be loaded and displayed one at a time.\n')
  fprintf('\tYou will neeed to verify the following features:\n');
  fprintf('\t\t 1) Correct Title (Name and Spelling).\n');
  fprintf('\t\t 2) Scale Bar Present (In Kilometers).\n');
  fprintf('\t\t 3) CReSIS and NSF Logos are Present.\n');
  fprintf('\t\t 4) Legend Present (Blue:Low,Red:High).\n');
  fprintf('\t\t 5) Legend Title Correct (Name and Spelling).\n');
  fprintf('\t\t 6) Legend Values (Ex. 3456 not 3455.84).\n');
  fprintf('\t\t 7) Location Inset Map Correct (Red Outline).\n');
  fprintf('\t\t Note: Errors Preview will just be a simple MATLAB plot,\n');
  fprintf('\t\t       simply confirm it loads and displays correctly.\n\n');
  start_previews = input('ARE YOU READY TO BEGIN VISUAL CHECKS? (Y/N): ','s');
  while ~strcmpi(start_previews,'y');
    if strcmpi(start_previews,'n')
      fprintf('WHY NOT ?! Enter Y to begin...\n');
      start_previews = input('ARE YOU READY TO BEGIN VISUAL CHECKS? (Y/N): ','s');
    end
    fprintf('INVALID INPUT. ENTER Y for yes or N for NO.\n');
    start_previews = input('ARE YOU READY TO BEGIN VISUAL CHECKS? (Y/N): ','s');
  end
  
  % Check Bottom Preview
  fprintf('Checking Bottom Preview Image:\n');
  try
    bottom_image = imread(bottom_preview_fn);
    image(bottom_image); axis tight; axis off;
    bottom_response = input('Is the image correct? (Y/N): ','s');
    while ~strcmpi(bottom_response,{'y','n'});
      fprintf('INVALID INPUT. ENTER Y for yes or N for NO.\n');
      bottom_response = input('Is the image correct? (Y/N): ','s');
    end
    if strcmpi(bottom_response,'n')
      error_messages{end+1} = 'Bottom Preview_Image contains errors.';
      fprintf('\n');
    else
      fprintf('Bottom Preview_Image Check Complete.\n\n');
      close;
    end
  catch ME
    error('Error loading Bottom Preview_Image.');
  end
  % Check Thickness Preview
  fprintf('Checking Thickness Preview Image:\n');
  try
    thickness_image = imread(thickness_preview_fn);
    image(thickness_image); axis tight; axis off;
    thickness_response = input('Is the image correct? (Y/N): ','s');
    while ~strcmpi(thickness_response,{'y','n'});
      fprintf('INVALID INPUT. ENTER Y for yes or N for NO.\n');
      thickness_response = input('Is the image correct? (Y/N): ','s');
    end
    if strcmpi(thickness_response,'n')
      error_messages{end+1} = 'Thickness Preview_Image contains errors.';
      fprintf('\n');
    else
      fprintf('Thickness Preview_Image Check Complete.\n\n');
      close;
    end
  catch ME
    error('Error loading Thickness Preview_Image.');
  end
  % Check Surface Preview
  fprintf('Checking Surface Preview Image:\n');
  try
    surface_image = imread(surface_preview_fn);
    image(surface_image); axis tight; axis off;
    surface_response = input('Is the image correct? (Y/N): ','s');
    while ~strcmpi(surface_response,{'y','n'});
      fprintf('INVALID INPUT. ENTER Y for yes or N for NO.\n');
      surface_response = input('Is the image correct? (Y/N): ','s');
    end
    if strcmpi(surface_response,'n')
      error_messages{end+1} = 'Surface Preview_Image contains errors.';
      fprintf('\n');
    else
      fprintf('Surface Preview_Image Check Complete.\n\n');
      close;
    end
  catch ME
    error('Error loading Bottom Preview_Image.');
  end
  % Check Flightline Preview
  fprintf('Checking Flightline Preview Image:\n');
  try
    flightline_image = imread(flightline_preview_fn);
    image(flightline_image); axis tight; axis off;
    flightline_response = input('Is the image correct? (Y/N): ','s');
    while ~strcmpi(flightline_response,{'y','n'});
      fprintf('INVALID INPUT. ENTER Y for yes or N for NO.\n');
      flightline_response = input('Is the image correct? (Y/N): ','s');
    end
    if strcmpi(flightline_response,'n')
      error_messages{end+1} = 'Flightline Preview_Image contains errors.';
      fprintf('\n');
    else
      fprintf('Flightline Preview_Image Check Complete.\n\n');
      close;
    end
  catch ME
    error('Error loading Flightline Preview_Image.');
  end
  % Check Errors Preview
  fprintf('Checking Error Preview Image:\n');
  try
    error_image = imread(errors_preview_fn);
    image(error_image); axis tight; axis off;
    error_response = input('Is the image correct? (Y/N): ','s');
    while ~strcmpi(error_response,{'y','n'});
      fprintf('INVALID INPUT. ENTER Y for yes or N for NO.\n');
      error_response = input('Is the image correct? (Y/N): ','s');
    end
    if strcmpi(error_response,'n')
      error_messages{end+1} = 'Error Preview_Image contains errors.';
      fprintf('\n');
    else
      fprintf('Error Preview_Image Check Complete.\n\n');
      close;
    end
  catch ME
    error('Error loading Errors Preview_Image.');
  end
  
  %% Clean Up and Report Errors
  
  if exist(zip_fn(1:length(zip_fn)-4),'dir')
    rmdir(zip_fn(1:length(zip_fn)-4),'s');
  end
  % Kill Block
  if ~isempty(error_messages)
    fprintf('ERROR MESSAGES:\n');
    for error_idx = 1:length(error_messages)
      fprintf('(%d) %s \n',error_idx,error_messages{error_idx});
    end
    fprintf('=================================\n');
    fprintf('Errors Found. See Messages Above.\n');
    fprintf('=================================\n\n');
  else
    fprintf('=================================\n');
    fprintf('No errors found. Script Complete.\n');
    fprintf('=================================\n\n');
  end
  
catch ME
  % If error at any time delete unzipped folder.
  fprintf('\n========== ERROR IN SCRIPT ==========\n');
  fid_s = fopen('all');
  if ~isempty(fid_s)
    fclose('all');
  end
  rmdir(zip_fn(1:length(zip_fn)-4),'s');
  rethrow(ME);
end
end