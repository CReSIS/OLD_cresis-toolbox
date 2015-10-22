function run_picker_training(user_dir,varargin)
% Runs the Training Picker.
%
% run_picker_training(user_dir,flag);
%
% user_dir: A path to a folder on the users local directory.
% flag: Optional Integer of "1" will copy and load the correct picks.
%
% Both UNIX and WINDOWS paths are supported.
%
% The tool will copy over the data and then run the picker in that new
% directory. This is used for training someone how to pick. See the CReSIS
% Wiki Page for more details.
%
% WARNING: If you are re-loading your previously copied picks the tool will
% ask if you want to overwrite the directory. Saying Y(Yes) will delete ALL
% of your picks and return to RAW data form. Say N(No) to preserve your
% picks.
%
% Example:
% -------------------
% LOAD THE RAW DATA FOR PICKING
% user_dir = 'C:\Users\SomeUser\SomeFolder\SomePickingFolder\';
% run_picker_training(user_dir);
%
% LOAD THE CORRECT PICKS
% user_dir = 'C:\Users\SomeUser\SomeFolder\SomePickingFolder\';
% run_picker_training(user_dir,1);
%
% Author: Kyle Purdon
%
% See also picker.m

%%
fprintf('=====================\n')
fprintf('   PICKER TRAINING   \n');
fprintf('=====================\n')

%% Check for Optional Input (Load Correct Picks Yes/No)
load_c = false;
if nargin > 1
  load_c = true;
end

%% Check if user_dir exists, if not create it.
if ~exist(user_dir,'dir')
  mkdir(user_dir);
end

%% Check if user_dir is empty, print warning if it's not.
overwrite = 'y';
if ~isempty(get_filenames(user_dir,'','','','recursive')) && ~load_c
  fprintf('WARNING: "user_dir" is NOT empty.\n')
  tmp = false;
  while tmp==false
    overwrite = input('Do you want to overwrite layerData?(Y/N): ','s');
    if strcmpi(overwrite,'n') || strcmpi(overwrite,'y')
      tmp = true;
    end
  end
end

%% Make Path/s to layer_data
if ispc
  layer_path = 'Z:\mdce\training\picker_training';
  c_layer_path = 'Z:\mdce\training\correct_picks';
else
  layer_path = '/cresis/scratch2/mdce/training/picker_training';
  c_layer_path = '/cresis/scratch2/mdce/training/correct_picks';
end

%% Remove end filesep on user_dir if it exists
if strcmp(user_dir(end),filesep)
  user_dir(end) = '';
end

%% Copy The LayerData
if load_c % Copy the correct picks (Always Overwrite)
  user_dir = strcat(user_dir,filesep,'correct');
  mkdir(user_dir);
  tic;
  fprintf('Copying Correct Picks (May take a few minutes)');
  if ispc
    % Windows system call to robocopy
    str = ['!robocopy "',c_layer_path,'" "',user_dir,'" ','/E /NC /NS /NDL /NFL /NJH /NJS /NP'];
    eval(str);
  else
    % Linux system call to cp
    str = ['!cp -r ',c_layer_path,'/* ',user_dir];
    eval(str);
  end
  fprintf('  Done (%.1f sec)\n', toc);
else % Copy data from mdce to user_dir
  if strcmpi(overwrite,'y')
    tic;
    rmdir(user_dir,'s');mkdir(user_dir);
    fprintf('Copying Layer Data (May take a few minutes)');
    if ispc
      % Windows system call to robocopy
      str = ['!robocopy "',layer_path,'" "',user_dir,'" ','/E /NC /NS /NDL /NFL /NJH /NJS /NP'];
      eval(str);
    else
      % Linux system call to cp
      str = ['!cp -r ',layer_path,'/* ',user_dir];
      eval(str);
    end
    fprintf('  Done (%.1f sec)\n', toc);
  end
end

%% Run the Training Picker
clear param;
param.radar_name = 'training';
param.season_name = 'picker_training';

post_dir = '';

source_data = {};
source_data{end+1} = fullfile(user_dir,post_dir,'CSARP_csarp-combined');
source_data{end+1} = fullfile(user_dir,post_dir,'CSARP_standard');
source_data{end+1} = fullfile(user_dir,post_dir,'CSARP_mvdr');
source_data{end+1} = fullfile(user_dir,post_dir,'CSARP_qlook');

layer_data = fullfile(user_dir,post_dir,'CSARP_layerData');

geotiff_fns = {};
geotiff_fns{end+1} = fullfile(ct_filename_gis(param,'greenland'),'Landsat-7','Greenland_natural_250m.tif');
geotiff_fns{end+1} = fullfile(ct_filename_gis(param,'canada'),'Landsat-7','Canada_250m.tif');

param.fast_load.en = false;
param.fast_load.recreate = false;
param.fast_load.tmp_file = '';

param.landmarks = [];

picker(source_data,layer_data,geotiff_fns,param);

return;