% script ice_mask_update.m
%
% This script updates the .mat ice mask file from the .bin ice mask file.
% The .mat file is what tomo.create_surfdata.m uses and .bin is what
% imb.run_slice_browser.m uses in slicetool_icemask.m. Since the .bin file
% is updated by users inside the slice browser, the .mat file needs to be
% updated occasionally so that if create_surfdata is run again, it will
% have these updates included.
%
% https://wiki.cresis.ku.edu/cresis/Automatic_3D_Surface_Tracking#add_icemask_surfacedem
%
% Author: John Paden

%% Set the path to the .bin ice mask file
ice_mask_fn = ct_filename_gis([],fullfile('canada','ice_mask','03_rgi50_ArcticCanadaNorth','03_rgi50_ArcticCanadaNorth.bin'));

%% Load the .mat file's XY variables since these are needed to load the .bin file
[ice_mask_fn_dir ice_mask_fn_name] = fileparts(ice_mask_fn);
ice_mask_mat_fn = fullfile(ice_mask_fn_dir,[ice_mask_fn_name '.mat']);
fprintf('Loading .mat ice mask from %s\n', ice_mask_mat_fn);
ice_mask = load(ice_mask_mat_fn,'X','Y');

%% Load the .bin ice mask file (simple flat binary file with a single matrix)
fprintf('Reading .bin ice mask from %s\n', ice_mask_fn);
[fid,msg] = fopen(ice_mask_fn,'r');
if fid < 1
  fprintf('Could not open file %s\n', ice_mask_bin_fn);
  error(msg);
end
ice_mask.mask = logical(fread(fid,[length(ice_mask.Y),length(ice_mask.X)],'uint8'));
fclose(fid);

%% Write the .bin's ice mask back into the .mat file.
fprintf('Writing .bin ice mask to %s\n', ice_mask_mat_fn);
save(ice_mask_mat_fn,'-append','-struct','ice_mask','mask');
