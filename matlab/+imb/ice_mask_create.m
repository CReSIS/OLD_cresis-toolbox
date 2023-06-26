% script ice_mask_create.m
%
% Script to create ice mask files for the imb.slice_browser's slicetool_icemask
%
% Author: John Paden

ice_mask_tif_fn = ct_filename_gis(param,fullfile('greenland','IceMask','GimpIceMask_90m_v1.1.tif'));

ice_mask_bin_fn = ct_filename_gis(param,fullfile('greenland','IceMask','GimpIceMask_90m_v1.1.bin'));

ice_mask.proj = geotiffinfo(ice_mask_tif_fn);
[ice_mask.mask, ice_mask.R, ~] = geotiffread(ice_mask_tif_fn);
ice_mask.X = ice_mask.R(3,1) + ice_mask.R(2,1)*(1:size(ice_mask.mask,2));
ice_mask.Y = ice_mask.R(3,2) + ice_mask.R(1,2)*(1:size(ice_mask.mask,1));



[ice_mask_fn_dir ice_mask_fn_name] = fileparts(ice_mask_bin_fn);
ice_mask_mat_fn = fullfile(ice_mask_fn_dir,[ice_mask_fn_name '.mat']);
save(ice_mask_mat_fn,'-v7.3','-struct','ice_mask','R','X','Y','proj','mask');

[fid,msg] = fopen(ice_mask_bin_fn,'w');
if fid < 1
  fprintf('Could not open file %s\n', ice_mask_bin_fn);
  error(msg);
end
fwrite(fid,uint8(ice_mask.mask),'uint8');
fclose(fid);
