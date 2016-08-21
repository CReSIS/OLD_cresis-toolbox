function create_surf(param,mdata)

  surf = [];

  if isfield(mdata,'twtt')
    surf = [];
    surf(1).x = repmat((1:64).',[1 size(mdata.twtt,2)]);
    surf(1).y = interp1(mdata.Time,1:length(mdata.Time),mdata.twtt);
    surf(1).name = 'surface';
  end
  
  if isfield(mdata,'bottom_surface')
    surf(2).x =  repmat((1:64).',[1 size(mdata.twtt,2)]);
    surf(2).y = mdata.bottom_surface.';
    surf(2).name = 'bottom';
  end
  
%   surf(3).x = repmat((1:64).',[1 size(mdata.twtt,2)]);
%   surf(3).y = NaN * zeros(size(surf(1).y));
%   surf(3).y(33,:) = interp1(mdata.Time,1:length(mdata.Time),mdata.Bottom);
%   surf(3).plot_name_values = {'color','magenta','marker','+'};
%   surf(3).name = 'Bcontrol';
%   surf(3).surf_surf = 1;
%   surf(3).active_surf = 2;
%   surf(3).control_surf = 3;
%   surf(3).mask_surf = 4;

  if isfield(mdata,'ice_mask')
    surf(3).x = repmat((1:64).',[1 size(mdata.twtt,2)]);
    surf(3).y = mdata.ice_mask;
    surf(3).name = 'ice mask';
  end
    
  if ~isempty(surf)
    surf_fn = sprintf('Data_%s_%03.0f.mat', ...
      param.day_seg,mdata.frm);
    surf_dir = ct_filename_out(param,[],'CSARP_surfData');
    if ~isdir(surf_dir)
      mkdir(surf_dir);
    end
    surf_file = fullfile(surf_dir,surf_fn);
    
    save(surf_file,'surf','-v7.3');    
  end
end