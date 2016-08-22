function create_surf(param,mdata)

  surf = [];

  if isfield(mdata,'twtt') && isfield(mdata,'Time')
    surf(end+1).x = repmat((1:64).',[1 size(mdata.twtt,2)]);
    surf(end).y = interp1(mdata.Time,1:length(mdata.Time),mdata.twtt);
    surf(end).plot_name_values = {'color','black','marker','x'};
    surf(end).name = 'surface';
    surf(end).surf_layer = 1;
    surf(end).active_layer = 1;
    surf(end).mask_layer = 3;
    surf(end).control_layer = 1;
  end
  
  if isfield(mdata,'bottom_surface') && isfield(mdata,'twtt')
    surf(end+1).x =  repmat((1:64).',[1 size(mdata.twtt,2)]);
    surf(end).y = mdata.bottom_surface.';
    surf(end).plot_name_values = {'color','blue','marker','^'};
    surf(end).name = 'bottom';
    surf(end).surf_layer = 1;
    surf(end).active_layer = 2;
    surf(end).mask_layer = 3;
    surf(end).control_layer = 4;
  end

  if isfield(mdata,'ice_mask') && isfield(mdata,'twtt')
    surf(end+1).x = repmat((1:64).',[1 size(mdata.twtt,2)]);
    surf(end).y = mdata.ice_mask;
    surf(end).plot_name_values = {'color','white','marker','x'};
    surf(end).name = 'ice mask';
    surf(end).surf_layer = 1;
    surf(end).active_layer = 2;
    surf(end).mask_layer = 3;
    surf(end).control_layer = 4;
  end
  
  if isfield(mdata,'twtt') && isfield(mdata,'Bottom')
    surf(end+1).x = repmat((1:64).',[1 size(mdata.twtt,2)]);
    surf(end).y = NaN * zeros(size(surf(1).y));
    surf(end).y(33,:) = interp1(mdata.Time,1:length(mdata.Time),mdata.Bottom);
    surf(end).plot_name_values = {'color','magenta','marker','+'};
    surf(end).name = 'Bcontrol';
    surf(end).surf_layer = 1;
    surf(end).active_layer = 2;
    surf(end).mask_layer = 3;
    surf(end).control_layer = 4;
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