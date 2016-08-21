function [success,surfTimes] = tomo_collate_task(param)

  surfTimes = [];
  
  fn_dir = ct_filename_out(param,param.surf_extract.out_dir);
  frm = param.frm;
  
  global gRadar
  gRadar.tmp_path = param.tmp_path;
  
  for img=1:3
    fn = fullfile(fn_dir,sprintf('Data_img_%02.0f_%s_%03.0f.mat',img, ...
      param.day_seg,frm));
    mdata{img} = load(fn);
    mdata{img}.frm = frm;
  end
  
    surf_extract.add_layers_flag = 1;
    surf_extract.ice_twtt_flag = 1;
    surf_extract.extract_flag = 1;
    surf_extract.theta_calibrated = 1;
    surf_extract.surf_flag = 1;
  
  if param.surf_extract.add_layers_flag
    mdata = tomo.data_loader_prep(param,mdata);
  end
  if param.surf_extract.ice_twtt_flag
    mdata = tomo.DEM_alignment(param,mdata);
  end
  if param.surf_extract.extract_flag
    mdata_combined = tomo.surface_extractor(param,mdata);
    if param.surf_extract.surf_flag
      tomo.create_surf(param,mdata_combined);
    end
  end
  
  success = true;
  
  return;
    