%% Setup

% run_example_str = 'demo_new_file';
% run_example_str = 'demo_exist_file';
run_example_str = 'demo_update_file';

if strcmpi(run_example_str,'demo_new_file')
  % Demonstrate creation of a new surf file
  param = read_param_xls(ct_filename_param('rds_param_2019_Antarctica_Ground.xls'),'20200107_01');
  mdata = echo_load(param,'music3D_paden',1);
  surf = tomo.surfdata(mdata,'surfData_paden');
end

if strcmpi(run_example_str,'demo_exist_file')
  % Demonstrate loading an existing surf file and manipulating it
  
  param = read_param_xls(ct_filename_param('rds_param_2009_Antarctica_TO.xls'),'20091224_01');
  fn_original = fullfile(ct_filename_out(param,'surfData',''),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
  fn = [tempname '.mat'];
  copyfile(fn_original,fn);
  
  surf = tomo.surfdata(fn);
  
  %% get surf
  my_surfdata.get_surf('bottom') % get the bottom surface
  my_surfdata.get_surf({'top','bottom'}) % get the top and bottom surfaces
  my_surfdata.get_surf([1 2 3])
  try
    my_surfdata.get_surf('123') % surface does not exist;
  catch ME
    warning(ME.getReport);
  end
  try
    my_surfdata.get_surf(20) % invalid index
  catch ME
    warning(ME.getReport);
  end
  try
    my_surfdata.get_surf(0) % invalid index
  catch ME
    warning(ME.getReport);
  end
  
  %% insert_surf
  
  % Create a test surf structure from an existing surface
  test_surf = my_surfdata.get_surf('bottom');
  test_surf.name = 'my_test_surf';
  
  % Insert test surface
  my_surfdata.insert_surf(test_surf);
  my_surfdata
  
  % Insert multiple layers
  for surf_idx = 1:length(my_surfdata.surf)
    test_surf = my_surfdata.surf(surf_idx);
    empty_surfdata.insert_surf(test_surf);
  end
  empty_surfdata
  
  % create a test_surf structure where (size of x) != (size of y)
  test_surf = my_surfdata.get_surf('bottom');
  test_surf.x = test_surf.x(:,1:300);
  test_surf.name = 'test';
  try
    my_surfdata.insert_surf(test_surf);
  catch ME
    warning(ME.getReport);
  end
  
  % create a test_surf which does not match existing surf x and y sizes
  test_surf = my_surfdata.get_surf('bottom');
  test_surf.x = test_surf.x(1:32,1:300);
  test_surf.y = test_surf.y(1:32,1:300);
  test_surf.name = 'test';
  try
    my_surfdata.insert_surf(test_surf);
  catch ME
    warning(ME.getReport);
  end
  
  % insert an invalid struct
  test_surf = struct();
  try
    my_surfdata.insert_surf(test_surf);
  catch ME
    warning(ME.getReport);
  end
  
  % insert surf struct with invalid name
  test_surf = my_surfdata.surf(1);
  test_surf.name = 123;
  try
    my_surfdata.insert_surf(test_surf);
  catch ME
    warning(ME.getReport);
  end
  
  % insert a duplicate
  test_surf = my_surfdata.surf(2);
  try
    my_surfdata.insert_surf(test_surf);
  catch ME
    warning(ME.getReport);
  end
  
  %% set_name
  my_surfdata.set('bottom', 'name', 'bottom2');
  my_surfdata.get_surf('bottom2')
  my_surfdata.set('bottom2', 'name', 'bottom');
  
  try
    my_surfdata.set('123','name','1234'); % surf not exist
  catch ME
    warning(ME.getReport);
  end
  
  %% set_surf and clear_references
  test_surf = my_surfdata.get_surf('bottom')
  test_surf = my_surfdata.clear_references(test_surf)
  
  my_surfdata.set_surf(test_surf);
  my_surfdata.get_surf('bottom') % after; check result
  
  % set a surf that is not valid
  try
    my_surfdata.set_surf(struct());
  catch ME
    warning(ME.getReport);
  end
  
  % set a surface where the surface name does not match an existing surface
  try
    test_surf = my_surfdata.get_surf('bottom');
    test_surf.name = 'other';
    my_surfdata.set_surf(test_surf);
  catch ME
    warning(ME.getReport);
  end
  
  %% remove_surf
  my_surfdata.load_surfdata(fn);
  test_surf = my_surfdata.get_surf('bottom');
  test_surf.name = 'test';
  my_surfdata.insert_surf(test_surf);
  my_surfdata
  my_surfdata.remove_surf('test');
  my_surfdata
  
  % See the index changes
  v_names = {'name', 'index', 'top', 'active', 'mask', 'gt', 'quality'};
  
  my_surfdata_name = {my_surfdata.surf.name}';
  indices = (1:length(my_surfdata.surf))';
  my_surfdata_top = {my_surfdata.surf.top}';
  my_surfdata_active = {my_surfdata.surf.active}';
  my_surfdata_mask = {my_surfdata.surf.mask}';
  my_surfdata_gt = {my_surfdata.surf.gt}';
  my_surfdata_quality = {my_surfdata.surf.quality}';
  
  table(my_surfdata_name, indices, my_surfdata_top, ...
    my_surfdata_active, my_surfdata_mask, ...
    my_surfdata_gt, my_surfdata_quality, ...
    'VariableNames', v_names)
  
  my_surfdata.remove_surf('bottom');
  
  my_surfdata_name = {my_surfdata.surf.name}';
  indices = (1:length(my_surfdata.surf))';
  my_surfdata_top = {my_surfdata.surf.top}';
  my_surfdata_active = {my_surfdata.surf.active}';
  my_surfdata_mask = {my_surfdata.surf.mask}';
  my_surfdata_gt = {my_surfdata.surf.gt}';
  my_surfdata_quality = {my_surfdata.surf.quality}';
  
  table(my_surfdata_name, indices, my_surfdata_top, ...
    my_surfdata_active, my_surfdata_mask, ...
    my_surfdata_gt, my_surfdata_quality, ...
    'VariableNames', v_names)
  
  test_surf.name = 'bottom';
  my_surfdata.insert_surf(test_surf);
  my_surfdata.set({'bottom','ice mask','bottom gt','bottom quality'}, ...
    'top','top','active','bottom','mask','ice mask','gt','bottom gt','quality','bottom quality');
  
  my_surfdata_name = {my_surfdata.surf.name}';
  indices = (1:length(my_surfdata.surf))';
  my_surfdata_top = {my_surfdata.surf.top}';
  my_surfdata_active = {my_surfdata.surf.active}';
  my_surfdata_mask = {my_surfdata.surf.mask}';
  my_surfdata_gt = {my_surfdata.surf.gt}';
  my_surfdata_quality = {my_surfdata.surf.quality}';
  
  table(my_surfdata_name, indices, my_surfdata_top, ...
    my_surfdata_active, my_surfdata_mask, ...
    my_surfdata_gt, my_surfdata_quality, ...
    'VariableNames', v_names)
  
  try;
    my_surfdata.remove_surf('123'); % error
  catch ME
    warning(ME.getReport)
  end
  
  %% save_surfdata
  my_surfdata.save_surfdata(fn);
end

% Demonstrate updating a file to the newest format
if strcmpi(run_example_str,'demo_update_file')
  
  param = read_param_xls(ct_filename_param('rds_param_2009_Antarctica_TO.xls'),'20091224_01');
  fn_original = fullfile(ct_filename_out(param,'surfData',''),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
  fn_new = [tempname '.mat'];
  
  echogram_fn = fullfile(ct_filename_out(param,'music3D',''),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
  
  tomo.surfdata.update_file(fn_original,fn_new,echogram_fn);
  
  surf = tomo.surfdata(fn_new);

end
