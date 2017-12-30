%% Setup

% Create filename which will be used to demonstrate surfdata class
param = struct('radar_name','mcords','season_name','2009_Antarctica_TO','day_seg','20091224_01');
frm = 16;
fn_original = fullfile(ct_filename_out(param,'surfData_v2',''),sprintf('Data_%s_%03d.mat',param.day_seg,frm));

fn = [tempname '.mat'];
copyfile(fn_original,fn);


%% Constructor

% No arguments
empty_surfdata = tomo.surfdata()

% Filename argument
my_surfdata = tomo.surfdata(fn);
my_surfdata

% Invalid filename example
try
  my_surfdata = tomo.surfdata(1); % error; 1 is a number not a string
catch ME
  warning(ME.getReport);
end

%% load_surfdata
my_surfdata.load_surfdata(fn);

% Invalid filename example
try
  my_surfdata.load_surf(1);
catch ME
  warning(ME.getReport);
end

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

