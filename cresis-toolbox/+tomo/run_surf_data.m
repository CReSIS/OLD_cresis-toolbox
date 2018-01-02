filename = 'X:\ct_data\rds\2014_Greenland_P3\CSARP_shane_music_surfData\20140401_03\Data_20140401_03_006.mat';

% constructor with no argument
s = surf_data

% with constructor
my_surf = surf_data(filename);
my_surf


% invalid filename
my_surf2 = surf_data(1); % error; 1 is a number not a string

%% load_surf
my_surf.load_surf(filename);

% invalid filename
my_surf.load_surf(1);

%% insert

% create a test surf structure
test_surf = my_surf.surf(2);
test_surf.name = 'my_test_surf';
my_surf.insert_surf(test_surf);
my_surf

% insert into empty surf_data
s.insert_surf(test_surf);
s

% insert another one
test_surf = my_surf.surf(3);
s.insert_surf(test_surf);
% insert another five
for index = 1:5
  s.insert_surf(my_surf.surf(index+3));
end
s

%% insert with size mismatch in x or y
% create a test_surf structure where (size of x) != (size of y)
%my_surf = surf_data;
test_surf.x = test_surf.x(:,1:300);
test_surf.name = 'my_test_surf_with_mismatch_xy';
my_surf.insert_surf(test_surf); % error

%% insert with size mismatch with previously inserted x or y
% create a test_surf structure where (size of x) == (size of y) but does
% not match with surface that are inserted before

test_surf = my_surf.surf(1);
test_surf.x = test_surf.x(:,1:300);
test_surf.y = test_surf.y(:,1:300);
s2 = surf_data;
s2.insert_surf(test_surf);

my_surf.insert_surf(test_surf); % error (size(x) = [64 3332] for all x in surf struct in my_surf)

test_surf = my_surf.surf(2);
s2.insert_surf(test_surf); % error (size(x) = [64 300] for all surf struct in s2)

% we don't test y because dimension of x and y must match before this test

%% other insert problems
% insert a invalid struct
empty_struct = struct();
my_surf.insert_surf(empty_struct); % error

% insert surf struct with invalid name
surf_with_invalid_name = my_surf.surf(1); 
surf_with_invalid_name.name = 123;
my_surf.insert_surf(surf_with_invalid_name); %  error

% insert a duplicate
duplicate_bottom_surf = my_surf.surf(2);
my_surf.insert_surf(duplicate_bottom_surf);  % error

%% get surf
my_surf.get_surf('bottom') % get the bottom surface
my_surf.get_surf('123') % surface does not exist; return []
my_surf.get_surf([1 2 3]) % wrong input type



%% set name
my_surf.set_name('my_test_surf', 'my_test_surf2');
my_surf.get_surf('my_test_surf2')

my_surf.set_name('123', '1234'); % surf not exist
my_surf.set_name([1 2 3], [1 2 3 4]); % wrong input type
my_surf.set_name([1 2 3], '123'); % wrong input type
my_surf.set_name('123', 12 ); % wrong input type

%% set surf
my_surf.get_surf('my_test_surf2') % before
modified_test_surf = my_surf.get_surf('my_test_surf2');
modified_test_surf.surf_layer = [];
modified_test_surf.active_layer = [];
modified_test_surf.mask_layer = [];
modified_test_surf.control_layer = [];
modified_test_surf.quality_layer = [];

my_surf.set_surf(modified_test_surf);
my_surf.get_surf('my_test_surf2') % after; check result

% set a surf that is not valid
my_surf.set_surf(struct()); % error
% set a surface where the surface is not in the object
test123_surface = modified_test_surf;
test123_surface.name = 'test123';
my_surf.set_surf(test123_surface); % error

%% remove
my_surf.remove_surf('my_test_surf2');
my_surf

%% remove; see the index changes
% before
v_names = {'name', 'index', 'surf_layer', 'active_layer', 'mask_layer', 'control_layer'};

my_surf_name = {my_surf.surf.name}';
indices = 1:length(my_surf.surf); indices = indices';
my_surf_surf_layer = {my_surf.surf.surf_layer}';
my_surf_active_layer = {my_surf.surf.active_layer}';
my_surf_mask_layer = {my_surf.surf.mask_layer}';
my_surf_control_layer = {my_surf.surf.control_layer}';

table(my_surf_name, indices, my_surf_surf_layer, ...
  my_surf_active_layer, my_surf_mask_layer, ...
  my_surf_control_layer, 'VariableNames', v_names)

%remove -- remove anything in the surf array
my_surf.remove_surf('bottom');

% after
my_surf_name = {my_surf.surf.name}';
indices = 1:length(my_surf.surf)'; indices = indices';
my_surf_surf_layer = {my_surf.surf.surf_layer}';
my_surf_active_layer = {my_surf.surf.active_layer}';
my_surf_mask_layer = {my_surf.surf.mask_layer}';
my_surf_control_layer = {my_surf.surf.control_layer}';

table(my_surf_name, indices, my_surf_surf_layer, ...
  my_surf_active_layer, my_surf_mask_layer, ...
  my_surf_control_layer, 'VariableNames', v_names)


%% remove a surface that is not in the struct array
my_surf.remove_surf('123'); % error

%% save
% save as my_surface.mat
my_surf.save_surf('my_surface'); 






