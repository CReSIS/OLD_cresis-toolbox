% script hd5_example
%
% Example of how to load the contents of a struct array.
% Specific example (see bottom of this file) is for a layer file where
% names of all layers are to be loaded:
%   layer(idx).name

if 0
  %%
  layer = [];
  layer.name = '1234';
  % layer(2).name = '3';
  delete('test.mat');
  save('test.mat','-v7.3','layer')
  a=hdf5info('test.mat');
  a.GroupHierarchy.Groups(1).Datasets
  a.GroupHierarchy.Groups(1).Datasets.Datatype
elseif 1
  %%
  num_layers = 5;
  layer = [];
  for idx=1:length(layer)
    layer(idx).group_name = 'standard';
  end
  for idx=1:length(layer)
    layer(idx).quality = 1;
  end
  layer(1).name = '012';
  layer(2).name = '345';
  layer(3).name = '678';
  layer(4).name = 'ABC';
  layer(5).name = 'DEF';
  for idx=1:length(layer)
    layer(idx).age = idx;
  end
  for idx=1:length(layer)
    layer(idx).description = 'afsfdafsd';
  end
  param = [];
  param.day_seg = '12345678_90'
  delete('test.mat');
  save('test.mat','-v7.3','param','layer')
  if 0
    a=hdf5info('test.mat');
    a.GroupHierarchy.Groups(2)
    a.GroupHierarchy.Groups(2).Attributes(1).Value.Data % struct
    a.GroupHierarchy.Groups(2).Attributes(2).Value(1).Data(1).Data % name
    a.GroupHierarchy.Groups(2).Attributes(2).Value(2).Data(1).Data % age
    for idx0 = 1:length(a.GroupHierarchy.Groups(2).Datasets)
      a.GroupHierarchy.Groups(2).Datasets(idx0)
      a.GroupHierarchy.Groups(2).Datasets(idx0).Datatype
    end
    for idx = 1:length(a.GroupHierarchy.Groups(1).Datasets)
      a.GroupHierarchy.Groups(1).Datasets(idx)
      for aidx = 1:length(a.GroupHierarchy.Groups(1).Datasets(idx).Attributes)
        a.GroupHierarchy.Groups(1).Datasets(idx).Attributes(aidx)
        if isa(a.GroupHierarchy.Groups(1).Datasets(idx).Attributes(aidx).Value,'hdf5.h5string')
          a.GroupHierarchy.Groups(1).Datasets(idx).Attributes(aidx).Value.Data
        end
      end
    end
  end
end

% % >>> f[struArray['name'][0, 0]].value.tobytes()[::2].decode()
% /layer/name

filename = 'test.mat';
dataset = '/layer/name';
fid = H5F.open(filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
dset_id = H5D.open(fid, dataset);
dspace_id = H5D.get_space(dset_id);
type_id = H5D.get_type(dset_id);
% class_id
% H5ML.get_constant_value('H5T_INTEGER') %0
% H5ML.get_constant_value('H5T_FLOAT') %1
% H5ML.get_constant_value('H5T_STRING') %3
% H5ML.get_constant_value('H5T_BITFIELD') %4
% H5ML.get_constant_value('H5T_OPAQUE') %5
% H5ML.get_constant_value('H5T_COMPOUND') %6
% H5ML.get_constant_value('H5T_REFERENCE') %7
% H5ML.get_constant_value('H5T_ENUM') %8
% H5ML.get_constant_value('H5T_VLEN')
% H5ML.get_constant_value('H5T_ARRAY')
% H5ML.get_constant_value('H5T_COMPLEX') %?

% Contents
if (H5ML.compare_values(H5T.get_class(H5D.get_type(dset_id)), ...
    H5ML.get_constant_value('H5T_REFERENCE')))
  
  % Read a "virtual" padding tile that is stored elsewhere in the file.
  plist = 'H5P_DEFAULT';
  % space = 'H5S_ALL'; % dspace_id
  %ref = H5D.read(dset_id, 'H5ML_DEFAULT', dspace_id, dspace_id, plist);
  ref = H5D.read(dset_id,'H5T_STD_REF_OBJ',dspace_id,dspace_id,plist);
  for col = 1:size(ref,2)
    deref_dset_id = H5R.dereference(dset_id, 'H5R_OBJECT', ref(:,col));
    data = H5D.read(deref_dset_id, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');
    attr_id = H5A.open(deref_dset_id, 'MATLAB_class')
    matlab_class = H5A.read(attr_id).'
    class_id = H5T.get_class(H5D.get_type(deref_dset_id))
    switch(class_id)
      case H5ML.get_constant_value('H5T_INTEGER')
        fprintf('Integer\n');
      case H5ML.get_constant_value('H5T_FLOAT')
        fprintf('Floating point\n');
      case H5ML.get_constant_value('H5T_STRING')
        fprintf('String\n');
      case H5ML.get_constant_value('H5T_BITFIELD')
        fprintf('Bitfield\n');
      case H5ML.get_constant_value('H5T_OPAQUE')
        fprintf('Opaque\n');
      case H5ML.get_constant_value('H5T_COMPOUND')
        fprintf('Compound'\n');
      case H5ML.get_constant_value('H5T_REFERENCE')
        fprintf('Reference\n');
      case H5ML.get_constant_value('H5T_ENUM')
        fprintf('Enumerated\n');
      case H5ML.get_constant_value('H5T_VLEN')
        fprintf('Variable length\n');
      case H5ML.get_constant_value('H5T_ARRAY')
        fprintf('Array\n');
    end
    
    char(data)
    H5D.close(deref_dset_id);
  end
  
else
  data = H5D.read(dset_id, 'H5ML_DEFAULT', dspace_id, dspace_id, 'H5P_DEFAULT');
  char(data)
end

class_id = H5T.get_class(type_id);
switch(class_id)
  case H5ML.get_constant_value('H5T_INTEGER')
    fprintf('Integer\n');
  case H5ML.get_constant_value('H5T_FLOAT')
    fprintf('Floating point\n');
  case H5ML.get_constant_value('H5T_STRING')
    fprintf('String\n');
  case H5ML.get_constant_value('H5T_BITFIELD')
    fprintf('Bitfield\n');
  case H5ML.get_constant_value('H5T_OPAQUE')
    fprintf('Opaque\n');
  case H5ML.get_constant_value('H5T_COMPOUND')
    fprintf('Compound'\n');
  case H5ML.get_constant_value('H5T_REFERENCE')
    fprintf('Reference\n');
  case H5ML.get_constant_value('H5T_ENUM')
    fprintf('Enumerated\n');
  case H5ML.get_constant_value('H5T_VLEN')
    fprintf('Variable length\n');
  case H5ML.get_constant_value('H5T_ARRAY')
    fprintf('Array\n');
end

H5T.close(type_id);
H5S.close(dspace_id);
H5D.close(dset_id);
H5F.close(fid);

% I have a struct array created by matlab and stored in v7.3 format mat file:
%
% struArray = struct('name', {'one', 'two', 'three'},
%                    'id', {1,2,3},
%                    'data', {[1:10], [3:9], [0]})
% save('test.mat', 'struArray', '-v7.3')
%
% >>> import h5py
% >>> f = h5py.File('test.mat')
% >>> list(f.keys())
% ['#refs#', 'struArray']
% >>> struArray = f['struArray']
% >>> struArray['name'][0, 0]  # this is the HDF5 reference
% <HDF5 object reference>
% >>> f[struArray['name'][0, 0]].value  # this is the actual data
% array([[111],
%        [110],
%        [101]], dtype=uint16)
% To read struArray(i).id:
%
% >>> f[struArray['id'][0, 0]][0, 0]
% 1.0
% >>> f[struArray['id'][1, 0]][0, 0]
% 2.0
% >>> f[struArray['id'][2, 0]][0, 0]
% 3.0
% Notice that Matlab stores a number as an array of size (1, 1), hence the final [0, 0] to get the number.
%
% To read struArray(i).data:
%
% >>> f[struArray['data'][0, 0]].value
% array([[  1.],
%        [  2.],
%        [  3.],
%        [  4.],
%        [  5.],
%        [  6.],
%        [  7.],
%        [  8.],
%        [  9.],
%        [ 10.]])
% To read struArray(i).name, it is necessary to convert the array of integers to string:
%
% >>> f[struArray['name'][0, 0]].value.tobytes()[::2].decode()
% 'one'
% >>> f[struArray['name'][1, 0]].value.tobytes()[::2].decode()
% 'two'
% >>> f[struArray['name'][2, 0]].value.tobytes()[::2].decode()
% 'three'










return
%%
fn = 'test.mat';

% Get the list of all layer names from .mat -7.3 file without having to
% load the file
dataset = '/layer/name';
fid = H5F.open(fn, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
dset_id = H5D.open(fid, dataset);
if (H5ML.compare_values(H5T.get_class(H5D.get_type(dset_id)), ...
    H5ML.get_constant_value('H5T_REFERENCE')))
  ref = H5D.read(dset_id,'H5T_STD_REF_OBJ','H5S_ALL','H5S_ALL','H5P_DEFAULT');
  layer_names = cell(1,size(ref,2));
  for col = 1:size(ref,2)
    deref_dset_id = H5R.dereference(dset_id, 'H5R_OBJECT', ref(:,col));
    data = H5D.read(deref_dset_id, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');
    layer_names{col} = char(data);
    H5D.close(deref_dset_id);
  end
else
  data = H5D.read(dset_id, 'H5ML_DEFAULT', dspace_id, dspace_id, 'H5P_DEFAULT');
  layer_names = {char(data)};
end
H5D.close(dset_id);
H5F.close(fid);

layer_names
