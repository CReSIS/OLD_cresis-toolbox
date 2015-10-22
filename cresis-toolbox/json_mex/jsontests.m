% MATLAB-JSON TESTS
% MATLAB-JSON FROM https://github.com/christianpanton/matlab-json/blob/master/bin/windows_precompiled.zip 3/3/2014

% cd 'C:\Users\Kyle\Documents\GitHub\matlab-json\bin\windows_precompiled';

fprintf('you should clear the memory before running these tests...\n ''dbcont'' to automatically do it ''dbquit'' to do it yourself.\n');
keyboard;
clear all

float_to_json_test = 12345.6789012345;
if ~strcmp(tojson(float_to_json_test),'12345.6789012345')
    warning('float to json failed');
end

json_to_float_test = '{"float_test":12345.6789012345}';
tmp = fromjson(json_to_float_test);
if ~(tmp.float_test == 12345.6789012345)
    warning('json to float failed');
end

string_to_json_test = 'thisIsATest';
if ~strcmp(tojson(string_to_json_test),'"thisIsATest"')
    warning('string to json failed');
end

json_to_string_test = '{"string_test":"thisIsATest"}';
tmp = fromjson(json_to_string_test);
if ~strcmp(tmp.string_test,'thisIsATest')
    warning('json to string failed');
end

cell_to_json_test = {'thisIsATest','thisIsATest'};
if ~isequal(tojson(cell_to_json_test),'[ "thisIsATest", "thisIsATest" ]')
    warning('cell to json failed');
end

json_to_cell_test = '{"cell_test":[ "thisIsATest", "thisIsATest" ]}';
tmp = fromjson(json_to_cell_test);
if ~isequal(tmp.cell_test,{'thisIsATest';'thisIsATest'})
    warning('json to cell failed');
end

array_to_json_test = [12345.6789012345,12345.6789012345,12345.6789012345];
if ~strcmp(tojson(array_to_json_test),'[ 12345.6789012345, 12345.6789012345, 12345.6789012345 ]')
    warning('array to json failed');
end

json_to_array_test = '{"array_test":[ 12345.6789012345, 12345.6789012345, 12345.6789012345 ]}';
tmp = fromjson(json_to_array_test);
if ~isequal(cell2mat(tmp.array_test)',[12345.6789012345,12345.6789012345,12345.6789012345])
    warning('json to array failed');
end

struct_to_json_test.sfloat = float_to_json_test;
struct_to_json_test.sstring = string_to_json_test;
struct_to_json_test.scell = cell_to_json_test;
struct_to_json_test.sarray = array_to_json_test;
if ~strcmp(tojson(struct_to_json_test),'{ "sfloat": 12345.6789012345, "sstring": "thisIsATest", "scell": [ "thisIsATest", "thisIsATest" ], "sarray": [ 12345.6789012345, 12345.6789012345, 12345.6789012345 ] }')
    warning('struct to json failed');
end

json_to_struct_test = '{ "sfloat": 12345.6789012345, "sstring": "thisIsATest", "scell": [ "thisIsATest", "thisIsATest" ], "sarray": [ 12345.6789012345, 12345.6789012345, 12345.6789012345 ] }';
tmp = fromjson(json_to_struct_test);
tmp.scell = tmp.scell';
tmp.sarray = cell2mat(tmp.sarray)';
if ~isequal(tmp,struct_to_json_test)
    warning('json to struct failed');
end

% MEMORY LEAK TESTING (AS IS THIS WILL CRASH MATLAB)
clear all;
[usr, sys] = memory;
preTestMem = usr.MemUsedMATLAB;
data = ones(10000000,1) * 12345.6789012345; % SHOULD BE ABLE TO DO THIS...
% data = ones(10000000,1);
tmp = tojson(data);
clear tmp data;
[usr, sys] = memory;
postTestMem = usr.MemUsedMATLAB;
fprintf('Before %f GB and After %f GB\n',preTestMem/1000000000,postTestMem/1000000000);

