function nums = arena_convert_range(str)
% Convert 0:2,4,5,7:10 to [0 1 2 3 4 7 8 9 10]

nums = eval(['[' str(:).' ']']);
