function TDMS_TODO_List
%
%   THINGS TO MAYBE HANDLE  ... (LOW PRIORITY)
%   ===============================================
%   - change variable names in code to allow for better self-documentation
%   
%   - allow merged reading of variables when decimated for faster read times
%       i.e if 5 channels are all doubles, read them all in at once and
%       then distribute to each preallocated array
%
%   - Reading subset of data for interleaved data
%   
%   - The ability to return all data casted as doubles, instead of their
%   original data class format (might be useful as calculations with
%   classes can be odd -> 3/4 => 1 for integers)
%
%   THINGS TO HANDLE ????   ... (VERY LOW PRIORITY)
%   =======================================================================
%   - Handle endianess, this should only involve a change in fread ...
%
%   - Handle corrupted end of file writing (allow soft error) (there is a
%   flag that can be set if this happens)
%
%   - add on waveform handling as well (look for properties and change data type)
%
%   - optional property casting
